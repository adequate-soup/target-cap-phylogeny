#!/bin/bash

set -euo pipefail

# Ensure slurm log folder exists
mkdir -p slurm_output

usage() {
    echo "Usage: sbatch $0 [sample_list.txt] -o <output_dir> -i <input_dir1> [input_dir2 ...]" >&2
    echo "  sample_list.txt  Optional. If omitted, samples are inferred from filenames in the input dirs" >&2
    echo "  -o               Output directory for concatenated reads" >&2
    echo "  -i               Space-separated list of input directories (must be the last flag)" >&2
    exit 1
}

# ---- Argument parsing ----
SAMPLE_LIST=""
OUTPUT_DIR=""
INPUT_DIRS=()

# Optional sample list: first argument, if it is not a flag
if [[ $# -ge 1 && "$1" != -* ]]; then
    SAMPLE_LIST="$1"
    shift
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o)
            [[ $# -ge 2 ]] || { echo "Error: -o requires an argument." >&2; usage; }
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -i)
            shift
            # Consume everything remaining that is not a flag
            while [[ $# -gt 0 && "$1" != -* ]]; do
                INPUT_DIRS+=("$1")
                shift
            done
            ;;
        *)
            echo "Error: Unknown argument '$1'." >&2
            usage
            ;;
    esac
done

if [[ -n "$SAMPLE_LIST" && ! -f "$SAMPLE_LIST" ]]; then
    echo "Error: Sample list file '$SAMPLE_LIST' does not exist." >&2
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: -o <output_dir> is required." >&2
    usage
fi

if [[ ${#INPUT_DIRS[@]} -eq 0 ]]; then
    echo "Error: -i requires at least one input directory." >&2
    usage
fi

mkdir -p "$OUTPUT_DIR"

# Join input dirs into a space-separated string (so it can be exported to parallel)
SOURCE_DIRS="${INPUT_DIRS[*]}"
CONCAT_DIR="$OUTPUT_DIR"

# Concatenation is I/O-bound, but use all cores for simplicity
PROCESSES=${SLURM_CPUS_ON_NODE:-$(nproc)}

# Minimum number of samples that must share a prefix before it is stripped
MIN_PREFIX_COUNT=3

module load gnu-parallel 2>/dev/null || true

# ---- Infer sample IDs from filenames (used only when no sample list is given) ----
# Step 1: take each R1 filename and strip:
#   - everything from "_R1" onward   (e.g. _R1_001.fastq.gz)
#   - a trailing lane tag            (e.g. _L001)
#   - a trailing sample-number tag   (e.g. _S12)
#   Example: Pap_Plant7_S12_L001_R1_001.fastq.gz -> Pap_Plant7
# Step 2: the prefix is the first underscore-delimited token (e.g. "Pap").
#   If at least MIN_PREFIX_COUNT samples share it, strip "<prefix>_" from the
#   OUTPUT name (Pap_Plant7 -> Plant7). The original ID is kept for matching input files.
#   If stripping would make two outputs collide, those samples keep their full name.
# Output: one line per sample, "<original_id><TAB><output_id>"
infer_samples() {
    local dir
    for dir in $SOURCE_DIRS; do
        if [[ -d "$dir" ]]; then
            find "$dir" -maxdepth 1 -name "*_R1*.fastq.gz" -printf '%f\n'
        else
            echo "Warning: input directory '$dir' does not exist." >&2
        fi
    done \
        | sed -E 's/_R1.*$//; s/_L[0-9]+$//; s/_S[0-9]+$//' \
        | sort -u \
        | awk -v min="$MIN_PREFIX_COUNT" '
            {
                ids[NR] = $0
                n = index($0, "_")
                if (n > 1) { pre[NR] = substr($0, 1, n - 1); cnt[pre[NR]]++ }
            }
            END {
                for (i = 1; i <= NR; i++) {
                    short = ids[i]
                    if (pre[i] != "" && cnt[pre[i]] >= min) {
                        short = substr(ids[i], length(pre[i]) + 2)
                        if (short == "") short = ids[i]
                    }
                    out[i] = short
                    seen[short]++
                }
                for (i = 1; i <= NR; i++) {
                    name = out[i]
                    if (seen[name] > 1) name = ids[i]   # avoid output collisions
                    print ids[i] "\t" name
                }
            }'
}

if [[ -n "$SAMPLE_LIST" ]]; then
    AUTO_MODE=0
    # List mode: output name is the ID exactly as given
    SAMPLES=$(grep -v '^\s*$' "$SAMPLE_LIST" | grep -v '^\s*#' | awk '{print $0 "\t" $0}' || true)
else
    AUTO_MODE=1
    SAMPLES=$(infer_samples)
fi

if [[ -z "$SAMPLES" ]]; then
    echo "Error: no samples found." >&2
    exit 1
fi

echo "Running on $(hostname) with $PROCESSES parallel concatenation jobs."
echo "Searching directories: $SOURCE_DIRS"
echo "Output directory: $CONCAT_DIR"
if [[ "$AUTO_MODE" -eq 1 ]]; then
    echo "No sample list given; inferred $(echo "$SAMPLES" | wc -l) sample(s) from filenames (matched ID -> output name):"
    echo "$SAMPLES" | awk -F'\t' '{print "  " $1 " -> " $2}'
fi

export SOURCE_DIRS CONCAT_DIR AUTO_MODE

# $1 = sample ID used to match input files, $2 = name used for output files
process_sample() {
    local sample_id="$1"
    local out_id="$2"

    local concat_r1="${CONCAT_DIR}/${out_id}_R1.fastq.gz"
    local concat_r2="${CONCAT_DIR}/${out_id}_R2.fastq.gz"

    # Skip if concatenated outputs already exist
    if [[ -f "$concat_r1" && -f "$concat_r2" ]]; then
        echo "Skipping sample ${out_id} - concatenated outputs already exist."
        return 0
    fi

    # With a user-supplied list, match the ID anywhere in the filename (original behavior).
    # With inferred IDs, anchor at the start of the filename and require "_" after the ID,
    # so "Sample1" does not also match "Sample10".
    local r1_pat r2_pat
    if [[ "$AUTO_MODE" -eq 1 ]]; then
        r1_pat="${sample_id}_*R1*.fastq.gz"
        r2_pat="${sample_id}_*R2*.fastq.gz"
    else
        r1_pat="*${sample_id}*_R1*.fastq.gz"
        r2_pat="*${sample_id}*_R2*.fastq.gz"
    fi

    # Collect matching R1 and R2 files across all input directories
    local r1_files=()
    local r2_files=()

    for dir in $SOURCE_DIRS; do
        if [[ -d "$dir" ]]; then
            while IFS= read -r -d '' f; do
                r1_files+=("$f")
            done < <(find "$dir" -maxdepth 1 -name "$r1_pat" -print0 2>/dev/null | sort -z)

            while IFS= read -r -d '' f; do
                r2_files+=("$f")
            done < <(find "$dir" -maxdepth 1 -name "$r2_pat" -print0 2>/dev/null | sort -z)
        fi
    done

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "Warning: No R1 files found for sample $sample_id in $SOURCE_DIRS." >&2
        return 0
    fi

    if [[ ${#r1_files[@]} -ne ${#r2_files[@]} ]]; then
        echo "Warning: Sample $sample_id has ${#r1_files[@]} R1 file(s) but ${#r2_files[@]} R2 file(s); skipping." >&2
        return 0
    fi

    # Write to temp files then rename, so an interrupted job never leaves a partial output
    if [[ ! -f "$concat_r1" ]]; then
        echo "Concatenating ${#r1_files[@]} R1 file(s) for sample: $out_id"
        cat "${r1_files[@]}" > "${concat_r1}.tmp"
        mv "${concat_r1}.tmp" "$concat_r1"
    fi

    if [[ ! -f "$concat_r2" ]]; then
        echo "Concatenating ${#r2_files[@]} R2 file(s) for sample: $out_id"
        cat "${r2_files[@]}" > "${concat_r2}.tmp"
        mv "${concat_r2}.tmp" "$concat_r2"
    fi
}
export -f process_sample

echo "$SAMPLES" | parallel -j "$PROCESSES" --colsep '\t' process_sample {1} {2}