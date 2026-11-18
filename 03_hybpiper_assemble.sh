#!/bin/bash

set -euo pipefail

mkdir -p slurm_output hybpiper_out

# 1. Environment & Paths
export PATH="$HOME/.conda/envs/hybpiper/bin:${PATH}"
export SHELL=/bin/bash   # GNU parallel needs bash to see exported functions

WORK_DIR="$(realpath hybpiper_out)"   # default output dir; override with -o
THREADS_PER_JOB=8
READS_DIR=""
TARGET_FILE=""
NO_CLOBBER=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--reads-dir|-t|--target-file|-o|--output-dir|-tj|--threads-per-job)
            [[ $# -lt 2 ]] && { echo "Error: $1 requires a value" >&2; exit 1; } ;;&
        -i|--reads-dir) READS_DIR="$2"; shift 2 ;;
        -t|--target-file) TARGET_FILE="$2"; shift 2 ;;
        -o|--output-dir) WORK_DIR="$2"; shift 2 ;;
        -tj|--threads-per-job) THREADS_PER_JOB="$2"; shift 2 ;;
        --no-clobber) NO_CLOBBER=1; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$READS_DIR" || ! -d "$READS_DIR" ]]; then
    echo "Error: reads directory '${READS_DIR}' not provided or does not exist" >&2
    exit 1
fi

if [[ -z "$TARGET_FILE" || ! -f "$TARGET_FILE" ]]; then
    echo "Error: target file '${TARGET_FILE}' not provided or not found." >&2
    exit 1
fi

READS_DIR="$(realpath "$READS_DIR")"
TARGET_FILE="$(realpath "$TARGET_FILE")"
mkdir -p "$WORK_DIR"
WORK_DIR="$(realpath "$WORK_DIR")"

# 2. CPU / Memory / Process Calculation
TOTAL_RAM_GB=0
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    TOTAL_RAM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
fi
# --mem=0 can leave SLURM_MEM_PER_NODE unset or 0, so fall back to physical RAM
if [[ "$TOTAL_RAM_GB" -lt 1 ]]; then
    TOTAL_RAM_GB=$(awk '/MemTotal/ {print int($2/1024/1024)}' /proc/meminfo)
fi
[[ "$TOTAL_RAM_GB" -lt 1 ]] && TOTAL_RAM_GB=1

CPUS_PER_NODE=${SLURM_CPUS_ON_NODE:-$(nproc)}
PROCESSES_PER_NODE=$(( CPUS_PER_NODE / THREADS_PER_JOB ))
[[ "$PROCESSES_PER_NODE" -lt 1 ]] && PROCESSES_PER_NODE=1
MEMORY_PER_JOB=$(( TOTAL_RAM_GB / PROCESSES_PER_NODE ))
[[ "$MEMORY_PER_JOB" -lt 1 ]] && MEMORY_PER_JOB=1   # informational only; no ulimit is applied

echo "===================================================="
echo "Node RAM:             ${TOTAL_RAM_GB} GB Total"
echo "CPUs Per Node:        ${CPUS_PER_NODE} CPUs"
echo "Workers Running:      ${PROCESSES_PER_NODE} concurrent processes"
echo "Per-Worker Alloc:     ${THREADS_PER_JOB} threads | ~${MEMORY_PER_JOB} GB RAM share (not enforced)"
echo "===================================================="

# 3. Helper Functions (exported so GNU parallel subshells can use them)

# Strip directory, .gz/.fastq/.fq, then the last _R1/_R2 token and anything after it.
# Sample_A_R1_001.fastq.gz -> Sample_A ; Sample_A.fq.gz -> Sample_A
get_sample_id() {
    local name id prev="" s
    name="$(basename "$1")"
    name="${name%.gz}"; name="${name%.fastq}"; name="${name%.fq}"
    if [[ "$name" =~ ^(.*)_R[12]([._].*)?$ ]]; then
        id="${BASH_REMATCH[1]}"
    else
        id="$name"
    fi
    # Strip trailing processing suffixes repeatedly until nothing changes
    while [[ "$id" != "$prev" ]]; do
        prev="$id"
        for s in _trimmed _trim _normalized _norm; do
            id="${id%"$s"*}"
        done
    done
    echo "$id"
}

# True if the file's *basename* is an R2 file (directory names are ignored)
is_r2() {
    [[ "$(basename "$1")" =~ _R2([._]|$) ]]
}

# Print the R2 path for an R1 file by editing only the basename; fail if not an R1 file
get_r2_path() {
    local dir base
    dir="$(dirname "$1")"
    base="$(basename "$1")"
    if [[ "$base" =~ ^(.*)_R1([._].*)$ ]]; then
        echo "${dir}/${BASH_REMATCH[1]}_R2${BASH_REMATCH[2]}"
    else
        return 1
    fi
}

# Completion is judged by the presence of HybPiper's end-of-assemble output file
sample_is_complete() {
    local sample_id="$1"
    [[ -s "${WORK_DIR}/${sample_id}/${sample_id}_genes_with_seqs.txt" ]]
}

# True if the sample already has a supercontig FASTA in the output directory
# (also checks the sample's own subdirectory in case HybPiper wrote it there)
has_supercontig() {
    local sample_id="$1"
    [[ -s "${WORK_DIR}/${sample_id}_supercontig.fasta" || \
       -s "${WORK_DIR}/${sample_id}/${sample_id}_supercontig.fasta" ]]
}

# 4. Worker Function
assemble() {
    local input_reads="$1"
    local sample_id
    sample_id="$(get_sample_id "$input_reads")"

    local hybpiper_reads=()
    local r2

    if r2="$(get_r2_path "$input_reads")"; then
        if [[ ! -f "$r2" ]]; then
            echo "Error: Matching R2 file '$r2' not found for $sample_id" >&2
            return 0
        fi
        hybpiper_reads=("$input_reads" "$r2")
        echo ">> [Worker $(hostname)] Assembling sample (Paired-End): $sample_id"
    else
        hybpiper_reads=("$input_reads")
        echo ">> [Worker $(hostname)] Assembling sample (Unpaired): $sample_id"
    fi

    # Inside an `if !` condition set -e is ignored, so each step needs an explicit || exit 1
    if ! (
        cd "$WORK_DIR" || exit 1

        hybpiper assemble \
            -r "${hybpiper_reads[@]}" \
            -t_dna "$TARGET_FILE" \
            --prefix "$sample_id" \
            --cpu "$THREADS_PER_JOB" \
            --bwa \
            --force_overwrite || exit 1

        hybpiper retrieve_sequences supercontig \
            -t_dna "$TARGET_FILE" \
            --single_sample_name "$sample_id" \
            --cpu "$THREADS_PER_JOB" || exit 1
    ); then
        echo ">> [Worker $(hostname)] FAILED sample $sample_id. Skipping..." >&2
        return 0
    fi
}

# 5. Gather Primary Reads Files (skip R2 files)
ALL_PRIMARY_FILES=()
while IFS= read -r -d '' file; do
    if is_r2 "$file"; then
        continue
    fi
    ALL_PRIMARY_FILES+=("$file")
done < <(find "$READS_DIR" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -print0 | sort -z)

if [[ ${#ALL_PRIMARY_FILES[@]} -eq 0 ]]; then
    echo "No matching FASTQ files found in '$READS_DIR'." >&2
    exit 1
fi

# 6. Check Unfinished Samples & Calculate Total File Size
check_unfinished_with_size() {
    local input_reads="$1"
    local sample_id
    sample_id="$(get_sample_id "$input_reads")"

    # --no-clobber: skip samples that already have a supercontig FASTA
    if [[ "${NO_CLOBBER:-0}" -eq 1 ]] && has_supercontig "$sample_id"; then
        return 0
    fi

    if ! sample_is_complete "$sample_id"; then
        local size1 size2 total_bytes r2
        size1=$(stat -c%s "$input_reads" 2>/dev/null || echo 0)
        total_bytes=$size1

        if r2="$(get_r2_path "$input_reads")" && [[ -f "$r2" ]]; then
            size2=$(stat -c%s "$r2" 2>/dev/null || echo 0)
            total_bytes=$(( size1 + size2 ))
        fi

        # Output format: SIZE_IN_BYTES <TAB> FILE_PATH
        printf "%d\t%s\n" "$total_bytes" "$input_reads"
    fi
}
export -f get_sample_id is_r2 get_r2_path sample_is_complete has_supercontig check_unfinished_with_size
export WORK_DIR NO_CLOBBER

# 7. Parallelize Check, Sort Ascending by Size, and Extract File Paths
readarray -t UNFINISHED_FILES < <(
    printf "%s\n" "${ALL_PRIMARY_FILES[@]}" | \
        parallel -j "$PROCESSES_PER_NODE" check_unfinished_with_size {} | \
        sort -n -k1,1 | \
        cut -f2-
)

if [[ ${#UNFINISHED_FILES[@]} -eq 0 ]]; then
    echo "All ${#ALL_PRIMARY_FILES[@]} samples completed successfully!"
    exit 0
fi

echo "Found ${#ALL_PRIMARY_FILES[@]} total samples. ${#UNFINISHED_FILES[@]} remaining to assemble (ordered smallest to largest)."

# Export symbols needed for subshells
export -f assemble
export TARGET_FILE THREADS_PER_JOB WORK_DIR

# 8. Parallel Execution on This Node
printf "%s\n" "${UNFINISHED_FILES[@]}" | \
    parallel --line-buffer -j "$PROCESSES_PER_NODE" assemble {}
