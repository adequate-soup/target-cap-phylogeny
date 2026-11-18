#!/bin/bash

# NOTE: the slurm_output/ directory must exist BEFORE you run sbatch,
# otherwise Slurm cannot open the log files and the job fails silently:
#   mkdir -p slurm_output && sbatch bbnorm_50x.sh -i /path/to/reads

set -euo pipefail

usage() {
    cat >&2 <<EOF
Usage: sbatch $0 -i <reads_dir> [-c <target_coverage>] [-td <trimmed_dir>] [-nd <normalized_dir>]
  -i,  --reads-dir         Directory with *.fastq.gz / *.fq.gz reads (required)
  -s,  --skip-norm         Skips BBNorm
  -c,  --target-coverage   BBNorm target depth (default: 50)
  -td, --trimmed-dir       Output dir for fastp-trimmed reads (default: <reads_dir>_trimmed)
  -nd, --normalized-dir    Output dir for normalized reads   (default: <reads_dir>_normalized)
EOF
    exit 1
}

READS_DIR=""
TARGET_COVERAGE=50
OUT_DIR_TRIM=""
OUT_DIR_NORM=""
SKIP_NORM=0

[[ $# -lt 1 ]] && usage

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--reads-dir|-c|--target-coverage|-td|--trimmed-dir|-nd|--normalized-dir)
            [[ $# -lt 2 ]] && { echo "Error: $1 requires a value" >&2; usage; } ;;&
        -s|--skip-norm) SKIP_NORM=1; shift 1 ;;
        -i|--reads-dir) READS_DIR="$2"; shift 2 ;;
        -c|--target-coverage) TARGET_COVERAGE="$2"; shift 2 ;;
        -td|--trimmed-dir) OUT_DIR_TRIM="$2"; shift 2 ;;
        -nd|--normalized-dir) OUT_DIR_NORM="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

if [[ -z "$READS_DIR" || ! -d "$READS_DIR" ]]; then
    echo "Error: reads directory '${READS_DIR}' not provided or does not exist" >&2
    exit 1
fi

module load fastp

# Absolute path, no trailing slash (otherwise "dir/_trimmed" would be created inside dir)
READS_DIR="$(realpath "$READS_DIR")"

OUT_DIR_TRIM="${OUT_DIR_TRIM:-${READS_DIR}_trimmed}"
OUT_DIR_NORM="${OUT_DIR_NORM:-${READS_DIR}_normalized}"

mkdir -p "$OUT_DIR_TRIM"
OUT_DIR_TRIM="$(realpath "$OUT_DIR_TRIM")"

REPORT_DIR="${OUT_DIR_TRIM}/fastp_reports"
mkdir -p "$REPORT_DIR"

BIND_PATHS="${READS_DIR},${OUT_DIR_TRIM}"
BBTOOLS_SIF=""
if [[ "$SKIP_NORM" = 0 ]]; then
    module load apptainer 2>/dev/null || true

    SIF_DIR="${SIF_DIR:-$HOME/apptainer_images}"
    mkdir -p "$SIF_DIR"
    BBTOOLS_SIF="${SIF_DIR}/bbtools_latest.sif"

    if [[ ! -f "$BBTOOLS_SIF" ]]; then
        echo "[Init] Pulling BBTools image to $BBTOOLS_SIF ..."
        apptainer pull "$BBTOOLS_SIF" "docker://staphb/bbtools:latest"
    fi

    mkdir -p "$OUT_DIR_NORM"
    OUT_DIR_NORM="$(realpath "$OUT_DIR_NORM")"
    BIND_PATHS+=",${OUT_DIR_NORM}"
fi
export APPTAINER_BINDPATH="$BIND_PATHS"
export SINGULARITY_BINDPATH="$BIND_PATHS"

# Detect total RAM (GB) and CPU cores available on node
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    TOTAL_RAM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))   # SLURM_MEM_PER_NODE is in MB
else
    TOTAL_RAM_GB=$(awk '/MemTotal/ {print int($2/1024/1024)}' /proc/meminfo)
fi

TOTAL_CORES="${SLURM_CPUS_ON_NODE:-$(nproc)}"
THREADS_PER_TASK=8
PARALLEL_JOBS=$(( TOTAL_CORES / THREADS_PER_TASK ))
[[ "$PARALLEL_JOBS" -lt 1 ]] && PARALLEL_JOBS=1

# RAM limit per job with a 5% safety buffer for system overhead
MEM_PER_JOB_GB=$(( (TOTAL_RAM_GB * 95 / 100) / PARALLEL_JOBS ))
[[ "$MEM_PER_JOB_GB" -lt 4 ]] && MEM_PER_JOB_GB=4

echo "===================================================="
echo "Input Directory:       $READS_DIR"
echo "Trim Output Directory: $OUT_DIR_TRIM"
echo "Norm Output Directory: $OUT_DIR_NORM"
echo "Target Coverage:       $TARGET_COVERAGE"
echo "BBTools Image:         $BBTOOLS_SIF"
echo "Total Node RAM:        ${TOTAL_RAM_GB} GB"
echo "Parallel Workers:      $PARALLEL_JOBS concurrent jobs ($THREADS_PER_TASK threads/job)"
echo "Heap Per Worker:       ${MEM_PER_JOB_GB} GB (-Xmx${MEM_PER_JOB_GB}g)"
echo "===================================================="

export READS_DIR OUT_DIR_TRIM OUT_DIR_NORM REPORT_DIR THREADS_PER_TASK \
       BBTOOLS_SIF MEM_PER_JOB_GB TARGET_COVERAGE

# Outputs are written to a temp name and renamed on success, so an interrupted
# run never leaves a truncated file that the skip-check mistakes for a finished one.

trim_sample() {
    local first_file="$1"
    local base_fname
    base_fname="$(basename "$first_file")"

    if [[ "$base_fname" =~ _R1 ]]; then
        local r2_name="${base_fname/_R1/_R2}"
        local r1="${READS_DIR}/${base_fname}"
        local r2="${READS_DIR}/${r2_name}"

        if [[ ! -f "$r2" ]]; then
            echo "Error: Paired file '$r2' not found for '$r1'. Skipping..." >&2
            return 1
        fi

        local out_r1="${OUT_DIR_TRIM}/${base_fname}"
        local out_r2="${OUT_DIR_TRIM}/${r2_name}"

        if [[ -f "$out_r1" && -f "$out_r2" ]]; then
            echo "[Skip] Trimmed outputs already exist for $base_fname"
            return 0
        fi

        echo "Trimming ${base_fname} and ${r2_name}..."

        local tmp_r1="${OUT_DIR_TRIM}/.tmp.${base_fname}"
        local tmp_r2="${OUT_DIR_TRIM}/.tmp.${r2_name}"
        local prefix="${REPORT_DIR}/${base_fname%%_R1*}"

        fastp \
            -i "$r1" -I "$r2" \
            -o "$tmp_r1" -O "$tmp_r2" \
            --detect_adapter_for_pe \
            --correction \
            -l 36 \
            --thread "$THREADS_PER_TASK" \
            -h "${prefix}_fastp.html" -j "${prefix}_fastp.json"

        mv "$tmp_r1" "$out_r1"
        mv "$tmp_r2" "$out_r2"
    else
        local unpaired="${READS_DIR}/${base_fname}"
        local out_unpaired="${OUT_DIR_TRIM}/${base_fname}"

        if [[ -f "$out_unpaired" ]]; then
            echo "[Skip] Trimmed output already exists for $base_fname"
            return 0
        fi

        echo "Trimming $base_fname..."

        local tmp_unpaired="${OUT_DIR_TRIM}/.tmp.${base_fname}"
        local prefix="${REPORT_DIR}/${base_fname%.gz}"

        fastp \
            -i "$unpaired" \
            -o "$tmp_unpaired" \
            -l 36 \
            --thread "$THREADS_PER_TASK" \
            -h "${prefix}_fastp.html" -j "${prefix}_fastp.json"

        mv "$tmp_unpaired" "$out_unpaired"
    fi
}
export -f trim_sample

run_bbnorm() {
    local first_file="$1"
    local base_fname
    base_fname="$(basename "$first_file")"
    if [[ "$base_fname" =~ _R1 ]]; then
        local r2_name="${base_fname/_R1/_R2}"
        local r1="${OUT_DIR_TRIM}/${base_fname}"
        local r2="${OUT_DIR_TRIM}/${r2_name}"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            echo "Error: trimmed pair '$r1' / '$r2' not found. Skipping..." >&2
            return 1
        fi

        local out_r1="${OUT_DIR_NORM}/${base_fname}"
        local out_r2="${OUT_DIR_NORM}/${r2_name}"

        if [[ -f "$out_r1" && -f "$out_r2" ]]; then
            echo "[Skip] Normalized outputs already exist for $base_fname"
            return 0
        fi

        echo "Normalizing ${base_fname} and ${r2_name} to $TARGET_COVERAGE depth..."

        local tmp_r1="${OUT_DIR_NORM}/.tmp.${base_fname}"
        local tmp_r2="${OUT_DIR_NORM}/.tmp.${r2_name}"

        apptainer exec "$BBTOOLS_SIF" bbnorm.sh \
            "-Xmx${MEM_PER_JOB_GB}g" \
            in1="$r1" \
            in2="$r2" \
            out1="$tmp_r1" \
            out2="$tmp_r2" \
            target="$TARGET_COVERAGE" \
            min=2 \
            threads="$THREADS_PER_TASK"

        mv "$tmp_r1" "$out_r1"
        mv "$tmp_r2" "$out_r2"
    else
        local unpaired="${OUT_DIR_TRIM}/${base_fname}"
        local out_unpaired="${OUT_DIR_NORM}/${base_fname}"

        if [[ ! -f "$unpaired" ]]; then
            echo "Error: trimmed file '$unpaired' not found. Skipping..." >&2
            return 1
        fi

        if [[ -f "$out_unpaired" ]]; then
            echo "[Skip] Normalized output already exists for $base_fname"
            return 0
        fi

        echo "Normalizing $base_fname to $TARGET_COVERAGE depth..."

        local tmp_unpaired="${OUT_DIR_NORM}/.tmp.${base_fname}"

        apptainer exec "$BBTOOLS_SIF" bbnorm.sh \
            "-Xmx${MEM_PER_JOB_GB}g" \
            in="$unpaired" \
            out="$tmp_unpaired" \
            target="$TARGET_COVERAGE" \
            min=2 \
            threads="$THREADS_PER_TASK"

        mv "$tmp_unpaired" "$out_unpaired"
    fi
}
export -f run_bbnorm

mapfile -t SAMPLE_FILES < <(find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) ! -name "*_R2*" | sort)

if [[ ${#SAMPLE_FILES[@]} -eq 0 ]]; then
    echo "No matching FASTQ files found in $READS_DIR" >&2
    exit 1
fi

echo "Found ${#SAMPLE_FILES[@]} sample file sets to process."

parallel --line-buffer -j "$PARALLEL_JOBS" trim_sample ::: "${SAMPLE_FILES[@]}"

if [[ "$SKIP_NORM" = 1 ]]; then
    echo "Trimming complete (BBNorm skipped)."
    exit 0
fi

parallel --line-buffer -j "$PARALLEL_JOBS" run_bbnorm ::: "${SAMPLE_FILES[@]}"

echo "Trimming and BBnorm downsampling pipeline completed successfully."
