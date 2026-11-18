#!/bin/bash

set -euo pipefail

IN_DIR=""
OUT_DIR=""
DS_FRAC=0.6

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--in-dir) IN_DIR="${2:?Missing value for $1}"; shift 2 ;;
        -o|--out-dir) OUT_DIR="${2:?Missing value for $1}"; shift 2 ;;
        -f|--downsample-fraction) DS_FRAC="${2:?Missing value for $1}"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$IN_DIR" || ! -d "$IN_DIR" ]]; then
    echo "Error: reads directory '${IN_DIR}' not provided or does not exist" >&2
    exit 1
fi

if ! awk -v f="$DS_FRAC" 'BEGIN { exit !(f > 0 && f <= 1) }'; then
    echo "Error: downsample fraction must be in (0, 1], got '${DS_FRAC}'" >&2
    exit 1
fi

IN_DIR="${IN_DIR%/}"
OUT_DIR="${OUT_DIR:-${IN_DIR}_75bp_SE}"
mkdir -p "$OUT_DIR"
OUT_DIR="$(realpath "$OUT_DIR")"

# Load modules on Savio (names may vary; check `module avail`)
module load seqtk 2>/dev/null || true
for cmd in seqtk pigz parallel; do
    command -v "$cmd" >/dev/null || { echo "Error: '$cmd' not found in PATH" >&2; exit 1; }
done

# Detect all allocated CPU cores on the exclusive node
TOTAL_CPUS=${SLURM_CPUS_ON_NODE:-$(nproc)}

COMPRESS_THREADS=4
JOBS=$((TOTAL_CPUS / COMPRESS_THREADS))
[[ "$JOBS" -lt 1 ]] && JOBS=1

echo "Running on $(hostname) with $TOTAL_CPUS cores ($JOBS concurrent samples, $COMPRESS_THREADS compression threads each)."

export IN_DIR OUT_DIR COMPRESS_THREADS DS_FRAC

process_sample() {
    set -euo pipefail   # not inherited by parallel's subshell

    local r1="$1"
    local sample
    sample=$(basename "$r1" _R1_trimmed.fastq.gz)
    local r2="${IN_DIR}/${sample}_R2_trimmed.fastq.gz"
    local out="${OUT_DIR}/${sample}_trimmed.fastq.gz"
    local tmp="${out}.tmp"

    if [[ ! -f "$r2" ]]; then
        echo "Skipping ${sample}: missing R2 file ($r2)" >&2
        return 0
    fi

    # Remove stale or corrupted tiny files (< 100 bytes) from previous failed runs
    if [[ -f "$out" ]] && [[ ! -s "$out" || $(stat -c%s "$out" 2>/dev/null || stat -f%z "$out" 2>/dev/null) -lt 100 ]]; then
        rm -f "$out"
    fi

    if [[ -s "$out" ]]; then
        echo "Skipping ${sample}: output already exists ($out)"
        return 0
    fi

    echo "Processing ${sample}..."
    trap 'rm -f "$tmp"' ERR
    {
        seqtk sample -s100 "$r1" "$DS_FRAC" | seqtk trimfq -q 0 -L 75 -
        seqtk sample -s101 "$r2" "$DS_FRAC" | seqtk trimfq -q 0 -L 75 -
    } | pigz -p "$COMPRESS_THREADS" > "$tmp"
    mv "$tmp" "$out"
}

export -f process_sample

mapfile -t R1_FILES < <(find "$IN_DIR" -maxdepth 1 -name "*_R1_trimmed.fastq.gz" | sort)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "No *_R1_trimmed.fastq.gz files found in $IN_DIR" >&2
    exit 1
fi

echo "Found ${#R1_FILES[@]} samples to process."
printf "%s\n" "${R1_FILES[@]}" | parallel --line-buffer --halt now,fail=1 -j "$JOBS" process_sample {}

echo "Done."
