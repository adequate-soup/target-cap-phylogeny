#!/bin/bash

set -eo pipefail

# --------------------------------------------------------------------------
# 0. Input Validation & Argument Parsing
# --------------------------------------------------------------------------
# Usage: sbatch superexon_worker_v3.sh [options] [-e sample_to_exclude ...]
#   -gc  generate an ASTRAL constraint tree (constraint.tre)
#   -rr  two-pass mode: build unconstrained gene trees, make the constraint
#        tree (implies -gc), then resubmit this script with -cp constraint.tre.
#        The resubmitted pass writes to 5_gene_trees_constrained/ and builds
#        the final species tree from those trees.
#   -cp  path to a constraint tree (relative paths are taken from $WORKDIR)
#   -bt  bootstrap collapse threshold      -tp  TreeShrink quantile
#   -tt  TrimAl gap threshold              -e   sample names to exclude (LAST option)
#
# THIS VERSION NEVER DELETES FILES. It is resumable at two levels:
#   - step level: a .skip marker in $SKIP_DIR skips a whole finished step
#   - gene level: each per-gene job skips genes whose output already exists
# To force a step to rerun, move its output directory aside and delete its
# .skip marker by hand.
#
# INPUT: per-gene FASTAs from filter_low_copy.sh (<gene>.fasta, headers like
# ">SAMPLE_paralog_1"). Override the location with:
#   sbatch --export=ALL,GENE_FASTA_DIR=/path/to/filtered_genes superexon_worker_v3.sh
# Optional: MAPPING_FILE = two-column "sample species" file (see step 6).

export CPUS=${SLURM_CPUS_ON_NODE:-$(nproc)}
export PARALLEL_SHELL=/bin/bash   # exported bash functions need bash

export HYBPIPER_ENV="~/.conda/envs/hybpiper"
export ROOTDIR="~/"

# Environment overrides are honoured (previously these were overwritten)
GENE_FASTA_DIR="${GENE_FASTA_DIR:-$ROOTDIR/filtered_genes}"   # Output of filter_low_copy.sh
MAPPING_FILE="${MAPPING_FILE:-$ROOTDIR/data/mapping.txt}"

export WORKDIR="$ROOTDIR/analysis_pipeline"
export SKIP_DIR="$WORKDIR/skip"
export LOG_DIR="$WORKDIR/logs"

mkdir -p "$WORKDIR" "$SKIP_DIR" "$LOG_DIR"
cd "$WORKDIR"

# Timeout ceilings per step (seconds)
export TIMEOUT_MAFFT=1800      # 30 min
export TIMEOUT_TRIMAL=600      # 10 min
export TIMEOUT_QUICKTREE=900   # 15 min
export TIMEOUT_TREESHRINK=900  # 15 min
export TIMEOUT_FINALTREE=7200  # 120 min

# Timestamp used in joblog/summary names so reruns never overwrite old logs
export RUN_TAG=$(date +%Y%m%d_%H%M%S)

EXCLUDE_SAMPLES=()   # sample names to drop before step 1
export TRIMAL_THRESH=0.5
export BS_COLLAPSE_THRESH=30 # Bootstrap support below which internal nodes are collapsed
export TREESHRINK_THRESH=0.05
export MIN_TAXA=4 # Minimum taxa (sequences) and alignment length required for a meaningful tree
export MIN_ALN_LEN=100
export CONSTRAINT_TREE_STR=""   # "--tree-constraint <file>" (intentionally unquoted when used)
GEN_CONSTRAINT_TREE=0
RERUN_WITH_CONSTRAINT_TREE=0
CONSTRAINT_TREE_FILE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -gc|--gen-constraint-tree) GEN_CONSTRAINT_TREE=1; shift 1 ;;
        -rr|--rerun-with-constraint) RERUN_WITH_CONSTRAINT_TREE=1; shift 1 ;;
        -cp|--constraint-tree-path)
            case "$2" in /*) CONSTRAINT_TREE_FILE="$2" ;; *) CONSTRAINT_TREE_FILE="$WORKDIR/$2" ;; esac
            CONSTRAINT_TREE_STR="--tree-constraint $CONSTRAINT_TREE_FILE"; shift 2 ;;
        -wd|--working-directory-name) WORKDIR="$ROOTDIR/$2"; shift 2 ;;
        -bt|--bs-collapse-thresh) BS_COLLAPSE_THRESH="$2"; shift 2 ;;
        -tp|--treeshrink-thresh) TREESHRINK_THRESH="$2"; shift 2 ;;
        -tt|--trimal-thresh) TRIMAL_THRESH="$2"; shift 2 ;;
        -e|--exclude-samples) shift 1; EXCLUDE_SAMPLES=("$@"); break ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# -rr needs a constraint tree to resubmit with
if [[ "$RERUN_WITH_CONSTRAINT_TREE" = 1 ]]; then GEN_CONSTRAINT_TREE=1; fi

if [[ -n "$CONSTRAINT_TREE_FILE" && ! -s "$CONSTRAINT_TREE_FILE" ]]; then
    echo "ERROR: constraint tree $CONSTRAINT_TREE_FILE not found or empty." >&2
    exit 1
fi

# A constrained pass must NOT reuse the unconstrained step-5 outputs/markers,
# otherwise every gene is "already done" and the rerun does nothing.
if [[ -n "$CONSTRAINT_TREE_STR" ]]; then CTAG="_constrained"; else CTAG=""; fi

# Conda activation
conda activate "$HYBPIPER_ENV"

for tool in parallel mafft trimal raxml-ng run_treeshrink.py astral4 nw_ed nw_labels nw_topology; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool not found in PATH (nw_* tools: 'conda install -c bioconda newick_utils')." >&2
        exit 1
    fi
done

# TreeShrink filtering helper (used in parallel).
# output.txt has one line per tree with ALL removed taxa tab-separated, so every
# field must be read. Output is written to a .part file and renamed only if the
# write succeeded, so a failed job never leaves a half-written *_clean.fa.
filter_treeshrink() {
    local gdir=$1
    local gene=$(basename "$gdir")
    local orig_aln="$2/${gene}_trimmed.fa"
    local out_aln="$3/${gene}_clean.fa"

    [[ -s "$out_aln" ]] && return 0
    [[ -s "$orig_aln" ]] || return 0

    if [[ -s "$gdir/output.txt" ]]; then
        awk 'NR==FNR{for(i=1;i<=NF;i++) remove[$i]; next} /^>/{header=substr($1,2); keep=!(header in remove)} keep' \
            "$gdir/output.txt" "$orig_aln" > "$out_aln.part" && mv "$out_aln.part" "$out_aln"
    else
        cp "$orig_aln" "$out_aln.part" && mv "$out_aln.part" "$out_aln"
    fi
}
export -f filter_treeshrink

# Sample-exclusion helper (used in parallel). Headers look like
# ">SAMPLE_paralog_N"; strip the suffix and compare EXACTLY (no regex).
filter_excluded() {
    local f=$1
    local gene=$(basename "$f" .fasta)
    local out="$FILTER_INDIR/${gene}.fasta"
    awk 'NR==FNR{ex[$1]; next}
         /^>/{h=substr($1,2); sub(/_paralog_[0-9]+$/,"",h); keep=!(h in ex)}
         keep' "$LOG_DIR/excluded_samples.txt" "$f" > "$out"
    if [ ! -s "$out" ]; then
        echo "$gene: no sequences remain after sample exclusion" >> "$LOG_DIR/0_empty_after_filter.log"
    fi
}
export -f filter_excluded

# --------------------------------------------------------------------------
# PROCESSING STEPS
# --------------------------------------------------------------------------

echo "======================================================"
echo "Run tag: $RUN_TAG"
echo "Input genes: $GENE_FASTA_DIR"
echo "TrimAl gap threshold (-gt, min fraction of non-gap): ${TRIMAL_THRESH}"
echo "TreeShrink quantile: ${TREESHRINK_THRESH}"
echo "Collapsing nodes below support thresh: ${BS_COLLAPSE_THRESH}"
echo "Constraint tree: ${CONSTRAINT_TREE_FILE:-none}"
if [[ ${#EXCLUDE_SAMPLES[@]} -gt 0 ]]; then
    echo "Excluding samples: ${EXCLUDE_SAMPLES[*]}"
fi
echo "======================================================"
if [[ "$GEN_CONSTRAINT_TREE" = 1 ]]; then
    echo -n "Generating a constraint tree"
    if [[ "$RERUN_WITH_CONSTRAINT_TREE" = 1 ]]; then
        echo -n " and re-running"
    fi
    echo ""
    echo "======================================================"
fi

RAW_INDIR="$GENE_FASTA_DIR"

if [[ ! -d "$RAW_INDIR" ]] || [[ -z "$(find "$RAW_INDIR" -maxdepth 1 -name '*.fasta' -size +0 -print -quit)" ]]; then
    echo "ERROR: no non-empty *.fasta files in $RAW_INDIR (set GENE_FASTA_DIR to the filter script's output folder)." >&2
    exit 1
fi

export FILTER_INDIR="0_sample_filtered"
export ALIGN_DIR="1_initial_alignments"
export TRIM_DIR="2_trimmed_alignments"
export QUICK_TREES="3_quick_trees"
export FILTERED_DIR="4_filtered_alignments"
export MULTI_GENE_TREES="5_gene_trees${CTAG}"

mkdir -p "$FILTER_INDIR" "$ALIGN_DIR" "$TRIM_DIR" \
         "$QUICK_TREES" "$FILTERED_DIR" "$MULTI_GENE_TREES"

export CONCAT_TREE="$MULTI_GENE_TREES/concat.tre"

# --------------------------
# 0. Remove excluded samples before any alignment happens
# --------------------------
if [[ ! -f $SKIP_DIR/0_sample_filter.skip ]] \
   || [[ -z "$(find "$FILTER_INDIR" -maxdepth 1 -name '*.fasta' -print -quit)" ]]; then
    if [[ ${#EXCLUDE_SAMPLES[@]} -gt 0 ]]; then
        echo "Removing ${#EXCLUDE_SAMPLES[@]} sample(s) from input FASTAs..."
        printf '%s\n' "${EXCLUDE_SAMPLES[@]}" > "$LOG_DIR/excluded_samples.txt"

        find "$RAW_INDIR" -maxdepth 1 -name "*.fasta" -size +0 \
            | parallel --joblog "$LOG_DIR/0_sample_filter.${RUN_TAG}.joblog" -j "$CPUS" filter_excluded {} || true
        awk '$7 != 0' "$LOG_DIR/0_sample_filter.${RUN_TAG}.joblog" > "$LOG_DIR/0_sample_filter.${RUN_TAG}.failed" || true
    else
        echo "No samples to exclude; symlinking input as-is..."
        find "$RAW_INDIR" -maxdepth 1 -name "*.fasta" -size +0 \
            | parallel -j "$CPUS" '[ -e "$FILTER_INDIR"/$(basename {}) ] || [ -L "$FILTER_INDIR"/$(basename {}) ] || ln -s {} "$FILTER_INDIR"/$(basename {})' || true
    fi
    touch "$SKIP_DIR/0_sample_filter.skip"
fi

INDIR="$FILTER_INDIR"

# --------------------------
# 1. MAFFT alignment
# --------------------------
if [[ ! -f $SKIP_DIR/1_mafft.skip ]]; then
    echo "Starting MAFFT alignments..."

    find "$INDIR" -name "*.fasta" -size +0 | parallel --joblog "$LOG_DIR/1_mafft.${RUN_TAG}.joblog" -j "$CPUS" '
        gene=$(basename {} .fasta)
        [ -s "$ALIGN_DIR/${gene}.fa" ] && exit 0
        if timeout "$TIMEOUT_MAFFT" mafft --thread 1 --auto {} > "$ALIGN_DIR/${gene}.fa.part" 2>> "$LOG_DIR/${gene}_mafft.log"; then
            mv "$ALIGN_DIR/${gene}.fa.part" "$ALIGN_DIR/${gene}.fa"
        else
            echo "$gene: mafft failed or timed out" >> "$LOG_DIR/1_mafft.failed"
        fi
    ' || true
    awk '$7 != 0' "$LOG_DIR/1_mafft.${RUN_TAG}.joblog" > "$LOG_DIR/1_mafft.${RUN_TAG}.failed_summary" || true
    touch "$SKIP_DIR/1_mafft.skip"
fi

# --------------------------
# 2. TrimAl Trimming (Relaxed for Target-Capture Flanks)
# --------------------------
if [[ ! -f $SKIP_DIR/2_trimal.skip ]]; then
    echo "Running TrimAl..."

    find "$ALIGN_DIR" -maxdepth 1 -name "*.fa" | parallel --joblog "$LOG_DIR/2_trimal.${RUN_TAG}.joblog" -j "$CPUS" '
        aln={}
        gene=$(basename "$aln" .fa)
        out="$TRIM_DIR/${gene}_trimmed.fa"
        [ -s "$out" ] && exit 0
        if [ -s "$aln" ]; then
            # write to .part, rename only on success, so a timeout never leaves a truncated alignment
            if timeout "$TIMEOUT_TRIMAL" trimal -in "$aln" -out "$out.part" -gt "$TRIMAL_THRESH" -fasta && [ -s "$out.part" ]; then
                mv "$out.part" "$out"
            else
                echo "$gene: trimal failed, timed out, or produced an empty alignment" >> "$LOG_DIR/2_trimal.empty"
            fi
        else
            echo "$gene: input alignment empty, skipped trimal" >> "$LOG_DIR/2_trimal.empty"
        fi
    ' || true
    awk '$7 != 0' "$LOG_DIR/2_trimal.${RUN_TAG}.joblog" > "$LOG_DIR/2_trimal.${RUN_TAG}.failed_summary" || true
    touch "$SKIP_DIR/2_trimal.skip"
fi

# --------------------------
# 3. Quick Gene Trees for TreeShrink
# --------------------------
if [[ ! -f $SKIP_DIR/3_treeshrink.skip ]]; then
    echo "Building quick gene trees & running TreeShrink..."
    find "$TRIM_DIR" -name "*_trimmed.fa" | parallel --joblog "$LOG_DIR/3_treeshrink.${RUN_TAG}.joblog" -j "$CPUS" '
        aln={}
        gene=$(basename "$aln" _trimmed.fa)
        gdir="$QUICK_TREES/${gene}"
        mkdir -p "$gdir"

        [ -f "$gdir/output.txt" ] && exit 0

        ntaxa=$(grep -c "^>" "$aln" || true)
        if [ -z "$ntaxa" ] || [ "$ntaxa" -lt "$MIN_TAXA" ]; then
            echo "$gene: only $ntaxa taxa (< $MIN_TAXA), skipped" >> "$LOG_DIR/3_treeshrink.skipped_low_taxa"
            exit 0
        fi

        alnlen=$(awk '\''/^>/ {if (seq) exit; next} {seq=seq$0} END {print length(seq)}'\'' "$aln")
        if [ -z "$alnlen" ] || [ "$alnlen" -lt "$MIN_ALN_LEN" ]; then
            echo "$gene: alignment only $alnlen bp (< $MIN_ALN_LEN), skipped" >> "$LOG_DIR/3_treeshrink.skipped_short_aln"
            exit 0
        fi

        if [ ! -f "$gdir/quick.raxml.bestTree" ]; then
            timeout "$TIMEOUT_QUICKTREE" raxml-ng --msa "$aln" --model HKY+G --threads 1 --prefix "$gdir/quick" --seed 12345 --search1 --redo $CONSTRAINT_TREE_STR > "$gdir/quick.log" 2>&1
        fi
        if [ ! -f "$gdir/quick.raxml.bestTree" ]; then
            echo "$gene: quick raxml-ng tree failed or timed out" >> "$LOG_DIR/3_treeshrink.failed"
            exit 0
        fi

        timeout "$TIMEOUT_TREESHRINK" run_treeshrink.py -m per-gene -t "$gdir/quick.raxml.bestTree" -a "$aln" -q "$TREESHRINK_THRESH" -o "$gdir" > "$gdir/treeshrink.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "$gene: treeshrink failed or timed out (continuing with unfiltered alignment)" >> "$LOG_DIR/3_treeshrink.failed"
        fi
    ' || true
    awk '$7 != 0' "$LOG_DIR/3_treeshrink.${RUN_TAG}.joblog" > "$LOG_DIR/3_treeshrink.${RUN_TAG}.failed_summary" || true
    touch "$SKIP_DIR/3_treeshrink.skip"
fi

# --------------------------
# 4. Remove Outliers from Alignments
# --------------------------
if [[ ! -f $SKIP_DIR/4_filter.skip ]]; then
    echo "Filtering outliers..."
    find "$QUICK_TREES" -mindepth 1 -maxdepth 1 -type d | parallel --joblog "$LOG_DIR/4_filter.${RUN_TAG}.joblog" -j "$CPUS" filter_treeshrink {} "$TRIM_DIR" "$FILTERED_DIR" || true
    awk '$7 != 0' "$LOG_DIR/4_filter.${RUN_TAG}.joblog" > "$LOG_DIR/4_filter.${RUN_TAG}.failed_summary" || true
    touch "$SKIP_DIR/4_filter.skip"
fi

# --------------------------
# 5. Final Gene Trees, Bootstrapping & Low-Support Branch Collapsing
# --------------------------
if [[ ! -f $SKIP_DIR/5_final_trees${CTAG}.skip ]]; then
    echo "Re-estimating gene trees (GTR+G) with bootstrapping & collapsing support <$BS_COLLAPSE_THRESH..."

    export BS_TREES="autoMRE{30}"

    find "$FILTERED_DIR" -name "*_clean.fa" | parallel --joblog "$LOG_DIR/5_final_trees${CTAG}.${RUN_TAG}.joblog" -j "$CPUS" '
        aln={}
        gene=$(basename "$aln" _clean.fa)
        gdir="$MULTI_GENE_TREES/${gene}"
        mkdir -p "$gdir"

        [ -s "$gdir/final.raxml.collapsed.tre" ] && exit 0

        ntaxa=$(grep -c "^>" "$aln" || true)
        if [ -z "$ntaxa" ] || [ "$ntaxa" -lt "$MIN_TAXA" ]; then
            echo "$gene: only $ntaxa taxa (< $MIN_TAXA), skipped" >> "$LOG_DIR/5_final_trees.skipped_low_taxa"
            exit 0
        fi

        alnlen=$(awk '\''/^>/ {if (seq) exit; next} {seq=seq$0} END {print length(seq)}'\'' "$aln")
        if [ -z "$alnlen" ] || [ "$alnlen" -lt "$MIN_ALN_LEN" ]; then
            echo "$gene: alignment only $alnlen bp (< $MIN_ALN_LEN), skipped" >> "$LOG_DIR/5_final_trees.skipped_short_aln"
            exit 0
        fi

        # --all = ML search + bootstrapping; writes final.raxml.support
        #timeout "$TIMEOUT_FINALTREE" raxml-ng --all --msa "$aln" --model GTR+G --threads 1 --prefix "$gdir/final" --seed 12345 --bs-trees "$BS_TREES" --redo $CONSTRAINT_TREE_STR > "$gdir/final.log" 2>&1
        timeout "$TIMEOUT_FINALTREE" raxml-ng --tree pars{3} --msa "$aln" --model GTR+G --threads 1 --prefix "$gdir/final" --seed 12345 --bs-trees "$BS_TREES" --redo $CONSTRAINT_TREE_STR > "$gdir/final.log" 2>&1
        if [ ! -f "$gdir/final.raxml.support" ]; then
            echo "$gene: final raxml-ng tree/bootstrap failed or timed out" >> "$LOG_DIR/5_final_trees.failed"
            exit 0
        fi

        # Collapse internal nodes with support < threshold into polytomies
        nw_ed "$gdir/final.raxml.support" "i & b<$BS_COLLAPSE_THRESH" o > "$gdir/final.raxml.collapsed.tre"
        if [ ! -s "$gdir/final.raxml.collapsed.tre" ]; then
            echo "$gene: nw_ed collapsing failed, falling back to uncollapsed support tree" >> "$LOG_DIR/5_final_trees.failed"
            cp "$gdir/final.raxml.support" "$gdir/final.raxml.collapsed.tre"
        fi
    ' || true
    awk '$7 != 0' "$LOG_DIR/5_final_trees${CTAG}.${RUN_TAG}.joblog" > "$LOG_DIR/5_final_trees${CTAG}.${RUN_TAG}.failed_summary" || true
    touch "$SKIP_DIR/5_final_trees${CTAG}.skip"
fi

if [[ -z "$(find "$MULTI_GENE_TREES" -mindepth 2 -maxdepth 2 -name final.raxml.collapsed.tre -size +0 -print -quit)" ]]; then
    echo "ERROR: no gene trees were produced in $MULTI_GENE_TREES; check $LOG_DIR." >&2
    exit 1
fi
cat "$MULTI_GENE_TREES"/*/final.raxml.collapsed.tre > "$CONCAT_TREE"

# --------------------------
# 6. ASTRAL species tree (multi-copy)
# --------------------------
# Skipped in the first pass of a two-pass (-rr) run; the constrained pass makes it.
if [[ "$RERUN_WITH_CONSTRAINT_TREE" = 0 ]]; then
    echo "Making final species tree with ASTRAL..."
    astral4 \
        -i "$CONCAT_TREE" \
        -a "$MAPPING_FILE" \
        -o "final_species${CTAG}.tre" \
        -t "$CPUS"
fi

# --------------------------
# 7. Generate ASTRAL constraint tree
# --------------------------
if [[ "$GEN_CONSTRAINT_TREE" = 1 ]]; then
    echo "Making constraint tree with ASTRAL..."
    astral4 \
        -i "$CONCAT_TREE" \
        -o "constraint_uncollapsed.tre" \
        -t "$CPUS"
    nw_ed constraint_uncollapsed.tre 'i & b<0.975' o | nw_topology -I - > constraint.tre
fi

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
N_EXPECTED=$(find "$INDIR" -name "*.fasta" -size +0 | wc -l)
N_TREES=$(find "$MULTI_GENE_TREES" -name "final.raxml.collapsed.tre" | wc -l)
echo "======================================================"
echo "Pipeline complete: $N_TREES / $N_EXPECTED gene trees produced."
echo "Check $LOG_DIR/*.failed_summary and *.skipped_* logs for genes that dropped out."

# Second pass: resubmit with the constraint tree (only when -rr was requested)
if [[ "$RERUN_WITH_CONSTRAINT_TREE" = 1 ]]; then
    if [[ ! -s "$WORKDIR/constraint.tre" ]]; then
        echo "ERROR: constraint.tre was not produced; not resubmitting." >&2
        exit 1
    fi
    echo "Re-running script with new constraint tree..."
    resub_args=(-cp "$WORKDIR/constraint.tre"
                -bt "$BS_COLLAPSE_THRESH"
                -tp "$TREESHRINK_THRESH"
                -tt "$TRIMAL_THRESH"
                -wd "$(basename "$WORKDIR")")
    if [[ ${#EXCLUDE_SAMPLES[@]} -gt 0 ]]; then
        resub_args+=(-e "${EXCLUDE_SAMPLES[@]}")   # must stay last
    fi
    sbatch "$0" "${resub_args[@]}"
fi