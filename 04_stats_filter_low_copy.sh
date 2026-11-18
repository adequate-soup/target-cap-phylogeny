#!/bin/bash

set -euo pipefail
export SHELL="$(type -p bash)"   # GNU parallel needs bash to run exported functions

export exception=""
gen_csv=0
summary_only=0
custom_thresholds_passed=0
fasta_dir="~/hybpiper_out"
min_gene_samples=75     # % of samples a gene must be present in
min_sample_genes=30     # % of retained genes a sample must have
coverage_threshold=25   # % of the gene's mean length a sample's copy must reach to count as "present"
cpus=${SLURM_CPUS_ON_NODE:-$(nproc)}
filtered_dir="filtered_genes"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -gs|--min-gene-samples) min_gene_samples="$2"; custom_thresholds_passed=1; shift 2 ;;
        -sg|--min-sample-genes) min_sample_genes="$2"; custom_thresholds_passed=1; shift 2 ;;
        -c|--coverage-threshold) coverage_threshold="$2"; custom_thresholds_passed=1; shift 2 ;;
        -s|--generate-stats-csv) gen_csv=1; summary_only=1; shift 1 ;;
        -o|--outdir) filtered_dir="$2"; shift 2 ;;
        -e|--exception-prefix) exception="$2"; custom_thresholds_passed=1; shift 2 ;;
        -f|--fasta-dir) fasta_dir="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Ensure outdir exists and set ABSOLUTE paths for all working directories
mkdir -p "$filtered_dir"
export filtered_dir="$(realpath "$filtered_dir")"
rm -rf "${filtered_dir:?}"/*

export fasta_dir="$(realpath "$fasta_dir")"

cd "$fasta_dir" || exit 1

if ! command -v parallel &> /dev/null; then
    echo "GNU Parallel could not be found. Please install it to run this script." >&2
    exit 1
fi

main() {
    # --- Stats file setup ---
    stats_file="$filtered_dir/stats.txt"
    echo "Filtering statistics" > "$stats_file"
    echo "===================" >> "$stats_file"
    echo "Run date: $(date)" >> "$stats_file"

    # --- Temp file cleanup trap ---
    skip_file=$(mktemp)
    all_gene_lengths_file=$(mktemp)
    present_file=$(mktemp)
    gene_counts_file=$(mktemp)
    genes_passed_file=$(mktemp)
    sample_names_file=$(mktemp)
    trap 'rm -f "$skip_file" "$all_gene_lengths_file" "$present_file" "$gene_counts_file" "$genes_passed_file" "$sample_names_file"' EXIT

    # --- Load samples (<SAMPLE>_supercontig.fasta) ---
    shopt -s nullglob
    samples=("$fasta_dir"/*_supercontig.fasta)
    shopt -u nullglob
    local num_samples=${#samples[@]}
    if [[ $num_samples -eq 0 ]]; then
        echo "No *_supercontig.fasta files found in $fasta_dir" >&2
        exit 1
    fi
    for f in "${samples[@]}"; do basename "$f" _supercontig.fasta; done > "$sample_names_file"
    echo "Found $num_samples unique samples to process."
    echo "Total samples: $num_samples" >> "$stats_file"

    # --- Load excluded genes (paralog warnings / chimeric contigs) ---
    find "$fasta_dir" -mindepth 2 -maxdepth 2 \
        \( -name '*genes_with_long_paralog_warnings.txt' \
           -o -name '*genes_derived_from_putative_chimeric_stitched_contig.csv' \) -print0 \
        | xargs -0 -r awk -F'[,[:space:]]+' 'NF && $1 != "" {print $1}' \
        | LC_ALL=C sort -u > "$skip_file"

    num_skips=$(wc -l < "$skip_file")
    echo "Number of genes flagged to skip: $num_skips"
    echo "Genes skipped (warnings): $num_skips" >> "$stats_file"

    # --- Parallel gene length calculation ---
    echo "Calculating gene lengths in parallel..."
    export skip_file

    _get_gene_lengths() {
        local fasta_file="$1"
        local sample_name
        sample_name=$(basename "$fasta_file" _supercontig.fasta)

        awk -v skipfile="$skip_file" -v sample="$sample_name" '
            BEGIN {
                while ((getline g < skipfile) > 0) skip[g]=1
                close(skipfile)
            }
            function emit() {
                if (gene_id != "" && !(gene_id in skip)) {
                    print gene_id, sample, "single_hit", length(seq)
                }
            }
            /^>/ {
                emit()
                header = $1
                sub(/^>/, "", header)
                # header may be SAMPLE-GENE; gene is the part after the last "-"
                n = split(header, parts, "-")
                gene_id = parts[n]
                seq = ""
                next
            }
            { seq = seq $0 }
            END { emit() }
        ' "$fasta_file"
    }
    export -f _get_gene_lengths

    # One line per (gene, sample): if a sample has several copies of a gene, keep the longest.
    parallel -j "$cpus" _get_gene_lengths {} ::: "${samples[@]}" \
        | LC_ALL=C sort -k1,1 -k2,2 -k4,4nr \
        | awk '!seen[$1,$2]++' > "$all_gene_lengths_file"

    # --- Handle standalone -s vs filtered execution ---
    if [[ $summary_only -eq 1 && $custom_thresholds_passed -eq 0 ]]; then
        echo "Summary mode (-s alone): Retaining 100% of detected genes and samples for baseline summary."
        echo "Mode: Unfiltered Raw Summary (-s standalone)" >> "$stats_file"

        awk '{print $1}' "$all_gene_lengths_file" | LC_ALL=C sort -u > "$genes_passed_file"
        num_genes_passed=$(wc -l < "$genes_passed_file")
        echo "Total unique genes evaluated: $num_genes_passed"
        echo "Genes evaluated: $num_genes_passed" >> "$stats_file"

        mapfile -t final_samples < "$sample_names_file"
    else
        echo "Applying filtering parameters (min_gene_samples=$min_gene_samples%, min_sample_genes=$min_sample_genes%, cov=$coverage_threshold%)..."
        echo "Mode: Filtered Run" >> "$stats_file"

        # 1. Mark a gene as "present" in a sample only if its length >= coverage_threshold %
        #    of that gene's mean length across samples.
        awk -v cov_thresh="$coverage_threshold" '
            FNR==NR { total_len[$1] += $4; count[$1]++; next }
            {
                avg = total_len[$1] / count[$1]
                if (avg > 0 && 100 * $4 / avg >= cov_thresh) print
            }
        ' "$all_gene_lengths_file" "$all_gene_lengths_file" > "$present_file"

        # 2. FILTER GENES FIRST: keep genes present in >= min_gene_samples % of ALL samples
        awk '{print $1}' "$present_file" | LC_ALL=C sort | uniq -c > "$gene_counts_file"
        min_samples=$(( (min_gene_samples * num_samples + 99) / 100 ))   # ceiling
        awk -v min="$min_samples" '$1>=min {print $2}' "$gene_counts_file" | LC_ALL=C sort > "$genes_passed_file"
        num_genes_passed=$(wc -l < "$genes_passed_file")
        echo "$num_genes_passed genes passed presence filter (at least $min_samples samples [$min_gene_samples%])."
        echo "Genes passed: $num_genes_passed (threshold: $min_samples samples [$min_gene_samples%])" >> "$stats_file"

        [[ "$num_genes_passed" -eq 0 ]] && { echo "No genes passed filter. Exiting."; exit 0; }
        cp "$genes_passed_file" "$filtered_dir/genes_passed_list.txt"

        # 3. THEN FILTER SAMPLES: judged only on the genes that passed step 2
        min_genes=$(( (min_sample_genes * num_genes_passed + 99) / 100 ))   # ceiling
        echo "Sample threshold: at least $min_genes of $num_genes_passed genes [$min_sample_genes%]" >> "$stats_file"

        removed_samples_file="$filtered_dir/removed_samples.txt"
        : > "$removed_samples_file"

        mapfile -t final_samples < <(
            awk -v min_g="$min_genes" \
                -v exc="$exception" \
                -v rem_file="$removed_samples_file" \
                -v tot="$num_genes_passed" \
                -v genes_file="$genes_passed_file" \
                -v names_file="$sample_names_file" '
                BEGIN {
                    while ((getline g < genes_file) > 0) valid_genes[g]=1
                    close(genes_file)
                    # seed every sample so samples with zero passing genes are still evaluated
                    while ((getline s < names_file) > 0) cnt[s]=0
                    close(names_file)
                }
                {
                    if (($1 in valid_genes) && !(($1,$2) in seen)) {
                        seen[$1,$2] = 1
                        cnt[$2]++
                    }
                }
                END {
                    for (s in cnt) {
                        if (length(exc) > 0 && index(s, exc) == 1) { print s; continue }
                        if (cnt[s] >= min_g) {
                            print s
                        } else {
                            print s " removed: only " cnt[s] "/" tot " genes (threshold: " min_g ")" >> rem_file
                        }
                    }
                }' "$present_file" | LC_ALL=C sort
        )
        echo "Samples removed: $(( num_samples - ${#final_samples[@]} ))" >> "$stats_file"
    fi

    final_samples_file="$filtered_dir/final_samples_list.txt"
    printf "%s\n" "${final_samples[@]}" > "$final_samples_file"
    echo "${#final_samples[@]} samples in filtered list."
    echo "Samples in final list: ${#final_samples[@]}" >> "$stats_file"

    [[ "${#final_samples[@]}" -eq 0 ]] && { echo "No samples found. Exiting."; exit 0; }

    # --- Generate coverage CSVs directly from parsed length data ---
    if [ $gen_csv -eq 1 ]; then
        echo "Creating coverage CSVs from evaluated dataset..."

        awk -v genes_file="$genes_passed_file" -v samples_file="$final_samples_file" '
            BEGIN {
                while ((getline g < genes_file) > 0) valid_genes[g]=1
                while ((getline s < samples_file) > 0) valid_samples[s]=1
            }
            ($1 in valid_genes) && ($2 in valid_samples) {
                print $1 "\t" $2 "\t" $4
            }
        ' "$all_gene_lengths_file" > "$filtered_dir/gene_sample_coverage.csv"

        # asorti requires gawk
        awk '
        {
            cov[$1,$2]=$3
            genes[$1]; samples[$2]
        }
        END{
            n=asorti(samples, sample_arr)
            printf "Gene"
            for(i=1;i<=n;i++) printf ",%s", sample_arr[i]
            print ""
            m=asorti(genes, gene_arr)
            for(i=1;i<=m;i++){
                g=gene_arr[i]
                printf "%s", g
                for(j=1;j<=n;j++){
                    s=sample_arr[j]
                    v=((g,s) in cov ? cov[g,s] : 0)
                    printf ",%s", v
                }
                print ""
            }
        }' "$filtered_dir/gene_sample_coverage.csv" > "$filtered_dir/gene_sample_coverage_wide.csv"

        rm -f "$filtered_dir/gene_sample_coverage.csv"
    fi

    # --- Early exit for -s flag ---
    if [[ $summary_only -eq 1 ]]; then
        echo "=========================================="
        echo " Summary stats completed (-s mode active) "
        echo "=========================================="
        cat "$stats_file"
        echo "Stats log saved to: $stats_file"
        echo "Coverage CSVs generated in: $filtered_dir/"
        echo "Exiting without extracting individual gene FASTA files."
        exit 0
    fi

    # --- Extract filtered sequences ---
    echo "Extracting filtered sequences..."
    export final_samples_file filtered_dir fasta_dir

    _extract_all_genes() {
        local gene="$1"
        local final_samples_file="$2"
        local outfile="$filtered_dir/${gene}.fasta"
        : > "$outfile"

        # One sequence per sample, taken from its supercontig file; header is just the sample name
        while IFS= read -r sample || [ -n "$sample" ]; do
            sample=$(echo "$sample" | tr -d '\r')
            [ -z "$sample" ] && continue
            local fasta_file="$fasta_dir/${sample}_supercontig.fasta"
            [ -f "$fasta_file" ] || continue

            awk -v gene="$gene" -v sample="$sample" '
                function emit() {
                    if (p && seq != "" && !done) {
                        print ">" sample
                        print seq
                        done = 1          # keep only the first matching copy
                    }
                }
                /^>/ {
                    emit()
                    h = $1
                    sub(/^>/, "", h)
                    n = split(h, parts, "-")
                    p = (parts[n] == gene)   # exact match on gene name
                    seq = ""
                    next
                }
                p { seq = seq $0 }
                END { emit() }
            ' "$fasta_file" >> "$outfile"
        done < "$final_samples_file"
    }
    export -f _extract_all_genes

    parallel -j "$cpus" _extract_all_genes {} "$final_samples_file" :::: "$genes_passed_file"

    echo "Done. Filtered genes are in '$filtered_dir/'"
    echo "Stats written to $stats_file"
}

main "$@"
