# Target-Capture Phylogenomics Pipeline (HybPiper to ASTRAL)

A set of Bash/Slurm scripts that take raw paired-end target-capture reads (e.g. Angiosperms353 / hybrid-capture data) through read cleaning, HybPiper assembly, paralog/low-copy gene filtering, and multi-copy-aware species-tree inference.

```
raw reads
   │  01_trim_norm.sh              fastp trimming  (+ optional BBNorm normalization)
   ˅
trimmed / normalized reads
   │  02_trunc_to_75.sh            (optional) downsample + truncate to 75 bp, merge R1/R2 → single-end
   ˅
reads
   │  03_hybpiper_assemble.sh      HybPiper assemble + supercontig retrieval, per sample
   ˅
<sample>/ + <sample>_supercontig.fasta
   │  04_stats_filter_low_copy.sh  gene/sample presence filtering, per-gene FASTAs, coverage stats
   ˅
filtered_genes/<gene>.fasta
   │  05_supercontig_pipeline.sh   MAFFT → TrimAl → TreeShrink → RAxML-NG → ASTRAL (optional constraint tree)
   ˅
final_species.tre
```

All scripts are designed to run on a Slurm cluster node (they read `SLURM_CPUS_ON_NODE` / `SLURM_MEM_PER_NODE` and fall back to `nproc` / `/proc/meminfo` otherwise). They are **resumable**: finished outputs are detected and skipped on rerun.

---

## Contents

| Script | Purpose |
|---|---|
| `01_trim_norm.sh` | Adapter/quality trimming with fastp; optional depth normalization with BBNorm (via Apptainer) |
| `02_trunc_to_75.sh` | Optional: randomly subsample reads and truncate to 75 bp; R1 + R2 are concatenated into one single-end file (use this to match older/degraded samples) |
| `03_hybpiper_assemble.sh` | Runs `hybpiper assemble` and `retrieve_sequences supercontig` for every sample in parallel |
| `04_stats_filter_low_copy.sh` | Filters genes and samples by presence/length, writes per-gene FASTAs and optional coverage matrix |
| `05_supercontig_pipeline.sh` | Alignment, trimming, outlier removal, gene trees, ASTRAL species tree, optional constrained second pass |

---

## Requirements

**Cluster / general**
- Slurm (`sbatch`), Bash ≥ 4 (uses `mapfile`)
- [GNU Parallel](https://www.gnu.org/software/parallel/) (all scripts)
- Standard tools: `awk` (**gawk** is needed for the wide CSV in script 04, which uses `asorti`), `find`, `sort`, `realpath`

**Per script**

| Script | Tools |
|---|---|
| 01 | `fastp` (loaded with `module load fastp`); for normalization, `apptainer` and internet access on first run to pull `docker://staphb/bbtools:latest` |
| 02 | `seqtk`, `pigz`, `parallel` |
| 03 | [HybPiper v2](https://github.com/mossmatters/HybPiper) (`hybpiper`), BWA, and their dependencies, in a conda env |
| 04 | `parallel`, `gawk` |
| 05 | In the same conda env as HybPiper: `mafft`, `trimal`, `raxml-ng`, `treeshrink` (`run_treeshrink.py`), `astral4`, `newick_utils` (`nw_ed`, `nw_labels`, `nw_topology`) |

Example environment:

```bash
conda create -n hybpiper -c bioconda -c conda-forge \
    hybpiper mafft trimal raxml-ng treeshrink newick_utils parallel seqtk pigz
# astral4 (ASTER package): https://github.com/chaoszhang/ASTER  -> make sure `astral4` is on PATH
```

---

## Before you run: edit these paths

Some paths are hard-coded as placeholders and must be changed for your system:

| Script | Variable | Default | Notes |
|---|---|---|---|
| 03 | `PATH` export (line near top) | `$HOME/.conda/envs/hybpiper/bin` | Point to your conda env's `bin` |
| 04 | `fasta_dir` | `$HOME/hybpiper_out` | Overridable with `-f` |
| 05 | `HYBPIPER_ENV` | `$HOME/.conda/envs/hybpiper` | Conda env to activate |
| 05 | `ROOTDIR` | `$HOME` | Base directory for everything below |
| 05 | `GENE_FASTA_DIR` | `$ROOTDIR/filtered_genes` | Output of script 04 (override via `--export`) |
| 05 | `MAPPING_FILE` | `$ROOTDIR/data/mapping.txt` | Sample → species map for ASTRAL |
| 05 | `WORKDIR` | `$ROOTDIR/analysis_pipeline` | All results and logs go here |

`01`, `03` and `05` write Slurm logs to `slurm_output/`; this directory must exist **before** `sbatch` is called, or the job will fail without a log:

```bash
mkdir -p slurm_output
```

> The scripts contain no `#SBATCH` headers. Pass resources on the command line (examples below) or add your own headers. They assume they own the node (they spawn as many workers as `CPUS / threads-per-job`), so request a full node, e.g. `--exclusive` or `--cpus-per-task=<all cores>`.

---

## Step-by-step usage

### 01 — Trim and normalize reads

```bash
mkdir -p slurm_output
sbatch --nodes=1 --exclusive --time=24:00:00 -o slurm_output/%x_%j.out \
    01_trim_norm.sh -i /path/to/raw_reads
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--reads-dir` | Directory of `*.fastq.gz` / `*.fq.gz` reads (**required**) | — |
| `-s`, `--skip-norm` | Trim only; skip BBNorm | off |
| `-c`, `--target-coverage` | BBNorm target depth | `50` |
| `-td`, `--trimmed-dir` | Trimmed output directory | `<reads_dir>_trimmed` |
| `-nd`, `--normalized-dir` | Normalized output directory | `<reads_dir>_normalized` |

Notes:
- Paired files are detected by `_R1` / `_R2` in the filename; files without `_R1` are treated as unpaired.
- fastp runs with adapter auto-detection for PE, overlap correction, and a minimum length of 36 bp. HTML/JSON QC reports go to `<trimmed_dir>/fastp_reports/`.
- Output filenames are identical to input filenames, just in the new directories.
- Each task uses 8 threads; heap per BBNorm job is computed from node RAM (95 % divided among concurrent jobs, minimum 4 GB).
- The BBTools image is cached in `$SIF_DIR` (default `~/apptainer_images`).

### 02 — (Optional) Downsample and truncate to 75 bp single-end

Useful for making short-read, uniform-length, uniform-depth datasets (e.g. to match older/degraded samples).

```bash
sbatch --nodes=1 --exclusive --time=04:00:00 -o slurm_output/%x_%j.out \
    02_trunc_to_75.sh -i /path/to/reads_trimmed -f 0.6
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--in-dir` | Directory with trimmed reads (**required**) | — |
| `-o`, `--out-dir` | Output directory | `<in_dir>_75bp_SE` |
| `-f`, `--downsample-fraction` | Fraction of reads to keep, in (0, 1] | `0.6` |

Notes:
- **Input files must be named `<sample>_R1_trimmed.fastq.gz` and `<sample>_R2_trimmed.fastq.gz`.** Output is `<sample>_trimmed.fastq.gz`, where R1 and R2 (each subsampled with fixed seeds, quality-trimmed with `seqtk trimfq -q 0` and truncated to 75 bp) are concatenated into a single single-end file.
- Samples with an existing non-trivial output are skipped; the job halts on the first failure (`--halt now,fail=1`).
- To calculate the total proportion of reads kept after running this script, take (read_length - 75) ⋅ downsample-fraction, (e.g., for 150 bp reads using 0.6 as the ds fraction: (150 - 75) ⋅ 0.6 = 0.3).

### 03 — HybPiper assembly

```bash
sbatch --nodes=1 --exclusive --time=48:00:00 -o slurm_output/%x_%j.out \
    03_hybpiper_assemble.sh \
        -i /path/to/reads_for_assembly \
        -t /path/to/targets.fasta \
        -o hybpiper_out
```

| Option | Description | Default |
|---|---|---|
| `-i`, `--reads-dir` | Directory of `*.fastq.gz` / `*.fq.gz` (**required**) | — |
| `-t`, `--target-file` | Target DNA FASTA (`-t_dna`) (**required**) | — |
| `-o`, `--output-dir` | Where sample folders and supercontig FASTAs are written | `./hybpiper_out` |
| `-tj`, `--threads-per-job` | Threads per HybPiper process | `8` |

Notes:
- Sample IDs are derived from filenames: the last `_R1`/`_R2` token and everything after it is removed, as are trailing `_trimmed`, `_trim`, `_normalized`, `_norm` suffixes. E.g. `Sample_A_R1_trimmed.fastq.gz` → `Sample_A`; `Sample_A_trimmed.fastq.gz` → `Sample_A`.
- R1 files with a matching R2 are assembled as paired-end; everything else is treated as single-end.
- Runs `hybpiper assemble --bwa --force_overwrite`, then `hybpiper retrieve_sequences supercontig`, so you get `<sample>/` directories and `<sample>_supercontig.fasta` in the output directory.
- A sample is considered complete if `<sample>/<sample>_genes_with_seqs.txt` exists and is non-empty. Rerunning only processes unfinished samples, **smallest first**.
- A failed sample is logged and skipped; it does not stop the run.
- The RAM-per-job figure printed at startup is informational only and is not enforced.

### 04 — Filter genes and samples

Run this from anywhere; it `cd`s into the HybPiper output directory.

```bash
sbatch --nodes=1 --exclusive --time=02:00:00 -o slurm_output/%x_%j.out \
    04_stats_filter_low_copy.sh -f hybpiper_out -o filtered_genes
```

| Option | Description | Default |
|---|---|---|
| `-f`, `--fasta-dir` | HybPiper output directory containing `*_supercontig.fasta` | `~/hybpiper_out` |
| `-o`, `--outdir` | Output directory (**its contents are deleted at start**) | `filtered_genes` |
| `-gs`, `--min-gene-samples` | % of samples a gene must be present in | `75` |
| `-sg`, `--min-sample-genes` | % of retained genes a sample must have | `30` |
| `-c`, `--coverage-threshold` | A copy counts as "present" if its length ≥ this % of the gene's mean length | `25` |
| `-e`, `--exception-prefix` | Samples whose name starts with this prefix are never removed | none |
| `-s`, `--generate-stats-csv` | Write coverage CSV(s); alone it summarizes unfiltered data and exits without writing gene FASTAs | off |

How it filters:
1. Genes flagged by HybPiper (`*genes_with_long_paralog_warnings.txt`, `*genes_derived_from_putative_chimeric_stitched_contig.csv`) are excluded.
2. For each gene/sample, the longest copy is kept; it counts as present only if it reaches `-c` % of that gene's mean length.
3. **Genes are filtered first** (present in ≥ `-gs` % of all samples), **then samples** are judged against the retained genes (≥ `-sg` %).
4. One FASTA per retained gene is written (header = sample name, one sequence per sample).

Outputs in `filtered_genes/`:
- `<gene>.fasta`: input for script 05
- `stats.txt`, `genes_passed_list.txt`, `final_samples_list.txt`, `removed_samples.txt`
- `gene_sample_coverage_wide.csv` (with `-s`; gene × sample sequence lengths)

Quick baseline stats without filtering:

```bash
bash 04_stats_filter_low_copy.sh -f hybpiper_out -o stats_only -s
```

### 05 — Gene trees and species tree

First create the sample → species mapping file used by ASTRAL (two columns, whitespace-separated: `sample<TAB>species`). Samples that belong to the same species share a species name.

```bash
sbatch --nodes=1 --exclusive --time=72:00:00 -o slurm_output/%x_%j.out \
    --export=ALL,GENE_FASTA_DIR=$PWD/filtered_genes,MAPPING_FILE=$PWD/data/mapping.txt \
    05_supercontig_pipeline.sh
```

| Option | Description | Default |
|---|---|---|
| `-gc` | Also build an ASTRAL constraint tree (`constraint.tre`) | off |
| `-rr` | Two-pass mode: unconstrained gene trees → constraint tree → automatic `sbatch` resubmission with the constraint (implies `-gc`) | off |
| `-cp <file>` | Use an existing constraint tree (absolute path, or relative to `$WORKDIR`) | none |
| `-bt <n>` | Collapse gene-tree nodes with bootstrap support below *n* | `30` |
| `-tp <q>` | TreeShrink quantile | `0.05` |
| `-tt <x>` | TrimAl `-gt` (minimum fraction of non-gap sites per column) | `0.5` |
| `-e <s1> <s2> …` | Samples to exclude; **must be the last option** | none |

Examples:

```bash
# Simple run, excluding two samples
sbatch ... 05_supercontig_pipeline.sh -e SampleX SampleY

# Two-pass run with constraint tree (the script resubmits itself)
sbatch ... 05_supercontig_pipeline.sh -rr -bt 30
```

Pipeline stages (each writes into `$WORKDIR`):

| Directory | Stage |
|---|---|
| `0_sample_filtered/` | Excluded samples removed (or symlinks to input) |
| `1_initial_alignments/` | MAFFT `--auto` alignments |
| `2_trimmed_alignments/` | TrimAl-trimmed alignments |
| `3_quick_trees/` | Quick RAxML-NG (HKY+G) trees + TreeShrink |
| `4_filtered_alignments/` | Alignments with TreeShrink outliers removed (`*_clean.fa`) |
| `5_gene_trees[_constrained]/` | Final GTR+G RAxML-NG trees with bootstraps; low-support nodes collapsed; `concat.tre` |
| `final_species[_constrained].tre` | ASTRAL species tree |
| `constraint.tre` | With `-gc`/`-rr`: ASTRAL topology with nodes under 0.975 support collapsed |
| `logs/`, `skip/` | Logs and per-step completion markers |

Genes with fewer than 4 taxa or alignments shorter than 100 bp are skipped at the tree-building steps (see `MIN_TAXA`, `MIN_ALN_LEN` near the top of the script). Per-step timeouts are also set there.

**Resuming / rerunning a step.** The script never deletes anything. Each finished step leaves a marker in `skip/`, and each gene is skipped if its output already exists. To force a step to rerun, move its output directory aside and delete the matching `skip/*.skip` file.

---

## End-to-end example

```bash
mkdir -p slurm_output
SB="sbatch --nodes=1 --exclusive -o slurm_output/%x_%j.out"

$SB --time=24:00:00 01_trim_norm.sh -i raw_reads -s                  # trim only
$SB --time=04:00:00 02_trunc_to_75.sh -i raw_reads_trimmed           # optional
$SB --time=48:00:00 03_hybpiper_assemble.sh -i raw_reads_trimmed -t targets.fasta -o hybpiper_out
$SB --time=02:00:00 04_stats_filter_low_copy.sh -f hybpiper_out -o filtered_genes
$SB --time=72:00:00 --export=ALL,GENE_FASTA_DIR=$PWD/filtered_genes,MAPPING_FILE=$PWD/mapping.txt \
    05_supercontig_pipeline.sh -rr
```

Chain jobs with `--dependency=afterok:<jobid>` if you want them to run unattended.

---

## Tips and troubleshooting

- **Job fails with no log:** `slurm_output/` didn't exist at submission time.
- **Nothing happens on rerun:** that's the resume logic working. Check for `[Skip]` messages, output files, or `.skip` markers.
- **Samples missing downstream:** check `removed_samples.txt` (script 04), HybPiper failures in the Slurm log (script 03), and `logs/*.failed_summary` / `*.skipped_*` (script 05).
- **Paired reads not detected:** the pairing logic relies on `_R1` / `_R2` in the filename.
- **Command not found inside the job:** make sure the conda env path settings above are correct and the env contains every tool listed in the requirements.

## Citation

If you use this pipeline, please cite the underlying tools: fastp, BBNorm/BBTools, HybPiper, BWA, MAFFT, TrimAl, TreeShrink, RAxML-NG, and ASTRAL/ASTER.
