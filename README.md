# CoMiCall

**Co**nsensus **Mi**tochondrial somatic variant **Ca**lling (**CoMiCall**) from clonal WGS data.

A Nextflow DSL2 pipeline that derives high-confidence mitochondrial single-nucleotide variants (MT-SNVs) from single-cell derived colony Whole Genome Sequencing data (BAMs files). Each sample (colony) is independently filtered to MT reads, realigned against a NUMT-masked reference, scanned read-by-read for strand-aware base calls, and aggregated into per-donor allelic / coverage tables. A final filtering step yields the confident MT-SNV set together with QC metrics.

## Pipeline

```
INDEX_BAM → FILTER_MT_READS → REALIGN_MT_READS → MAKE_TABLES ─┐
                                                              ├─ GATHER_TABLES → FILTER_MUTS
                              INDEX_BAM ── ESTIMATE_MT_CN ────┘
```

The MT contig (`chrM` or `MT`) is auto-detected from each BAM header / reference `.fai` — no `--region` flag is needed.

## Usage

```bash
nextflow run main.nf \
    -params-file params/params_test.json \
    -profile conda,local
```

For an HPC run with containers and LSF:

```bash
nextflow run main.nf \
    -params-file params/params.json \
    -profile singularity,lsf \
    -c config/user.config
```

## Profiles

| Profile        | Effect                                              |
| -------------- | --------------------------------------------------- |
| `conda`        | Use Conda environments under `envs/`                |
| `singularity`  | Use Apptainer / Singularity bioconda containers     |
| `local`        | Run on the current machine                          |
| `lsf`          | Submit jobs to LSF (`team274`, `normal` queue)      |

## Parameters

All parameters can be set via `--<name> <value>` on the CLI or through a JSON file passed with `-params-file`.

### Input / output (required)

| Param             | Description                                                |
| ----------------- | ---------------------------------------------------------- |
| `input_folder`    | Directory containing per-colony deduplicated `*.bam` files |
| `individual`      | Donor identifier (used in output filenames)                |
| `output_folder`   | Output directory                                           |

### Reference (required)

| Param         | Description                                                                                  |
| ------------- | -------------------------------------------------------------------------------------------- |
| `ref`         | Faidx-indexed reference FASTA containing `chrM` (or `MT`); used by `MAKE_TABLES`             |
| `masked_ref`  | NUMT-masked, BWA-indexed full-genome FASTA (see `code/mask_numts.sh`); used by `REALIGN_MT_READS` |

### Quality / depth (defaults shown)

| Param           | Default | Description                                                       |
| --------------- | ------- | ----------------------------------------------------------------- |
| `min_mq`        | `30`    | Minimum mapping quality                                           |
| `min_bq`        | `30`    | Minimum base quality                                              |
| `max_softclip`  | `0.1`   | Max fraction of softclipped bases per read (`MAKE_TABLES` filter) |

### `FILTER_MUTS` thresholds (defaults shown)

| Param                       | Default | Description                                                       |
| --------------------------- | ------- | ----------------------------------------------------------------- |
| `min_strand_ratio`          | `0.1`   | Lower bound on per-call fwd/rev strand ratio                      |
| `max_strand_ratio`          | `0.9`   | Upper bound on per-call fwd/rev strand ratio                      |
| `af_threshold`              | `0.05`  | AF cutoff splitting high-AF and low-AF stages                     |
| `min_callable_coverage`     | `10`    | Minimum callable coverage at the site                             |
| `max_sb_pval`               | `0.05`  | Minimum strand-bias p-value (calls below this are dropped)        |
| `max_prevalence_low_AF`     | `0.25`  | Max cross-sample prevalence allowed for low-AF calls              |
| `min_gs_count_low_AF`       | `1`     | Min `gs>1_count` evidence required for low-AF calls               |
| `min_overlap_count_low_AF`  | `1`     | Min `overlap_count` evidence required for low-AF calls            |
| `min_n_mutations`           | `1`     | Min number of MT-SNVs per sample for outlier-trim retention       |
| `max_prevalence_germline`   | `0.9`   | Drop variants above this prevalence as germline (Stage III)       |

## Outputs

Published under `${output_folder}/${individual}/`:

- `${individual}.allelic_table.tsv.gz` — pooled per-donor allelic table
- `${individual}.coverage_table.tsv.gz` — pooled per-donor coverage table
- `${individual}.mt_cn.tsv.gz` — per-colony mtDNA copy-number estimates
- `${individual}.confident_allelic_table.tsv.gz` — post Stage I+II filtered calls
- `${individual}.final_allelic_table.tsv.gz` — confident calls with germline removed
- `${individual}.metrics.txt` — long-format `donor / metric / value` QC summary

## Stub / dry run

```bash
nextflow run main.nf -params-file params/params_test.json -profile local -stub
```

## Help

```bash
nextflow run main.nf --help
```
