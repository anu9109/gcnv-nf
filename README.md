# gcnv-nf

Nextflow pipeline for germline CNV detection in a single sample using two independent callers, followed by merging, annotation, and prioritization.

## Overview

The pipeline runs two CNV callers in parallel, merges their outputs, annotates with AnnotSV, filters for clinically relevant events, and generates read-depth coverage plots for each priority event.

```
BAM + cohort BAMs
       │
       ├──► A: cn.mops          (CNMOPS subworkflow)
       │         └─ CNVs → VCF
       │
       ├──► B: GATK gCNV        (GATK_GCNV subworkflow)
       │         └─ genotyped segments → VCF
       │
       ├──► C: SURVIVOR         merge VCFs from A + B
       │
       ├──► D: AnnotSV          annotate merged VCF
       │
       ├──► E: Filter           keep "full" events where ranking ≥ 4 OR SUPP > 1
       │
       └──► F: Coverage plots   one PDF per priority event
```

The GATK gCNV branch (B) expects a pre-built reference genome (fasta + index + dict + interval lists)
at `--ref_path`. That genome is produced separately by the standalone `PREPARE_GENOME` workflow — see
[Standalone: Genome Preparation](#standalone-genome-preparation) below.

## Usage

```bash
nextflow run main.nf \
    --sample_id   SAMPLE_ID \
    --bam_file    /path/to/sample.bam \
    --bams_list   /path/to/cohort_bams.txt \
    --ref_path    /path/to/genome_dir \
    --depth_file  /path/to/sample.regions.bed.gz \
    --outdir      /path/to/outdir \
    -profile hpc_slurm
```

`--ref_path` points to a directory containing the pre-built genome files expected by GATK gCNV:
`gr37_clean.fasta`, `.fasta.fai`, `.dict`, `.interval_list`, and `_annotated.interval_list`.

### Standalone: Genome Preparation

If you need to build the `--ref_path` genome from a new/different raw fasta, run the `PREPARE_GENOME`
entry workflow once, ahead of the main pipeline:

```bash
nextflow run main.nf -entry PREPARE_GENOME \
    --genome_fasta /path/to/raw_genome.fasta \
    --outdir       /path/to/outdir \
    -profile hpc_slurm
```

This subsets the fasta to chromosomes 1-22/X/Y, indexes it, creates a sequence dictionary, and builds
GATK interval lists (with GC/mappability/segmental-duplication annotations), publishing everything
directly to `outdir/`. Point `--ref_path` at that directory for subsequent main pipeline runs.

### Re-running visualization only

If you already have cn.mops and GATK gCNV VCFs from a previous run, skip straight to merging/annotation/plots:

```bash
nextflow run main.nf --rerun_viz \
    --sample_id  SAMPLE_ID \
    --depth_file /path/to/sample.regions.bed.gz \
    --outdir     /path/to/outdir \
    -profile hpc_slurm
```

## Parameters

| Parameter | Description |
|---|---|
| `sample_id` | [Required] Sample identifier |
| `bam_file` | [Required] BAM file for the sample |
| `ref_path` | [Required] Directory with the pre-built genome files for GATK gCNV (see [Standalone: Genome Preparation](#standalone-genome-preparation)) |
| `depth_file` | [Required] Per-base or binned depth BED file (for coverage plots) |
| `outdir` | [Required] Output directory |
| `annotsv_bin` | [Required] Path to AnnotSV bin |
| `annotsv_dir` | [Required] Path to AnnotSV annotations |
| `genome_fasta` | [Required for `-entry PREPARE_GENOME` only] Raw genome fasta to subset/preprocess into `ref_path` |
| `bams_list` | [Optional] Text file listing BAM paths for the cn.mops cohort (default: `ref_path/cnmops_reference_samples.txt`) |
| `genome_chrs` | [Optional] GATK gCNV: genome_chrs.txt or equivalent (default under `ref_path`) |
| `mappability_bed` | [Optional] GATK gCNV: hg19.nochr.k100.umap.single.merged.bed.gz or equivalent (default under `ref_path`) |
| `segmental_duplication_bed` | [Optional] GATK gCNV: hg19.nochr.SegDups.elements.no_gl.bed or equivalent (default under `ref_path`) |
| `gc_file` | [Optional] GATK gCNV: gr37_clean.gc_content.bed or equivalent (default under `ref_path`) |
| `map_file` | [Optional] GATK gCNV: hg19.nochr.k100.umap.single.merged.bed.gz or equivalent (default under `ref_path`) |
| `scatter_count` | [Optional] Number of genomic shards for GATK gCNV (default: 65) |
| `model_ploidy_outdir` | [Optional] Pre-trained GATK ploidy model directory (default under `ref_path`) |
| `model_cnvs_outdir` | [Optional] Pre-trained GATK CNV model directory (default under `ref_path`) |
| `bedtools` | [Optional] Path to bedtools executable (default: `bedtools` resolved from `$PATH`) |

## Outputs

| Path | Contents |
|---|---|
| `outdir/` | Pre-built genome files from the standalone `PREPARE_GENOME` workflow (`gr37_clean.*`) |
| `outdir/cnmops/` | cn.mops segmentation, CNVs, CNVRs, and VCF |
| `outdir/gatk_gcnv/` | Genotyped segments and intervals VCFs, denoised copy ratios |
| `outdir/survivor/` | Merged VCF |
| `outdir/annotsv/` | Full AnnotSV-annotated TSV |
| `outdir/*.priority.tsv` | Filtered events (AnnotSV ranking ≥ 4 or SUPP > 1, full records only) |
| `outdir/plots/` | Per-event read-depth coverage PDFs |

## Filtering criteria (step E)

Events are retained if both conditions are met:
- `AnnotSV type` == `full`
- `AnnotSV ranking` ≥ 4 (Likely Pathogenic or Pathogenic) **OR** `SUPP` > 1 (supported by both callers)

## Dependencies

- [Nextflow](https://www.nextflow.io/) ≥ 22.10
- Singularity (containers defined in `nextflow.config`; most are pulled on demand via Docker/[Wave](https://seqera.io/wave/))
- [cn.mops](https://bioconductor.org/packages/cn.mops/) — `ghcr.io/anu9109/cnmops` container
- [seqtk](https://github.com/lh3/seqtk) — `biocontainers/seqtk` container (standalone genome prep only)
- [GATK](https://gatk.broadinstitute.org/) 4.6.1.0 — `broadinstitute/gatk` container
- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) + bcftools — Wave-built container
- [AnnotSV](https://lbbe-software.github.io/AnnotSV/) + bedtools — run natively on the host (not containerized)
- R with `ggplot2`, `data.table`, `tidyverse` (for coverage plots) — `ghcr.io/anu9109/gcnv_viz` container

