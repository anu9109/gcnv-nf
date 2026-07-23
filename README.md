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

## Usage

```bash
nextflow run main.nf \
    --sample_id   SAMPLE_ID \
    --bam_file    /path/to/sample.bam \
    --bams_list   /path/to/cohort_bams.txt \
    --reference   /path/to/GRCh37.fa \
    --depth_file  /path/to/sample.regions.bed.gz \
    --outdir      /path/to/outdir \
    -profile hpc_slurm
```

## Parameters

| Parameter | Description |
|---|---|
| `sample_id` | Sample identifier |
| `bam_file` | BAM file for the sample |
| `bams_list` | Text file listing BAM paths for the cn.mops cohort |
| `reference` | GRCh37 reference FASTA |
| `depth_file` | Per-base or binned depth BED file (for coverage plots) |
| `outdir` | Output directory |
| `scatter_count` | Number of genomic shards for GATK gCNV (default: 65) |
| `model_ploidy_outdir` | Pre-trained GATK ploidy model directory |
| `model_cnvs_outdir` | Pre-trained GATK CNV model directory |

## Outputs

| Path | Contents |
|---|---|
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
- Singularity (containers defined in `nextflow.config`)
- [cn.mops](https://bioconductor.org/packages/cn.mops/) — via `cnmops.sif`
- [GATK](https://gatk.broadinstitute.org/) 4.6.1.0 — via `gatk_4.6.1.0.sif`
- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
- [AnnotSV](https://lbbe-software.github.io/AnnotSV/)
- R with `ggplot2`, `data.table`, `tidyverse` (for coverage plots)
