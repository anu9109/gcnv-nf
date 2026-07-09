#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPARE_SAMPLE_LIST } from './modules/cnmops.nf'
include { RUN_CNMOPS } from './modules/cnmops.nf'
include { CNMOPS_TO_VCF } from './modules/cnmops.nf'


workflow {
    
    CNMOPS(
        params.sample_id, 
        file(params.bam_file),
        file(params.bams_list),
        params.outdir
    )
    
    //GATK_GCNV()
}



// Subworkflow: CNMOPS CNV Detection
//
// This subworkflow performs CNV detection using the cn.mops algorithm.
// It processes individual samples by normalizing against a cohort of control samples.
workflow CNMOPS {
    take:
        sample_id      // val: sample identifier
        bam_file       // val: path to sample BAM file
        bams_list      // val: path to file with list of all BAM files
        outdir         // val: output directory path
        //cnmops_sif     // val: path to cnmops singularity container
        //run_cnmops_r   // val: path to run_cnmops.R script
        //cnmops_to_vcf_r // val: path to cnmops_to_vcf_per_sample.R script
        //annotsv_bin    // val: path to annotsv binary
        //annotsv_dir    // val: path to annotsv directory
        //bedtools_bin   // val: path to bedtools binary

    main:
        // Step 1: Prepare sample list with controls
        PREPARE_SAMPLE_LIST(
            sample_id,
            bam_file,
            outdir,
            bams_list
        )

        // Step 2: Run cn.mops analysis
        RUN_CNMOPS(
            PREPARE_SAMPLE_LIST.out.bams_txt,
            //cnmops_sif,
            //run_cnmops_r,
            sample_id,
            outdir
        )

        // Step 3: Convert results to VCF format
        CNMOPS_TO_VCF(
            //cnmops_sif,
            //cnmops_to_vcf_r,
            sample_id,
            RUN_CNMOPS.out.sample_cnvs,
            outdir
        )

    emit:
        segmentation = RUN_CNMOPS.out.sample_seg
        cnvs = RUN_CNMOPS.out.sample_cnvs
        cnvr = RUN_CNMOPS.out.sample_cnvr
        vcf = CNMOPS_TO_VCF.out.sample_cnvs_vcf
}

/*
workflow.onComplete = {
    log.info(
        workflow.success
            ? "\nDone!\n cn.mops ran successfully. See results in: ${params.outdir}\n"
            : "\nOops .. something went wrong\n"
    )
}
*/


/*
//
// Subworkflow: GATK gCNV Analysis
//
// This subworkflow performs germline CNV detection using GATK's gCNV caller.
// It includes genome preprocessing, read count collection, ploidy determination,
// CNV calling, postprocessing, and joint cohort segmentation.

include { COLLECT_READ_COUNTS } from './modules/gatk.nf'
include { FILTER_GENOME } from './modules/gatk.nf'
include { SCATTER_GENOME } from './modules/gatk.nf'
include { DETERMINE_PLOIDY } from './modules/gatk.nf'
include { CALL_CNVS } from './modules/gatk.nf'
include { POSTPROCESS_CNVS } from './modules/gatk.nf'
include { JOINT_CNVS_SEGMENTATION } from './modules/gatk.nf'

workflow GATK_GCNV {
    take:
        bams_channel             // channel of [sample_id, bam_file] tuples
        scatter_count            // val: number of intervals
        genome_outdir            // val: output directory for genome
        sample_outdir            // val: output directory for samples
        cnvs_outdir              // val: output directory for CNV calls
        model_outdir             // val: output directory for models
        scatter_outdir           // val: output directory for scattered intervals
        outdir                   // val: final output directory
        gr37_fasta_in            // val: path to reference genome fasta
        seqtk                    // val: path to seqtk
        genome_chrs              // val: path to genome chromosomes file
        gatk_sif                 // val: path to gatk singularity container
        mappability_bed          // val: path to mappability BED
        segmental_duplication_bed // val: path to segmental duplication BED
        hg19_ploidy_priors       // val: path to ploidy priors file
        model_prefix             // val: prefix for model output files
        cnvs_prefix              // val: prefix for CNV output files
        pedigree                 // val: path to pedigree file

    main:

        // Step 2: Collect read counts from samples
        COLLECT_READ_COUNTS(
            bams_channel,
            sample_outdir,
            gatk_sif,
            PREPROCESS_GENOME_FASTA.out.interval_list,
            PREPROCESS_GENOME_FASTA.out.fasta
        )

        // Step 3: Filter genome intervals based on read count outliers
        COLLECT_READ_COUNTS.out.read_counts.collect().set { read_count_list }
        FILTER_GENOME(
            read_count_list,
            genome_outdir,
            gatk_sif,
            PREPROCESS_GENOME_FASTA.out.annotated_interval_list,
            PREPROCESS_GENOME_FASTA.out.interval_list
        )

        // Step 4: Scatter genome intervals into shards
        SCATTER_GENOME(
            scatter_count,
            scatter_outdir,
            gatk_sif,
            FILTER_GENOME.out.filtered_interval_list
        )

        // Step 5: Determine ploidy model - cohort mode
        DETERMINE_PLOIDY(
            COLLECT_READ_COUNTS.out.read_counts.collect(),
            model_outdir,
            model_prefix,
            gatk_sif,
            FILTER_GENOME.out.filtered_interval_list,
            hg19_ploidy_priors
        )

        // Step 6: Create channel for interval IDs
        Channel
            .from(1..scatter_count)
            .map { String.format("%04d", it) }
            .set { interval_ids }

        // Step 7: Call germline CNVs in cohort mode
        COLLECT_READ_COUNTS.out.read_counts.collect().set { all_read_counts }

        CALL_CNVS(
            all_read_counts,
            interval_ids.flatten(),
            scatter_count,
            cnvs_outdir,
            cnvs_prefix,
            gatk_sif,
            model_outdir,
            model_prefix,
            scatter_outdir
        )

        // Step 8: Create sample index and sample ID mapping
        bams_channel
            .toList()
            .flatMap { samples ->
                samples.eachWithIndex { sample, idx ->
                    tuple(idx, sample[0])
                }
            }
            .set { sample_info }

        // Step 9: Postprocess CNVs per sample
        sample_info
            .combine(Channel.value(interval_ids.toList()))
            .set { sample_intervals }

        POSTPROCESS_CNVS(
            sample_info,
            interval_ids.toList(),
            cnvs_outdir,
            cnvs_prefix,
            scatter_count,
            model_outdir,
            model_prefix,
            genome_outdir,
            gatk_sif,
            pedigree
        )

        // Step 10: Create joint CNV segmentation
        POSTPROCESS_CNVS.out.genotyped_intervals.collect().set { vcfs }


    emit:
        genome_fasta = PREPROCESS_GENOME_FASTA.out.fasta
        read_counts = COLLECT_READ_COUNTS.out.read_counts
        ploidy_calls = DETERMINE_PLOIDY.out.ploidy_calls
        cnv_calls = CALL_CNVS.out.cnv_calls
        cnv_models = CALL_CNVS.out.cnv_model
        genotyped_vcfs = POSTPROCESS_CNVS.out.genotyped_intervals
        denoised_copy_ratios = POSTPROCESS_CNVS.out.denoised_copy_ratios
}

*/








