#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPARE_SAMPLE_LIST } from './modules/cnmops.nf'
include { RUN_CNMOPS } from './modules/cnmops.nf'
include { CNMOPS_TO_VCF } from './modules/cnmops.nf'
include { PREPROCESS_GENOME_FASTA } from './modules/gatk.nf'
include { COLLECT_READ_COUNTS } from './modules/gatk.nf'
include { FILTER_GENOME } from './modules/gatk.nf'
include { SCATTER_GENOME } from './modules/gatk.nf'
include { DETERMINE_PLOIDY_CASE } from './modules/gatk.nf'
include { CALL_CNVS_CASE } from './modules/gatk.nf'
include { POSTPROCESS_CNVS } from './modules/gatk.nf'
//include { JOINT_CNVS_SEGMENTATION } from './modules/gatk.nf'
include { SURVIVOR_MERGE } from './modules/survivor.nf'
include { ANNOTSV } from './modules/annotsv.nf'


workflow {
    
    // A: CNMOPS
    outdir_cnmops = "${params.outdir}/cnmops"
    CNMOPS(
        params.sample_id, 
        file(params.bam_file),
        file(params.bams_list),
        outdir_cnmops
    )
    
    // B: GATK_GCNV
    outdir_gatk = "${params.outdir}/gatk_gcnv"
    bams_channel = Channel.of(tuple(params.sample_id, file(params.bam_file)))
    genome_outdir = file("${outdir_gatk}/genome")
    scatter_count = params.scatter_count as int
    interval_ids = Channel
        .from(1..scatter_count)
        .map { String.format("%04d", it) }
    bams_channel
        .combine(interval_ids)
        .map { row -> tuple(row[0], row[1], row[2]) }
        .set { sample_id_intervals_ch }
    pedigree = file("${outdir_gatk}/pedigree.txt")
    GATK_GCNV(
        bams_channel,
        outdir_gatk,
        params.reference,
        genome_outdir,
        params.seqtk,         
        params.genome_chrs, 
        params.mappability_bed,
        params.segmental_duplication_bed, 
        scatter_count,
        sample_id_intervals_ch, 
        params.model_ploidy_outdir,
        params.model_cnvs_outdir,
        interval_ids,
        pedigree
    )

    // C: SURVIVOR
    SURVIVOR_MERGE(
        params.sample_id,
        CNMOPS.out.vcf,
        GATK_GCNV.out.genotyped_segments_vcf
    )

    // D: AnnotSV
    ANNOTSV(
        params.sample_id,
        SURVIVOR_MERGE.out.merged_vcf
    )
}



// Subworkflow: CNMOPS CNV Detection

// This subworkflow performs CNV detection using the cn.mops algorithm.
// It processes individual samples by normalizing against a cohort of control samples.
workflow CNMOPS {


    take:
        sample_id
        bam_file 
        bams_list      // val: path to file with list of all BAM files
        outdir


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



//
// Subworkflow: GATK gCNV Analysis
//
// This subworkflow performs germline CNV detection using GATK's gCNV caller.
// It includes genome preprocessing, read count collection, ploidy determination,
// CNV calling, postprocessing, and joint cohort segmentation.
workflow GATK_GCNV {


    take:
        bams_channel // channel of [sample_id, bam_file] tuples
        outdir            
        gr37_fasta_in            
        genome_outdir            
        seqtk                    
        genome_chrs              
        mappability_bed          
        segmental_duplication_bed 
        scatter_count
        sample_id_intervals_ch            
        model_ploidy_outdir
        model_cnvs_outdir
        interval_ids
        pedigree


    main:

        // Step 1: Preprocess genome fasta
        PREPROCESS_GENOME_FASTA(
            gr37_fasta_in,
            genome_outdir,
            seqtk,
            genome_chrs, 
            mappability_bed,
            segmental_duplication_bed
        )

        // Step 2: Collect read counts from samples
        COLLECT_READ_COUNTS(
            bams_channel,
            outdir,
            PREPROCESS_GENOME_FASTA.out.interval_list,
            PREPROCESS_GENOME_FASTA.out.ref_fasta,
            PREPROCESS_GENOME_FASTA.out.fasta_index,
            PREPROCESS_GENOME_FASTA.out.dict
        )

        // Step 3: Filter genome intervals based on read count outliers
        COLLECT_READ_COUNTS.out.sample_read_counts.collect().set { read_count_list }
        FILTER_GENOME(
            read_count_list,
            genome_outdir,
            PREPROCESS_GENOME_FASTA.out.annotated_interval_list,
            PREPROCESS_GENOME_FASTA.out.interval_list
        )

        // Step 4: Scatter genome intervals into shards
        SCATTER_GENOME(
            scatter_count,
            FILTER_GENOME.out.filtered_interval_list
        )

    
        // Step 5: Determine ploidy model - case mode
        DETERMINE_PLOIDY_CASE(
            bams_channel,
            COLLECT_READ_COUNTS.out.sample_read_counts,
            model_ploidy_outdir
        )

        // Step 6: Call germline CNVs in case mode
        CALL_CNVS_CASE(
            sample_id_intervals_ch,
            COLLECT_READ_COUNTS.out.sample_read_counts.first(),
            DETERMINE_PLOIDY_CASE.out.ploidy_calls.first(),
            scatter_count, 
            model_cnvs_outdir
        )
    
        // Step 7: Postprocess CNVs
        POSTPROCESS_CNVS(
            bams_channel,
            CALL_CNVS_CASE.out.cnv_calls_dir.collect(), 
            model_cnvs_outdir,
            PREPROCESS_GENOME_FASTA.out.dict,
            DETERMINE_PLOIDY_CASE.out.ploidy_calls,
            interval_ids.collect(),
            scatter_count
        )  

        /*
        // Step 8: Joint cohort segmentation
        JOINT_CNVS_SEGMENTATION(
            bams_channel,
            POSTPROCESS_CNVS.out.genotyped_segments_vcf.first(),
            POSTPROCESS_CNVS.out.genotyped_segments_vcf_index.first(),
            PREPROCESS_GENOME_FASTA.out.ref_fasta,
            PREPROCESS_GENOME_FASTA.out.fasta_index,
            PREPROCESS_GENOME_FASTA.out.dict,
            FILTER_GENOME.out.filtered_interval_list,
            pedigree
        )
        */
    

    emit:
        genome_fasta = PREPROCESS_GENOME_FASTA.out.ref_fasta
        read_counts = COLLECT_READ_COUNTS.out.sample_read_counts
        genotyped_segments_vcf = POSTPROCESS_CNVS.out.genotyped_segments_vcf
        genotyped_intervals_vcf = POSTPROCESS_CNVS.out.genotyped_intervals_vcf
        denoised_copy_ratios = POSTPROCESS_CNVS.out.denoised_copy_ratios
}












