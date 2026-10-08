#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPARE_SAMPLE_LIST } from './modules/cnmops.nf'
include { RUN_CNMOPS } from './modules/cnmops.nf'
include { CNMOPS_TO_VCF } from './modules/cnmops.nf'
include { SUBSET_GENOME_FASTA } from './modules/gatk.nf'
include { PREPROCESS_GENOME_FASTA } from './modules/gatk.nf'
include { COLLECT_READ_COUNTS } from './modules/gatk.nf'
include { FILTER_GENOME } from './modules/gatk.nf'
include { SCATTER_GENOME } from './modules/gatk.nf'
include { DETERMINE_PLOIDY_CASE } from './modules/gatk.nf'
include { CALL_CNVS_CASE } from './modules/gatk.nf'
include { POSTPROCESS_CNVS } from './modules/gatk.nf'
include { FILTER_GATK } from './modules/gatk.nf'
//include { JOINT_CNVS_SEGMENTATION } from './modules/gatk.nf'
include { SURVIVOR_MERGE } from './modules/survivor.nf'
include { ANNOTSV } from './modules/annotsv.nf'
include { FILTER_PRIORITY_EVENTS } from './modules/priority_viz.nf'
include { PLOT_EVENT_COVERAGE } from './modules/priority_viz.nf'
include { GENERATE_REPORT } from './modules/report.nf'


workflow {
    if (params.get('rerun_viz', false)) {
        RERUN_VIZ()
    } else {
    
        // A: CNMOPS
        CNMOPS(
            params.sample_id, 
            file(params.bam_file),
            file(params.bams_list)
        )
        
        // B: GATK_GCNV
        bams_channel = Channel.of(tuple(params.sample_id, file(params.bam_file)))
        scatter_count = params.scatter_count as int
        interval_ids = Channel
            .from(1..scatter_count)
            .map { String.format("%04d", it) }
        bams_channel
            .combine(interval_ids)
            .map { row -> tuple(row[0], row[1], row[2]) }
            .set { sample_id_intervals_ch }
        pedigree = file("${params.outdir}/gatk_gcnv/pedigree.txt")
        GATK_GCNV(
            bams_channel,
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
            GATK_GCNV.out.genotyped_segments_filtered_vcf
        )

        // D: AnnotSV
        ANNOTSV(
            params.sample_id,
            SURVIVOR_MERGE.out.merged_vcf
        )

        // E: Filter AnnotSV results for pathogenic/likely pathogenic events and for events supported by more than one caller
        FILTER_PRIORITY_EVENTS(
            params.sample_id,
            ANNOTSV.out.annotated_tsv
        )

        // F: Plot read-depth coverage for each priority event
        PLOT_EVENT_COVERAGE(
            params.sample_id,
            file(params.depth_file),
            FILTER_PRIORITY_EVENTS.out.priority_tsv,
            file(params.gc_file),
            file(params.map_file)
        )

        // G: Generate HTML report
        GENERATE_REPORT(
            params.sample_id,
            CNMOPS.out.vcf,
            GATK_GCNV.out.genotyped_segments_filtered_vcf,
            SURVIVOR_MERGE.out.merged_vcf,
            ANNOTSV.out.annotated_tsv,
            FILTER_PRIORITY_EVENTS.out.priority_tsv,
            PLOT_EVENT_COVERAGE.out.coverage_plots.collect()
        )
    }
}

workflow RERUN_VIZ {

    // Point these at your existing output files
    depth_file = file(params.depth_file)
    cnmops_vcf = file("${params.outdir}/cnmops/${params.sample_id}.vcf")
    gatk_vcf   = file("${params.outdir}/gatk_gcnv/${params.sample_id}_genotyped-segments-filtered.vcf.gz")

    // C: SURVIVOR
    SURVIVOR_MERGE(
        params.sample_id,
        cnmops_vcf,
        gatk_vcf
    )

    // D: AnnotSV
    ANNOTSV(
        params.sample_id,
        SURVIVOR_MERGE.out.merged_vcf
    )

    // E: Filter priority events
    FILTER_PRIORITY_EVENTS(
        params.sample_id,
        ANNOTSV.out.annotated_tsv
    )

    // F: Plot coverage
    PLOT_EVENT_COVERAGE(
        params.sample_id,
        depth_file,
        FILTER_PRIORITY_EVENTS.out.priority_tsv,
        file(params.gc_file),
        file(params.map_file)
    )

    // G: Generate HTML report
    GENERATE_REPORT(
        params.sample_id,
        cnmops_vcf,
        gatk_vcf,
        SURVIVOR_MERGE.out.merged_vcf,
        ANNOTSV.out.annotated_tsv,
        FILTER_PRIORITY_EVENTS.out.priority_tsv,
        PLOT_EVENT_COVERAGE.out.coverage_plots.collect()
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


    main:
        // Step 1: Prepare sample list with controls
        PREPARE_SAMPLE_LIST(
            sample_id,
            bam_file,
            bams_list
        )

        // Step 2: Run cn.mops analysis
        RUN_CNMOPS(
            PREPARE_SAMPLE_LIST.out.bams_txt,
            sample_id
        )

        // Step 3: Convert results to VCF format
        CNMOPS_TO_VCF(
            sample_id,
            RUN_CNMOPS.out.sample_cnvs
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
// It includes read count collection, ploidy determination,
// CNV calling, postprocessing, and joint cohort segmentation.
// Expects a pre-built genome (see standalone PREPARE_GENOME workflow) at params.genome_path.
workflow GATK_GCNV {


    take:
        bams_channel // channel of [sample_id, bam_file] tuples        
        scatter_count
        sample_id_intervals_ch            
        model_ploidy_outdir
        model_cnvs_outdir
        interval_ids
        pedigree


    main:

        // Pre-built genome files produced by the standalone PREPARE_GENOME workflow
        ref_fasta               = file("${params.ref_path}/gr37_clean.fasta")
        fasta_index             = file("${params.ref_path}/gr37_clean.fasta.fai")
        dict                    = file("${params.ref_path}/gr37_clean.dict")
        interval_list           = file("${params.ref_path}/gr37_clean.interval_list")
        annotated_interval_list = file("${params.ref_path}/gr37_clean_annotated.interval_list")

        // Step 1: Collect read counts from samples
        COLLECT_READ_COUNTS(
            bams_channel,
            interval_list,
            ref_fasta,
            fasta_index,
            dict
        )

        // Step 2: Filter genome intervals based on read count outliers
        COLLECT_READ_COUNTS.out.sample_read_counts.collect().set { read_count_list }
        FILTER_GENOME(
            read_count_list,
            annotated_interval_list,
            interval_list
        )

        // Step 3: Scatter genome intervals into shards
        SCATTER_GENOME(
            scatter_count,
            FILTER_GENOME.out.filtered_interval_list
        )

    
        // Step 4: Determine ploidy model - case mode
        DETERMINE_PLOIDY_CASE(
            bams_channel,
            COLLECT_READ_COUNTS.out.sample_read_counts,
            model_ploidy_outdir
        )

        // Step 5: Call germline CNVs in case mode
        CALL_CNVS_CASE(
            sample_id_intervals_ch,
            COLLECT_READ_COUNTS.out.sample_read_counts.first(),
            DETERMINE_PLOIDY_CASE.out.ploidy_calls.first(),
            scatter_count, 
            model_cnvs_outdir
        )
    
        // Step 6: Postprocess CNVs
        POSTPROCESS_CNVS(
            bams_channel,
            CALL_CNVS_CASE.out.cnv_calls_dir.collect(), 
            model_cnvs_outdir,
            dict,
            DETERMINE_PLOIDY_CASE.out.ploidy_calls,
            interval_ids.collect(),
            scatter_count
        )  

        // Step 7: Filter out diploid (ALT=.) segments
        FILTER_GATK(
            bams_channel.map { sample_id, bam -> sample_id },
            POSTPROCESS_CNVS.out.genotyped_segments_vcf,
            POSTPROCESS_CNVS.out.genotyped_segments_vcf_index
        )

        /*
        // Step 8: Joint cohort segmentation
        JOINT_CNVS_SEGMENTATION(
            bams_channel,
            POSTPROCESS_CNVS.out.genotyped_segments_vcf.first(),
            POSTPROCESS_CNVS.out.genotyped_segments_vcf_index.first(),
            ref_fasta,
            fasta_index,
            dict,
            FILTER_GENOME.out.filtered_interval_list,
            pedigree
        )
        */
    

    emit:
        genome_fasta = ref_fasta
        read_counts = COLLECT_READ_COUNTS.out.sample_read_counts
        genotyped_segments_vcf = POSTPROCESS_CNVS.out.genotyped_segments_vcf
        genotyped_intervals_vcf = POSTPROCESS_CNVS.out.genotyped_intervals_vcf
        denoised_copy_ratios = POSTPROCESS_CNVS.out.denoised_copy_ratios
        genotyped_segments_filtered_vcf = FILTER_GATK.out.genotyped_segments_filtered_vcf
}


//
// Standalone workflow: Genome Preparation
//
// Run independently from the main pipeline to preprocess a new genome fasta into the
// files expected at params.genome_path (ref fasta, index, dict, interval lists).
// Usage: nextflow run main.nf -entry PREPARE_GENOME --genome_fasta /path/to/genome.fasta
workflow PREPARE_GENOME {
    SUBSET_GENOME_FASTA(
        file(params.genome_fasta)
    )
    PREPROCESS_GENOME_FASTA(
        SUBSET_GENOME_FASTA.out.subset_fasta
    )
}


workflow.onComplete {
    if (workflow.success) {
        def workDir = new File(workflow.workDir.toString())
        log.info "Pipeline completed successfully. Removing work directory: ${workDir}"
        workDir.deleteDir()
    }
}










