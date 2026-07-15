process PREPARE_SAMPLE_LIST {

    publishDir "${params.outdir}/cnmops", mode: 'copy'

    input:
        val sample_id
        val bam_file
        val bams_list

    output:
        path "bams.txt", emit: bams_txt

    script:
    """
    echo "1. Creating output directory if it does not exist"
    mkdir -p ${params.outdir}/cnmops

    echo "2. Gathering samples to normalize sample of interest against"
    shuf -n 19 ${bams_list} > bams.txt
    echo -e "${sample_id}\\t${bam_file}" >> bams.txt
    """
}

process RUN_CNMOPS {

    tag "cn.mops on ${sample_id}"
    publishDir "${params.outdir}/cnmops", mode: 'copy'

    input:
        path bams_txt
        val sample_id

    output:
        path "${sample_id}_cohort_segmentation.tsv", emit: sample_seg
        path "${sample_id}_cohort_cnvs.tsv", emit: sample_cnvs
        path "${sample_id}_cohort_cnvr.tsv", emit: sample_cnvr

    script:
    """
    echo "3. Running cn.mops"
    run_cnmops.R ${bams_txt} 8 ${sample_id}
    """
}

process CNMOPS_TO_VCF {

    publishDir "${params.outdir}/cnmops", mode: 'copy'
    
    input:
        val sample_id
        path sample_cnvs

    output:
        path "${sample_id}.vcf", emit: sample_cnvs_vcf

    script:
    """
    echo "4. Create cn.mops VCF"
    cnmops_to_vcf_per_sample.R ${sample_id} ${sample_cnvs}
    """
}

