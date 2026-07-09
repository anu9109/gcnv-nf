process PREPARE_SAMPLE_LIST {

    publishDir "${params.outdir}", mode: 'copy'

    input:
        val sample_id
        val bam_file
        val outdir
        val bams_list

    output:
        path "bams.txt", emit: bams_txt

    script:
    """
    echo "1. Creating output directory if it does not exist"
    mkdir -p ${outdir}

    echo "2. Gathering samples to normalize sample of interest against"
    shuf -n 19 ${bams_list} > bams.txt
    echo -e "${sample_id}\\t${bam_file}" >> bams.txt
    """
}

process RUN_CNMOPS {

    tag "cn.mops on ${params.sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path bams_txt
        //val cnmops_sif
        //val run_cnmops_r
        val sample_id
        val outdir

    output:
        path "${sample_id}_cohort_segmentation.tsv", emit: sample_seg
        path "${sample_id}_cohort_cnvs.tsv", emit: sample_cnvs
        path "${sample_id}_cohort_cnvr.tsv", emit: sample_cnvr

    script:
    """
    echo "3. Running cn.mops"
    run_cnmops.R ${bams_txt} 32 ${sample_id}
    """
}

process CNMOPS_TO_VCF {

    publishDir "${params.outdir}", mode: 'copy'
    
    input:
        //val cnmops_sif
        //val cnmops_to_vcf_r
        val sample_id
        path sample_cnvs
        val outdir

    output:
        path "${sample_id}.vcf", emit: sample_cnvs_vcf

    script:
    """
    echo "4. Create cn.mops VCF"
    cnmops_to_vcf_per_sample.R ${sample_id} ${sample_cnvs}
    """
}

