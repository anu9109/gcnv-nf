process SURVIVOR_MERGE {

    publishDir "${params.outdir}/survivor", mode: 'copy'

    input:
        val sample_id
        path cnmops_vcf       // uncompressed VCF 
        path gatk_vcf         // compressed VCF 

    output:
        path "${sample_id}.cnvs.merged.vcf", emit: merged_vcf

    script:
    """
    # uncompress GATK VCF
    bcftools view \\
        --output-type v \\
        ${gatk_vcf} \\
        --threads 4 \\
        -o ${sample_id}_genotyped-segments.vcf

    # Create VCF list file for SURVIVOR
    echo "${cnmops_vcf}" > ${sample_id}.vcfs.txt
    echo "${sample_id}_genotyped-segments.vcf" >> ${sample_id}.vcfs.txt

    # Run SURVIVOR merge
    ${params.survivor_bin} merge \\
        ${sample_id}.vcfs.txt \\
        10000 \\
        1 \\
        1 \\
        1 \\
        0 \\
        50 \\
        ${sample_id}.cnvs.merged.vcf
    """
}


