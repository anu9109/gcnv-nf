process ANNOTSV {

    tag "AnnotSV on ${sample_id}"
    publishDir "${params.outdir}/annotsv", mode: 'copy'

    input:
        val sample_id
        path merged_vcf       // merged VCF from SURVIVOR_MERGE

    output:
        path "${sample_id}.cnvs.merged.annot.tsv", emit: annotated_tsv

    script:
    """
    #${params.annotsv} \\
    #    -SVinputFile ${merged_vcf} \\
    #    -SVinputInfo 1 \\
    #    -svtBEDcol 5 \\
    #    -genomeBuild GRCh37 \\
    #    -typeOfAnnotation both \\
    #    -outputDir . \\
    #    -outputFile ${sample_id}.cnvs.merged.annot.tsv \\
    #    -bedtools ${params.bedtools}

    # Run AnnotSV annotation
    ${params.annotsv_bin} \\
        --annotsv ${params.annotsv_dir} \\
        --bedtools ${params.bedtools} \\
        --outdir \${PWD} \\
        --slop 1000 \\
        --size_filter 40 \\
        --min_callers 1 \\
        --inp merged ${merged_vcf}
    """
}