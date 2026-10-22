process ANNOTSV {

    tag "AnnotSV on ${sample_id}"
    publishDir "${params.outdir}/annotsv", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3
    memory {
        // AnnotSV loads its annotation databases into memory
        def base = 16.GB
        return base * Math.pow(2, task.attempt - 1)
    }

    input:
        val sample_id
        path merged_vcf       // merged VCF from SURVIVOR_MERGE

    output:
        path "${sample_id}.cnvs.merged.annot.tsv", emit: annotated_tsv

    script:
    """
    export ANNOTSV="${params.annotsv_dir}"
    ${params.annotsv_bin} \\
        -SVinputFile ${merged_vcf} \\
        -SVinputInfo 1 \\
        -svtBEDcol 5 \\
        -genomeBuild GRCh37 \\
        -typeOfAnnotation both \\
        -outputDir . \\
        -outputFile ${sample_id}.cnvs.merged.annot.tsv \\
        -bedtools ${params.bedtools}
    """
}