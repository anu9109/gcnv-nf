process GENERATE_REPORT {

    tag "CNV report for ${sample_id}"
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
        val  sample_id
        path cnmops_vcf
        path gatk_vcf
        path merged_vcf
        path annotated_tsv
        path priority_tsv
        path coverage_plots   // collection of PDFs (may be empty)

    output:
        path "${sample_id}.cnv_report.html", emit: report_html

    script:
    """
    create_report.py \\
        ${sample_id} \\
        ${cnmops_vcf} \\
        ${gatk_vcf} \\
        ${merged_vcf} \\
        ${annotated_tsv} \\
        ${priority_tsv} \\
        ${coverage_plots}
    """
}
