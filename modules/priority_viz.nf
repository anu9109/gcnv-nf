process FILTER_PRIORITY_EVENTS {

    tag "Filter for priority events on ${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val sample_id
        path annotated_tsv

    output:
        path "${sample_id}.cnvs.merged.annot.priority.tsv", emit: priority_tsv

    script:
    """
    awk -F'\\t' '
        NR == 1 { print; next }
        {
            annotsv_type = \$18
            info         = \$12
            ranking      = \$72 + 0

            # Extract SUPP value from INFO field (e.g. SUPP=2)
            supp = 0
            if (match(info, /SUPP=([0-9]+)/, arr)) {
                supp = arr[1] + 0
            }

            if (annotsv_type == "full" && (ranking >= 4 || supp > 1)) {
                print
            }
        }
    ' ${annotated_tsv} > ${sample_id}.cnvs.merged.annot.priority.tsv
    """
}


process PLOT_EVENT_COVERAGE {

    tag "Coverage plots for priority events in ${sample_id}"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
        val sample_id
        path depth_file
        path priority_tsv

    output:
        path "*.pdf", emit: coverage_plots

    script:
    """
    plot_event_coverage.R ${sample_id} ${depth_file} ${priority_tsv}
    """
}
