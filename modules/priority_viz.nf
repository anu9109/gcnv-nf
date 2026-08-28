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
    cat > filter.awk << 'AWKEOF'
    BEGIN { FS = "\t"; OFS = "\t" }
    NR == 1                                                    { print; next }
    \$16 == "full" && tolower(\$53) ~ /pathogenic/ { print }
    AWKEOF
    awk -f filter.awk ${annotated_tsv} > ${sample_id}.cnvs.merged.annot.priority.tsv
    """
}


process PLOT_EVENT_COVERAGE {

    tag "Coverage plots for priority events in ${sample_id}"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
        val sample_id
        path depth_file
        path priority_tsv
        path gc_file
        path map_file

    output:
        path "*.png", emit: coverage_plots

    script:
    """
    plot_event_coverage.R ${sample_id} ${depth_file} ${priority_tsv} ${gc_file} ${map_file} ${params.genome_build}
    """
}