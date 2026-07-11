
process PREPROCESS_GENOME_FASTA {

    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
        val gr37_fasta_in
        val genome_outdir
        val seqtk
        val genome_chrs
        val mappability_bed
        val segmental_duplication_bed

    output:
        path "gr37_clean.fasta", emit: ref_fasta
        path "gr37_clean.fasta.fai", emit: fasta_index
        path "gr37_clean.dict", emit: dict
        path "gr37_clean.interval_list", emit: interval_list
        path "gr37_clean_annotated.interval_list", emit: annotated_interval_list

    script:
    """
    echo "1. Subsetting genome fasta to only chromosomes 1-22, X and Y.."
    ${seqtk} subseq ${gr37_fasta_in} ${genome_chrs} > gr37_clean.fasta
    
    echo "2. Indexing genome fasta.."
    samtools faidx gr37_clean.fasta

    echo "3. Creating genome dictionary file.."
    gatk CreateSequenceDictionary -R gr37_clean.fasta

    echo "4. Preprocessing genome fasta into GATK .interval_list format.."
    gatk PreprocessIntervals -R gr37_clean.fasta \\
        --padding 0 \\
        -imr OVERLAPPING_ONLY \\
        -O gr37_clean.interval_list

    echo "5. Annotating intervals with GC %, mappability and segmental duplication content.."
    gatk AnnotateIntervals -L gr37_clean.interval_list \\
        -R gr37_clean.fasta \\
        --mappability-track ${mappability_bed} \\
        --segmental-duplication-track ${segmental_duplication_bed} \\
        -imr OVERLAPPING_ONLY \\
        -O gr37_clean_annotated.interval_list
    """
}


process COLLECT_READ_COUNTS {

    tag "Collect read counts from ${sample_id}"
    publishDir "${params.sample_outdir}", mode: 'copy' 

    input:
        tuple val(sample_id), val(bam_file)
        val sample_outdir
        path interval_list
        path ref_fasta

    output:
        path "*.read_counts.tsv", emit: sample_read_counts

    script:
    """
    gatk CollectReadCounts \\
        -L ${interval_list} \\
        -R ${ref_fasta} \\
        -imr OVERLAPPING_ONLY \\
        -I ${bam_file} \\
        --format TSV \\
        -O ${sample_id}.read_counts.tsv
    
    echo "Read counts collected: ${sample_id}.read_counts.tsv"
    """
}

process FILTER_GENOME {

    tag "Filter genome intervals"
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
        path(read_count_files)
        val genome_outdir
        path annotated_interval_list
        path interval_list

    output:
        path "grch37_clean.annotated.filtered.interval_list", emit: filtered_interval_list

    script:
    def readcount_args = read_count_files.collect { "-I ${it}" }.join(' ')
    """
    echo "Filtering genome intervals based on sample read count distribution.."

    gatk FilterIntervals \\
        -L ${interval_list} \\
        --annotated-intervals ${annotated_interval_list} \\
        -imr OVERLAPPING_ONLY \\
        -O grch37_clean.annotated.filtered.interval_list \\
        ${readcount_args}

    echo "Sample specific filtered interval list created."
    """
}

