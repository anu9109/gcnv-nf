
process PREPROCESS_GENOME_FASTA {

    publishDir "${params.genome_outdir}", mode: 'copy' 

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
        -R ${fasta} \\
        -imr OVERLAPPING_ONLY \\
        -I ${bam_file} \\
        --format TSV \\
        -O ${sample_id}.read_counts.tsv
    
    echo "Read counts collected: ${sample_id}.read_counts.tsv"
    """
}

process FILTER_GENOME {

    tag "Filter genome intervals"
    publishDir "${genome_outdir}", mode: 'copy'

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


process SCATTER_GENOME {

    tag "Scatter genome intervals into ${scatter_count} shards"
    publishDir "${scatter_outdir}", mode: 'copy'

    input:
        val scatter_count
        val scatter_outdir
        path filtered_interval_list

    output:
        path "temp_*", emit: scattered_intervals, type: 'dir'

    script:
    """
    echo "Scattering genome intervals into ${scatter_count} shards"

    gatk IntervalListTools \\
        --INPUT ${filtered_interval_list} \\
        --SUBDIVISION_MODE INTERVAL_COUNT \\
        --SCATTER_COUNT ${scatter_count} \\
        --OUTPUT .

    echo "Scattered intervals created"
    ls -la temp_*/
    """
}


process DETERMINE_PLOIDY {
    
    tag "Determine germline contig ploidy"
    errorStrategy 'retry'
    maxRetries 3
    memory {
        // Base memory (8 GB), doubled for each retry attempt
        def base = 8.GB
        return base * Math.pow(2, task.attempt - 1)
    }
    publishDir "${model_outdir}", mode: 'copy'

    input: 
        path(read_count_files)
        val model_outdir
        val model_prefix
        path filtered_interval_list
        path ploidy_priors

    output:
        path "${model_prefix}-calls", emit: ploidy_calls

    script:
    def readcount_args = read_count_files.collect { "-I ${it}" }.join(' ')
    """
    echo "Determining germline contig ploidy from ${read_count_files.size()} samples"

    gatk DetermineGermlineContigPloidy \\
        -L ${filtered_interval_list} \\
        -imr OVERLAPPING_ONLY \\
        ${readcount_args} \\
        --contig-ploidy-priors ${ploidy_priors} \\
        --output . \\
        --output-prefix ${model_prefix}

    echo "Ploidy model created"
    ls -la ${model_prefix}-*
    """
}

process DETERMINE_PLOIDY_CASE {
    
    errorStrategy 'retry'
    maxRetries 3
    memory {
        // Base memory (8 GB), doubled for each retry attempt
        def base = 8.GB
        return base * Math.pow(2, task.attempt - 1)
    }

    input: 
    tuple val(sample_id), val(bam_file)

    script:
    """
    gatk DetermineGermlineContigPloidy \\
        --model ${params.model_outdir}/${params.model_prefix}-model\\
        -I ${params.sample_outdir}/${sample_id}.read_counts.tsv \\
        -O ${params.sample_outdir} \\
        --output-prefix ${sample_id}_ploidy
    """
}

process CALL_CNVS_CASE {

    errorStrategy 'retry'
    maxRetries 2
    memory {
        def base = 32.GB
        return base * Math.pow(2, task.attempt - 1)
    }

    tag "Calling CNVs on sample ${sample_id}"

    input:
    tuple(val(sample_id), val(bam_file), val (interval_id))

    script:
    """
    mkdir -p ${params.sample_outdir}/cnvs
    gatk GermlineCNVCaller \\
        --run-mode CASE \\
        --contig-ploidy-calls ${params.model_outdir}/${params.model_prefix}-calls \\
        --model ${params.cnvs_outdir}/${params.cnvs_prefix}_${interval_id}_of${params.scatter_count}-model \\
        --input ${params.sample_outdir}/${sample_id}.read_counts.tsv \\
        --output ${params.sample_outdir}/cnvs \\
        --output-prefix ${sample_id}_case_cnvs_${interval_id}_of${params.scatter_count}
    """
}


process POSTPROCESS_CNVS {

    tag "Postprocess CNVs for sample ${sample_id}"
    publishDir "${cnvs_outdir}", mode: 'copy'

    input:
        tuple val(sample_index), val(sample_id)
        val interval_ids
        val cnvs_outdir
        val cnvs_prefix
        val scatter_count
        val model_outdir
        val model_prefix
        val genome_outdir
        val pedigree

    output:
        path "${sample_id}_genotyped-intervals.vcf.gz", emit: genotyped_intervals
        path "${sample_id}_denoised_copy_ratios.tsv", emit: denoised_copy_ratios

    script:
    def model_shard_args = interval_ids.collect { "--model-shard-path ${cnvs_outdir}/${cnvs_prefix}_${it}_of_${scatter_count}-model" }.join(' ')
    def calls_shard_args = interval_ids.collect { "--calls-shard-path ${cnvs_outdir}/${cnvs_prefix}_${it}_of_${scatter_count}-calls" }.join(' ')
    """
    echo "Postprocessing CNVs for sample ${sample_id} (index ${sample_index})"

    gatk postprocessGermlineCNVCalls \\
        --sample-index ${sample_index} \\
        --allosomal-contig X \\
        --allosomal-contig Y \\
        --contig-ploidy-calls ${model_outdir}/${model_prefix}-calls \\
        --output-genotyped-intervals ${sample_id}_genotyped-intervals.vcf.gz \\
        --output-genotyped-segments ${sample_id}_denoised_copy_ratios.tsv \\
        --sequence-dictionary ${genome_outdir}/gr37_clean.dict \\
        ${model_shard_args} \\
        ${calls_shard_args}

    echo "Postprocessing completed for ${sample_id}"
    """
}





