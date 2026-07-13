
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
    publishDir "${params.outdir}", mode: 'copy' 

    input:
        tuple val(sample_id), val(bam_file)
        val outdir
        path interval_list
        path ref_fasta
        path fasta_index
        path ref_dict

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

process SCATTER_GENOME {

    tag "Scatter genome intervals into ${scatter_count} shards"
    publishDir "${params.outdir}/scatter", mode: 'copy'

    input:
        val scatter_count
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

process DETERMINE_PLOIDY_CASE {
    
    tag "Determine germline contig ploidy on sample ${sample_id}"
    errorStrategy 'retry'   
    maxRetries 3
    memory {
        // Base memory (8 GB), doubled for each retry attempt
        def base = 8.GB
        return base * Math.pow(2, task.attempt - 1)
    }

    publishDir "${params.outdir}", mode: 'copy'

    input: 
        tuple val(sample_id), val(bam_file)
        path sample_read_counts
        val model_ploidy_outdir

    output:
        path "${sample_id}_ploidy-calls", emit: ploidy_calls

    script:
    """
    export PYTENSOR_FLAGS="base_compiledir=\$PWD/.pytensor"
    mkdir -p \$PWD/.pytensor
    gatk DetermineGermlineContigPloidy \\
        --model ${model_ploidy_outdir}-model \\
        -I ${sample_read_counts} \\
        -O . \\
        --output-prefix ${sample_id}_ploidy
    """
}


process CALL_CNVS_CASE {

    tag "Calling CNVs on sample ${sample_id}"
    errorStrategy 'retry'
    maxRetries 2
    memory {
        def base = 32.GB
        return base * Math.pow(2, task.attempt - 1)
    }

    publishDir "${params.outdir}/cnvs", mode: 'copy'

    input:
        tuple val(sample_id), val(bam_file), val (interval_id)
        path sample_read_counts
        path ploidy_calls
        val scatter_count
        val model_cnvs_outdir

    output:
        path "${sample_id}_cnv_${interval_id}_of_${scatter_count}-calls", emit: cnv_calls_dir
        path "${sample_id}_cnv_${interval_id}_of_${scatter_count}-model", emit: cnv_model_dir

    script:
    """
    gatk GermlineCNVCaller \\
        --run-mode CASE \\
        --input ${sample_read_counts} \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --model ${model_cnvs_outdir}_${interval_id}_of_${scatter_count}-model \\
        --output . \\
        --output-prefix ${sample_id}_cnv_${interval_id}_of_${scatter_count}
    """
}


process POSTPROCESS_CNVS {

    tag "Postprocess CNVs for sample ${sample_id}"
    publishDir "${params.outdir}/cnvs", mode: 'copy'

    input:
        tuple val(sample_id), val(bam_file)
        path cnv_calls_dir
        path cnv_model_dir
        path dict 
        path ploidy_calls

    output:
        path "${sample_id}_genotyped-intervals.vcf.gz", emit: genotyped_intervals
        path "${sample_id}_genotyped-segments.vcf.gz", emit: genotyped_segments
        path "${sample_id}_denoised_copy_ratios.tsv", emit: denoised_copy_ratios

    script:
    def calls_list = cnv_calls_dir instanceof List ? cnv_calls_dir : [cnv_calls_dir]
    def models_list = cnv_model_dir instanceof List ? cnv_model_dir : [cnv_model_dir]
    def calls_shard_args = calls_list.sort().collect { "--calls-shard-path ${it}" }.join(' ')
    def model_shard_args = models_list.sort().collect { "--model-shard-path ${it}" }.join(' ')
    """
    echo "Postprocessing CNVs for sample ${sample_id}"

    gatk PostprocessGermlineCNVCalls \\
        --sample-index 0 \\
        --allosomal-contig X \\
        --allosomal-contig Y \\
        --contig-ploidy-calls ${ploidy_calls} \\
        --output-genotyped-intervals ${sample_id}_genotyped-intervals.vcf.gz \\
        --output-genotyped-segments ${sample_id}_genotyped-segments.vcf.gz \\
        --output-denoised-copy-ratios ${sample_id}_denoised_copy_ratios.tsv \\
        --sequence-dictionary ${dict} \\
        ${model_shard_args} \\
        ${calls_shard_args}
    """
}

process JOINT_CNVS_SEGMENTATION {

    tag "Joint CNV segmentation for sample ${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input: 
        tuple val(sample_id), val(bam_file)
        path segment_vcf
        path ref_fasta
        path fasta_index
        path filtered_interval_list
        val pedigree

    output: 
        path "${sample_id}_clustered.vcf.gz", emit: clustered_vcf

    script:
    """
    gatk JointGermlineCNVSegmentation \\
        -R ${ref_fasta} \\
        -V ${segment_vcf} \\
        --model-call-intervals ${filtered_interval_list} \\
        --pedigree ${pedigree} \\
        -O ${sample_id}_clustered.vcf.gz
    """
}