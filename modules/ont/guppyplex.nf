process guppyplex {
    container 'staphb/artic:1.8.5'
    tag "🧹🧹🧹 ⌨️ GUPPYPLEX — cat *.fastq | grep \"quality\" | pray 🙏 (demultiplexing $sample) 🧹🧹🧹"


    // Publish merged FASTQ for downstream processes
    publishDir "${params.outDir}/guppyplex_fastq", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample), path(fastqDir)

    // Output a single FASTQ per barcode
    output:
    tuple val(sample), path("${sample}.fastq.gz"), emit: fastq

    script:
    """
    # Run artic guppyplex
    artic guppyplex \
        --min-length ${params.minLength} \
        --max-length ${params.maxLength} \
        --prefix ${sample} \
        --directory ${fastqDir} \
        &> /dev/null
    
    gzip -c ${sample}_${fastqDir}.fastq > ${sample}.fastq.gz
    """
}
