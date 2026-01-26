process kraken {
    container 'staphb/kraken2:2.1.2-no-db'

    tag "check $sample"


    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_report",
    pattern: "*kraken2.report.txt",
    mode: "copy",
    overwrite: "true"
    )

    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_fastq",
    pattern: "*unclassified.fastq",
    mode: "copy",
    overwrite: "true"
    )

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    path(krakenDb)

    output:
    tuple val(sample), path("*1.unclassified.fastq"), path("*2.unclassified.fastq"), emit: dehosted
    tuple val(sample), path("*.kraken2.report.txt"), emit: taxon


    script:
    """
    kraken2 \\
        --threads $params.thread \\
        --memory-mapping \\
        --use-names \\
        --paired $fastq_1 $fastq_2 \\
        --classified-out ${sample}#.classified.fastq \\
        --unclassified-out ${sample}#.unclassified.fastq \\
        --db $krakenDb \\
        --report ${sample}.kraken2.report.txt
    """
}