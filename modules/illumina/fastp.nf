process fastp {
    container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'

    tag "trimming $sample"


    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_trimmed",
    pattern: "*trimmed_*fastq",
    mode: 'copy',
    overwrite: 'true'
    )

    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_report",
    pattern: "{*json, *html}",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("*trimmed_1.fastq"), path("*trimmed_2.fastq"), emit: trimmed
    tuple val(sample), path("*.json"), emit: fastP_json
    tuple val(sample), path("*.html"), emit: fastP_html


    script:
    """
    fastp \\
        --thread $params.thread \\
        --qualified_quality_phred $params.fastpPhred \\
        --in1 $fastq_1 \\
        --in2 $fastq_2 \\
        --out1 ${sample}.trimmed_1.fastq \\
        --out2 ${sample}.trimmed_2.fastq \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html
    """

}