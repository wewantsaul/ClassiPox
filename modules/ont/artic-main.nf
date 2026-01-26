process artic {
    container 'staphb/artic:1.8.5'
    tag "🦠🦠🦠 ARTIC — Stretching those amplicons 😏 (processing $sample) 🦠🦠🦠"

    publishDir "${params.outDir}/artic_consensus", mode: 'copy', pattern: '*.consensus.fasta'
    publishDir "${params.outDir}/artic_bam",       mode: 'copy', pattern: '*.primertrimmed.rg.sorted.bam*'

    input:
    tuple val(sample), path(guppyplex_fastq)
    path bed
    path reference
    path model

    output:
    tuple val(sample), path("${sample}.consensus.fasta"), emit: consensus
    tuple val(sample), 
          path("${sample}.primertrimmed.rg.sorted.bam"), 
          path("${sample}.primertrimmed.rg.sorted.bam.bai"), 
          emit: sortedBam
        
    script:
    """
    artic minion \
        --model-dir $model \
        --model r1041_e82_400bps_sup_v430 \
        --normalise ${params.normalise} \
        --threads 32 \
        --bed $bed \
        --ref $reference \
        --read-file $guppyplex_fastq \
        $sample \
        &> /dev/null
    """
}
