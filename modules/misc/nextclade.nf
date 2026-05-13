process nextclade {
    container 'nextstrain/nextclade:3.18.0'
    tag "🧬🧬🧬 Nextclade — Judging mutations 🧬🧬🧬"

    publishDir "${params.outDir}/nextclade", mode: 'copy'

    input:
    path consensus_fastas  // List of all *.consensus.fa

    output:
    path "nextclade.tsv", emit: tsv
    path "nextclade.csv", emit: csv
    path "allseq.fasta",  emit: combined

    script:
    """
    # Step 1: Concatenate all FASTA files
    cat $consensus_fastas > allseq.fasta

    # Step 2: Download dataset (cached)
    nextclade dataset get --name hMPXV --output-dir dataset/
    #  nextclade dataset get --name sars-cov-2 --output-dir dataset/

    # Step 3: Run nextclade on combined FASTA
    nextclade run \
        -D dataset/ \
        -t nextclade.tsv \
        --output-csv nextclade.csv \
        allseq.fasta
    """
}
