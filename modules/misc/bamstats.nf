#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process bamstats {
    container 'wewantsaul/seqkit-pysam:latest' 
    tag "📊📊📊 BAMstats — depth so low even Kraken2 looks confident 😌 📊📊📊"


    publishDir "${params.outDir}/bamstats", mode: 'copy', overwrite: true

    input:
    path bams
    path nextclade_tsv

    output:
    path "*_depth.png",              emit: depth_plots
    path "average_depth.png",        emit: avg_plot
    path "average_depth_summary.tsv",emit: depth_summary
    path "summary.txt",              emit: final_summary

    script:
    """
    set -euo pipefail

    mkdir -p bam_dir

    # ------------------------------------------------------------
    # 1. Link BAMs + BAI (already exist from artic)
    # ------------------------------------------------------------
    for file in $bams; do
        cp -L "\$file" bam_dir/
    done

    # ------------------------------------------------------------
    # 2. Generate depth plots (bamplotter.py must be in PATH or bin/)
    # ------------------------------------------------------------
    bamplotter.py bam_dir -o . > /dev/null

    # ------------------------------------------------------------
    # 3. Merge depth + Nextclade
    # ------------------------------------------------------------
    merge_summary.py \
        --depth average_depth_summary.tsv \
        --nextclade $nextclade_tsv \
        --output summary.txt > /dev/null
    """
}
