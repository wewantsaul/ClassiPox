process lofreqVariant {
    container 'ufuomababatunde/lofreq:v2.1.5'

    tag "variants and consensus: $sample"


    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_variant",
    pattern: "*.norm.filter.vcf",
    mode: 'copy',
    overwrite: 'true'
    )

    publishDir (
    path: "${params.outDir}/${task.process.replaceAll(":","_")}_consensus",
    pattern: "*.consensus.fasta",
    mode: 'copy',
    overwrite: 'true'
    )



    input:
    tuple val(sample), path(bam), path(bam_bai)
    path(reference)

    output:
    tuple val(sample), path("*.norm.filter.vcf"), emit: vcf
    tuple val(sample), path("*consensus.fasta"), emit: consensus


    script:
    """
    # Insert qualities of the indels into the BAM files
    lofreq indelqual --dindel \\
        -f ${reference} \\
        -o ${sample}.indelqual.bam \\
        $bam

    # Index the new BAM files
    samtools index ${sample}.indelqual.bam


    # Call indels aside from other variants
    lofreq call --call-indels \\
        -f ${reference} \\
        -o ${sample}.vcf \\
        ${sample}.indelqual.bam

    bgzip -c ${sample}.vcf > ${sample}.vcf.gz
    tabix -p vcf ${sample}.vcf.gz


    # Filter called variants with at least 10x depth and 20% frequency
    bcftools view \\
        -i 'DP>=10 & AF>=0.2' \\
        ${sample}.vcf.gz > ${sample}.pass.vcf

    bgzip -c ${sample}.pass.vcf > ${sample}.pass.vcf.gz
    tabix -p vcf ${sample}.pass.vcf.gz


    # Normalize indels
    bcftools norm \\
        -f ${reference}\\
        ${sample}.pass.vcf.gz -Ob \\
        -o ${sample}.norm.bcf

    # filter/remove clusters of indels separated by 50 base pairs thereby allowing only one to pass
    bcftools filter \\
        --IndelGap 50 \
        ${sample}.norm.bcf -Ob \\
        -o ${sample}.norm.filter.bcf

    tabix ${sample}.norm.filter.bcf

    bcftools convert \\
        -O v \\
        -o ${sample}.norm.filter.vcf \\
        ${sample}.norm.filter.bcf

    bgzip -c ${sample}.norm.filter.vcf > ${sample}.norm.filter.vcf.gz
    tabix -p vcf ${sample}.norm.filter.vcf.gz


    # Create consensus sequence with mix sites represented by their equivalent IUPAC codes (e.g. A or G as R) 
    cat ${reference}| bcftools consensus --iupac-codes ${sample}.norm.filter.vcf.gz > ${sample}.consensus.fasta

    sed -i '1s/.*/>${sample}/' ${sample}.consensus.fasta
    """
}