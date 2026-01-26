// enable dsl2
nextflow.enable.dsl=2


// import modules
include {fastp} from '../modules/illumina/fastp.nf'
include {fastqcRawPE} from '../modules/misc/fastqc.nf'
include {fastqcTrimmedPE} from '../modules/misc/fastqc.nf'
include {minimap2PE} from '../modules/misc/minimap2.nf'
include {sam2bam} from '../modules/misc/samtools.nf'
include {sortIndexMinimap} from '../modules/misc/samtools.nf'
include {lofreqVariant} from '../modules/illumina/variantCall.nf'
include {kraken} from '../modules/misc/dehost.nf'


workflow illuminaShotgun {
    
    main:
        Channel
            .fromFilePairs("${params.inDir}/*{,.trimmed}_{R1,R2,1,2}{,_001}.{fastq,fq}{,.gz}", flat:true)
            .ifEmpty{error "Cannot find any reads matching: ${params.inDir}"}
            .set{ch_sample}
            
        fastqcRawPE(ch_sample)
        fastp(ch_sample)
        fastqcTrimmedPE(fastp.out.trimmed)

        kraken(fastp.out.trimmed, params.krakenDb)
        minimap2PE(kraken.out.dehosted, params.refGenome)
        sam2bam(minimap2PE.out.sam)
        sortIndexMinimap(sam2bam.out.bam)

        lofreqVariant(sortIndexMinimap.out.bamBai, params.refGenome)

}