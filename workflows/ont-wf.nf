nextflow.enable.dsl = 2

include { guppyplex } from '../modules/ont/guppyplex.nf'
include { artic } from '../modules/ont/artic-main.nf'
include { nextclade } from '../modules/misc/nextclade.nf'
include { bamstats } from '../modules/misc/bamstats.nf'


workflow nanoporeArtic {

    take:
        samplesheet_map // Map<String, String> or null

    main:
        // ========== 0. CHECK IF SAMPLE SHEET EXISTS ==========
        def hasSamplesheet = samplesheet_map && !samplesheet_map.isEmpty()

        if (hasSamplesheet) {
            log.info "Using provided CSV samplesheet for barcode mapping (${samplesheet_map.size()} samples)."
            log.info "Will only process: ${samplesheet_map.keySet().join(', ')}"
        } else {
            log.info "No samplesheet provided – will process ALL detected barcodes."
        }

        // ========== 1. BUILD INPUT CHANNEL FROM DIRECTORIES ==========
        def barcodeDirs = file("${params.inDir}/barcode*", type: 'dir', maxdepth: 1)

        Channel
            .fromPath(barcodeDirs)
            .filter(~/.*barcode\d+$/)
            .map { dir ->
                def dirName = dir.getName()
                def barcodeMatch = (dirName =~ /barcode(\d+)/)

                if (!barcodeMatch) return null

                def bcNum = barcodeMatch[0][1].padLeft(2, '0')
                def bcKey = "barcode${bcNum}"
                def sample

                if (hasSamplesheet) {
                    sample = samplesheet_map[bcKey]
                    if (!sample) {
                        log.debug "Skipping ${bcKey} - not in samplesheet"
                        return null
                    }
                } else {
                    sample = bcKey
                }

                // Count FASTQ files in directory
                def count = 0
                dir.listFiles().each { f ->
                    if (f.isFile() && (f.name.endsWith('.fastq') || f.name.endsWith('.fastq.gz'))) {
                        count++
                    }
                }

                if (count == 0) {
                    log.warn "No FASTQ files found in ${dirName}"
                    return null
                }

                log.info "Found ${count} FASTQ files in ${dirName} (sample: ${sample})"
                return [sample, dir]
            }
            .filter { it != null }
            .ifEmpty { error "No barcode directories with FASTQ files found in ${params.inDir}" }
            .set { ch_barcode_dirs }

        // ========== 2. RUN GUPPYPLEX ==========
        guppyplex(ch_barcode_dirs)

        // ========== 3. RUN ARTIC ==========
        artic(
            guppyplex.out.fastq,
            params.bedFile,
            params.refGenome,
            params.model
        )

        // ========== 4. CONSENSUS SEQUENCES ==========
        artic.out.consensus
            .map { it[1] }
            .collect()
            .set { all_consensus_fastas }

        nextclade(all_consensus_fastas)

        // ========== 5. BAM STATS ==========
        artic.out.sortedBam
            .flatMap { sample, bam, bai -> [bam, bai] }
            .collect()
            .set { all_bam_files }

        bamstats(all_bam_files, nextclade.out.tsv)
}