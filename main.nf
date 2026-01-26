#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ================================================================
// ========================== HELP MESSAGE =========================
// ================================================================
def helpMessage = """
🧬🧬 💻 **classiPox Pipeline** 💻 🧬🧬

Usage:
    nextflow run classiPox [options]

🔧 **Options**
    --help                  Show this help message and exit
    --inDir                 Input directory  
                               • Illumina: reads/  
                               • ONT:      fastq_pass/
    --outDir                Output directory (default: results/)
    --refGenome             Reference genome FASTA (default: clade_iib)
    --krakenDb              Kraken database (Illumina shotgun only)
    --bedFile               ARTIC primer BED (default: clade_iib)
    --samplesheet           CSV (for nanopore ARTIC only):
                            e.g.    sample_id,barcode  
                                    MPOX25-01001,1
                                    MPOX25-01002,2
    --illuminaShotgun       Run the Illumina shotgun workflow
    --nanoporeArtic         Run the Nanopore ARTIC workflow

⚠️ **Important**
You must choose **exactly one** workflow:  
    --illuminaShotgun **or** --nanoporeArtic
""".stripIndent()

if (params.help) {
    println helpMessage
    exit 0
}

// ================================================================
// ====================== PARAMETER DEFAULTS ======================
// ================================================================
params.outDir = params.outDir ?: 'results'

// ================================================================
// ====================== VALIDATE WORKFLOW ======================
// ================================================================
def workflows = [params.illuminaShotgun, params.nanoporeArtic].findAll { it }
if (workflows.size() != 1) {
    println "ERROR: Choose exactly one workflow: --illuminaShotgun or --nanoporeArtic"
    println helpMessage
    exit 1
}

// ================================================================
// =================== VALIDATE REQUIRED PARAMS ===================
// ================================================================
if (!params.inDir) {
    error "Missing required parameter: --inDir"
}

// ================================================================
// ===================== PRETTY STARTUP BANNER ====================
// ================================================================
def banner = """
╔═══════════════════════════════════════════════════════════════════╗
             🦠 classiPox — MPXV Classification Pipeline           
╚═══════════════════════════════════════════════════════════════════╝

🚀 **Workflow:**   ${params.nanoporeArtic ? 'Nanopore ARTIC' : 'Illumina Shotgun'}
📂 **Input:**      ${params.inDir}
📦 **Output:**     ${params.outDir}
🧫 **Reference:**  ${params.refGenome}
🧩 **BED File:**   ${params.bedFile}

""".stripIndent()
println banner



// ================================================================
// ==================== SAMPLESHEET HANDLING =====================
// ================================================================
def samplesheetMap = null

if (params.nanoporeArtic && params.samplesheet) {
    def sheet = file(params.samplesheet)
    if (!sheet.exists()) {
        error "Samplesheet not found: ${params.samplesheet}"
    }
    if (!params.samplesheet.endsWith('.csv')) {
        error "Samplesheet must be a CSV file (got: ${params.samplesheet})"
    }

    def rawText = sheet.getText('UTF-8')
    if (rawText.startsWith('\uFEFF')) {
        rawText = rawText.substring(1)
    }
    rawText = rawText.replaceAll('\r\n', '\n').trim()

    def cleanSheet = file("${workflow.workDir}/samplesheet_clean.csv")
    cleanSheet.text = rawText

    def map = [:]
    cleanSheet.splitCsv(header: true, sep: ',').each { row ->
        def sid = row.sample_id?.toString()?.trim()
        def bc  = row.barcode?.toString()?.trim()
        if (!sid || !bc) {
            error "Invalid row in samplesheet: ${row}"
        }
        def barcodeKey = bc.padLeft(2, '0')
        map["barcode${barcodeKey}"] = sid
    }

    samplesheetMap = map
    log.info "Loaded ${samplesheetMap.size()} samples from ${params.samplesheet}"
    log.debug "Samplesheet mapping: ${samplesheetMap}"
}

// ================================================================
// ========================= IMPORT WORKFLOWS =====================
// ================================================================
include { illuminaShotgun } from './workflows/illumina-wf.nf'
include { nanoporeArtic   } from './workflows/ont-wf.nf'

// ================================================================
// ============================ WORKFLOW ==========================
// ================================================================
workflow {
    if (params.illuminaShotgun) {
        illuminaShotgun()
    } else if (params.nanoporeArtic) {
        nanoporeArtic(samplesheetMap)
    }
}
