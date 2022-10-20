include {checkVersion as FastQCVersion} from './modules/FastQC/fastqc.nf'
include {checkQ} from './modules/FastQC/fastqc.nf'
include {VersionBWASam} from './modules/BWA/bwa.nf'
include {alignFiles; extractChr20; } from './modules/BWA/bwa.nf'
include {deDuplicate} from './modules/BWA/bwa.nf'
include {reCalibrate} from './modules/BWA/bwa.nf'
include {checkVersion as GATKVersion} from './modules/haplotyper/haplotyper.nf'
include {HCGVCF; VARCALL} from './modules/haplotyper/haplotyper.nf'
include {checkVersion as VEPVersion} from './modules/VEP/VEP.nf'
include {VariantPredictor} from './modules/VEP/VEP.nf'
include {NFVersion} from './modules/ToolVersioning/ToolVersioning.nf'
include {ToolVersions} from './modules/ToolVersioning/ToolVersioning.nf'

workflow {

  main:

    ch_fastqs = Channel.fromFilePairs(params.mainIn+"*_{1,2}.fastq")
    ch_check = Channel.fromPath(params.mainIn+"*.fastq")
    ch_fasta = Channel.fromPath(params.mainIn+"*.fasta")
    ch_vcfs = Channel.fromPath(params.mainIn+"*.vcf.gz")
    FastQCVersion()
    VersionBWASam()
    GATKVersion()
    VEPVersion()
    NFVersion()

    checkQ(ch_check)
    alignFiles(ch_fasta, ch_fastqs)
    extractChr20(alignFiles.out.alignedBam)
    deDuplicate(extractChr20.out.alignedChr20Bam)
    reCalibrate(deDuplicate.out.deduplicate, ch_fasta, ch_vcfs)
    HCGVCF(ch_fasta, reCalibrate.out.recalibrate)
    VARCALL(ch_fasta, HCGVCF.out.gvcfs)
    VariantPredictor(VARCALL.out.variants)

    tools = FastQCVersion.out.mix(VersionBWASam.out, GATKVersion.out).mix(VEPVersion.out, NFVersion.out).collect()
    ToolVersions(tools)
}

workflow.onComplete{
    println "Status: ${ workflow.success ? 'OK' : 'failed' }"
    println """Completed at: $workflow.complete
               Duration: $workflow.duration
               WorkDir:  $workflow.workDir
             """
}

workflow.onError{
    println "Stopped: ${workflow.errorMessage}"
}