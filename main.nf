include {checkVersion as FastQCVersion} from './modules/FastQC/fastqc.nf'
include {checkQ} from './modules/FastQC/fastqc.nf'
include {VersionBWASam} from './modules/BWA/bwa.nf'
include {alignFiles; extractChr20; } from './modules/BWA/bwa.nf'
include {deDuplicate} from './modules/BWA/bwa.nf'
include {reCalibrate} from './modules/BWA/bwa.nf'
include {Samtools_sort} from './modules/qualimap_bamqc/qualimap_bamqc.nf'
include {Qualimap_bamqc} from './modules/qualimap_bamqc/qualimap_bamqc.nf'
include {checkVersion as GATKVersion} from './modules/haplotyper/haplotyper.nf'
include {HCGVCF; VARCALL} from './modules/haplotyper/haplotyper.nf'
include {snpEff} from './modules/snpEff/snpEff.nf'
include {checkVersion as VEPVersion} from './modules/VEP/VEP.nf'
// include {VariantPredictor} from './modules/VEP/VEP.nf'
include {NFVersion} from './modules/ToolVersioning/ToolVersioning.nf'
include {ToolVersions} from './modules/ToolVersioning/ToolVersioning.nf'

workflow {

  main:

    ch_fastqs = Channel.fromFilePairs(params.mainIn+"*_{1,2}.fastq")
    ch_check = Channel.fromPath(params.mainIn+"*.fastq")
    ch_fasta = Channel.fromPath(params.mainIn+"*.fasta")
    ch_vcfs = Channel.fromPath(params.mainIn+"*.vcf.gz")
    //ch_recalbam = Channel.fromPath(params.recalbam)

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
    ch_recalbam = reCalibrate.out.recalibrate
    ch_recalbam.view()
    Samtools_sort(ch_recalbam)
    ch_bamsorted = Samtools_sort.out.bamsorted.view()
    Qualimap_bamqc(ch_bamsorted)
    HCGVCF(ch_fasta, reCalibrate.out.recalibrate)
    VARCALL(ch_fasta, HCGVCF.out.gvcfs)
    // VariantPredictor(VARCALL.out.variants)
    ch_haplotyper_vcf = VARCALL.out.variants
    ch_dataDir = params.dataDir
    ch_config = params.config
    snpEff(ch_haplotyper_vcf, ch_dataDir, ch_config)



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
