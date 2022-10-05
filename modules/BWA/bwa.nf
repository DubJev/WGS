#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process alignFiles {

  tag "${sample_name}"
  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_fasta
    tuple val(sample_name), path(reads)
  output:
    path "${sample_name}_aligned.bam", emit: alignedBam
    val sample_name, emit: sample_name
  script:
    """
    bwa index ${ch_fasta}
    bwa mem -M ${ch_fasta} ${reads} | samtools sort > ${sample_name}_aligned.bam
    """
}

process deDuplicate {

  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_bam
    val (sample_name)
  output:
    path "${sample_name}_deduplicate.bam", emit: deduplicate
    path "${sample_name}_dMetrics.txt"
    val sample_name, emit: sample_name
  script:
    """
    gatk MarkDuplicates -I ${ch_bam} -M ${sample_name}_dMetrics.txt -O ${sample_name}_deduplicate.bam
    """
}

process reCalibrate {

  publishDir params.alignOutputs, mode: 'copy'

  input:
    path de_bam
    path ch_fasta
    path ch_vcfs
    val sample_name
  output:
    path "${sample_name}_recalibrated.bam"
    
  script:
    """
    samtools faidx ${ch_fasta}
    gatk CreateSequenceDictionary R=${ch_fasta}
    gatk IndexFeatureFile -I ${ch_vcfs}
    gatk AddOrReplaceReadGroups I=${de_bam} O=output.bam RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
    gatk BaseRecalibrator -I output.bam -R ${ch_fasta} --known-sites ${ch_vcfs} -O reTabela.table 
    gatk ApplyBQSR -R ${ch_fasta} -I ${de_bam} --bqsr-recal-file reTabela.table -O ${sample_name}_recalibrated.bam

    """
}

workflow alignFiles_wf {
  take:
    ch_fasta
    ch_fastqs
  main:
    alignFiles(ch_fasta, ch_fastqs)
  emit:
    alignedBam = alignFiles.out.alignedBam
    sample_name = alignFiles.out.sample_name
}

workflow deDuplicate_wf {
  take:
    ch_bam
    sample_name
  main:
    deDuplicate(ch_bam, sample_name)
  emit:
    deduplicate = deDuplicate.out.deduplicate

}

workflow reCalibrate_wf {
  take:
    de_bam
    ch_fasta
    ch_vcfs
    sample_name
  main:
    reCalibrate(de_bam, ch_fasta, ch_vcfs, sample_name)

}

workflow {
  ch_fasta = Channel.fromPath(params.alignInputs+"*.fasta")
  ch_fastqs = Channel.fromFilePairs(params.alignInputs+"*_{1,2}.fastq")
  ch_vcfs = Channel.fromPath(params.alignInputs+"*.vcf.gz")
  alignFiles_wf(ch_fasta, ch_fastqs)
  deDuplicate_wf(alignFiles_wf.out.alignedBam, alignFiles_wf.out.sample_name)
  reCalibrate_wf(deDuplicate_wf.out, ch_fasta, ch_vcfs, alignFiles_wf.out.sample_name)
}
