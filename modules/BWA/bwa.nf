#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process alignFiles {

  tag "${sample_name}"
  publishDir params.alignOutputs, mode: 'copy', saveAs: {filename -> "${sample_name}_aligned.bam"}
  
  input:
    path ch_fasta
    tuple val(sample_name), path(reads)
  output:
    path "aligned.bam", emit: alignedBam
    val sample_name, emit: sample_name
  script:
    """
    bwa index ${ch_fasta}
    bwa mem -M ${ch_fasta} ${reads} | samtools sort > aligned.bam
    """
}

process deDuplicate {

  publishDir params.alignOutputs, mode: 'copy', saveAs: { filename -> filename.endsWith(".bam") ? "${sample_name}_deduplicate.bam" : "${sample_name}_dMetrics.txt"}
  
  input:
    path ch_bam
    val (sample_name)
  output:
    path "deduplicate.bam", emit: deduplicate
    path "dMetrics.txt"
  script:
    """
    gatk MarkDuplicates -I ${ch_bam} -M dMetrics.txt -O deduplicate.bam
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
}

workflow {
  ch_fasta = Channel.fromPath(params.alignInputs+"*.fasta")
  ch_fastqs = Channel.fromFilePairs(params.alignInputs+"*_{1,2}.fastq", siza:-1)
  alignFiles_wf(ch_fasta, ch_fastqs)
  deDuplicate_wf(alignFiles_wf.out.alignedBam, alignFiles_wf.out.sample_name)
}
