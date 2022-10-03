#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process alignFiles {

  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_fasta
    path ch_fastq1
    path ch_fastq2
  output:
    path "aligned.bam", emit: alignedBam
  script:
    """
    bwa index ${ch_fasta}
    bwa mem -M ${ch_fasta} ${ch_fastq1} ${ch_fastq2} | samtools sort > aligned.bam
    """
}

process deDuplicate {

  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_bam
  output:
    path "duplicate.bam"
    path "dMetrics.txt"
  script:
    """
    gatk MarkDuplicates -I ${ch_bam} -M dMetrics.txt -O duplicate.bam
    """
}

workflow alignFiles_wf {
  take:
    ch_fasta
    ch_fastq1
    ch_fastq2
  main:
    alignFiles(ch_fasta, ch_fastq1, ch_fastq2)
  emit:
    alignedBam = alignFiles.out.alignedBam
}

workflow deDuplicate_wf {
  take:
    ch_bam
  main:
    deDuplicate(ch_bam)
}

workflow {
  ch_fasta = Channel.fromPath(params.fastas)
  ch_fastq1 = Channel.fromPath(params.fastq1)
  ch_fastq2 = Channel.fromPath(params.fastq2)
  alignFiles_wf(ch_fasta, ch_fastq1, ch_fastq2)
  deDuplicate_wf(alignFiles_wf.out)
}
