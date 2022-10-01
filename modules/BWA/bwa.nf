#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process indexFiles {
  
  publishDir params.indexOutputs, mode: 'copy'
  
  input:
    path ch_fasta
  output:
    path "example.fasta.*"
  script:
    """
    bwa index ${ch_fasta}
    """
}

process alignFiles {

  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_fasta
    path ch_fastq1
    path ch_fastq2
  output:
    path "aligned.sam"
  script:
    """
    bwa mem ${ch_fasta} ${ch_fastq1} ${ch_fastq2} > aligned.sam
    """
}

workflow makeFai {
  take:
    ch_fasta
  main:
    indexFiles(ch_fasta)
}

workflow makeSam {
  take:
    ch_fasta
    ch_fastq1
    ch_fastq2
  main:
    alignFiles(ch_fasta, ch_fastq1, ch_fastq2)
}

workflow {
  ch_fasta = Channel.fromPath(params.fastas)
  ch_fastq1 = Channel.fromPath(params.fastq1)
  ch_fastq2 = Channel.fromPath(params.fastq2)
  makeFai(ch_fasta)
//  makeSam(ch_fasta, ch_fastq1, ch_fastq2)
}
