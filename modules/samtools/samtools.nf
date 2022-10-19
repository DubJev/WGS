#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process samView {

  label "bwa"
  publishDir params.viewOutputs, mode: 'copy'
  
  input:
    path ch_bam
  output:
    path alignmentsBam
    path headerBam
    path "samtools_version.yml", emit: version
  script:
    """
    samtools view ${ch_bam} > alignmentsBam
    samtools view -H ${ch_bam} > headerBam

    export SAMTOOLS_VER=\$(samtools --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
    echo samtools: \$SAMTOOLS_VER > samtools_version.yml
    """
}

workflow.onError{
    println "Stopped: ${workflow.errorMessage}"
}

workflow viewSam {
  take:
    ch_bam
  main:
    samView(ch_bam)
}

workflow {
  ch_bam = Channel.fromPath(params.bams+"*.bam")
  viewSam(ch_bam)
}
