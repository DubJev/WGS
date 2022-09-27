#!/usr/bin/env nextflow

params.viewOutput = "tests/samtools/outputs"
params.refBam = "example.bam"

Channel.value("$workflow.launchDir/tests/samtools/inputs/$params.refBam").set { ch_bam }

process samView {

  publishDir params.viewOutput, mode: 'copy'
  container 'djevtic/samtools'
  
  input:
    path refBam from ch_bam
  output:
    path "alignmentsBam" into ch_out_view
    path "headerBam" into ch_out_headerview

  """
  samtools view $params.refBam > alignmentsBam
  samtools view -H $params.refBam > headerBam
  """
}

