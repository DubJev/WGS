#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process snpEff {

  label "snpEff"
  tag "${sample_name}"
  publishDir params.annvcf, mode: 'copy'

  input:
  tuple val(sample_name), path(variants)
  path dataDir
  path config

  output:
  tuple val(sample_name), path ("${sample_name}_ann20.vcf"), emit: annvcf

  script:
  """
  java -Xmx8g -jar /home/biodocker/bin/snpEff/snpEff.jar -dataDir ${dataDir} -c ${config} humanchr20 ${variants} > ${sample_name}_ann20.vcf
  """
  }

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow snpEff_wf {
  take:
  variants
  dataDir
  config
  main:
  snpEff(variants, dataDir, config)
  emit:
  annvcf = snpEff.out.annvcf

}

workflow {
  // ch_vcf = Channel.fromPath(params.vcf)
  variants = Channel.fromPath(params.variants)
  ch_dataDir = Channel.fromPath(params.dataDir)
  ch_config = Channel.fromPath(params.config)
  snpEff(variants, ch_dataDir, ch_config)
}
