#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process HCGVCF {
  
  label "gatk"
  tag "${sample_name}"
  publishDir params.hcOutputs, mode: 'copy'
  
  input:
    path ch_fasta
    tuple val(sample_name), path(re_bam)
  output:
    tuple val(sample_name), path("${sample_name}.g.vcf"), emit: gvcfs
  script:
    """
    samtools faidx ${ch_fasta}
    samtools index ${re_bam}
    gatk CreateSequenceDictionary R=${ch_fasta}
    gatk HaplotypeCaller -R ${ch_fasta} -I ${re_bam} -ERC GVCF -O ${sample_name}.g.vcf
    """
}

process VARCALL {

  label "gatk"
  publishDir params.hcOutputs, mode: 'copy'

  input:
    path ch_fasta
    tuple val(sample_name), path(gvcf)
  output:
    tuple val(sample_name), path("${sample_name}.vcf"), emit: variants
  script:
    """
    samtools faidx ${ch_fasta}
    gatk CreateSequenceDictionary R=${ch_fasta}
    gatk GenotypeGVCFs -R ${ch_fasta} -V ${gvcf} -O ${sample_name}.vcf
    """
}

process checkVersion {

  label "gatk"
  publishDir params.verOutputs, mode: 'copy'
  output:
    path "gatk_version.yml", emit: gatkversion
  script:
  """
    GATK_VER=\$(gatk -version |  sed -n -e '1p' | grep -Eo [0-9][.]+[0-9]*[.]+[0-9]*[.]+[0-9]*)
    echo GATK: \$GATK_VER > gatk_version.yml
    """
}

workflow.onError{
    println "Stopped: ${workflow.errorMessage}"
}

workflow HCGVCF_WF {
  take:
    ch_fasta
    ch_bams
  main:
    HCGVCF(ch_fasta, ch_bams)
  emit:
    gvcfs = HCGVCF.out.gvcfs
}

workflow VARCALL_WF {
  take:
    ch_fasta
    ch_bam
  main:
    VARCALL(ch_fasta, ch_bam)
  emit:
    variants = VARCALL.out.variants
}

workflow checkVersion_wf {
  main:
    checkVersion()
  emit:
    gatkversion = checkVersion.out.gatkversion
}

workflow {
  ch_fasta = Channel.fromPath(params.hcInputs+"*.fasta")
  ch_bam = Channel.fromFilePairs(params.hcInputs+"*.bam", size:-1)
  HCGVCF_WF(ch_fasta, ch_bam)
  VARCALL_WF(ch_fasta, HCGVCF_WF.out)
  checkVersion_wf()
  
}
