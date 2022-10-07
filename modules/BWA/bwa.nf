#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process alignFiles {

  tag "${sample_name}"
  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    path ch_fasta
    tuple val(sample_name), path(reads)
  output:
    tuple val(sample_name), path("${sample_name}_aligned.bam"), emit: alignedBam
    path "version.yml"
  script:
    """
    bwa index ${ch_fasta}
    bwa mem -M ${ch_fasta} ${reads} | samtools sort > ${sample_name}_aligned.bam

    BWA_VER="bwa: "\$(echo \$(bwa 2>&1) | sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*-[A-Za-z]+[0-9]*)
    echo \$BWA_VER > bwa_version.yml
    SAMTOOLS_VER=\$(samtools --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]*[0-9]*)
    echo samtools: \$SAMTOOLS_VER > samtools_version.yml
    cat samtools_version.yml bwa_version.yml > version.yml
    """
}

process extractChr20 {

  publishDir params.alignOutputs, mode: 'copy'

  input:
    tuple val(sample_name), path(largeBam)
  output:
    tuple val(sample_name), path("${sample_name}_aligned_chr20.bam"), emit: alignedChr20Bam
  script:
    """
    samtools index ${largeBam}
    samtools view -h ${largeBam} 20 > ${sample_name}_aligned_chr20.bam
    """
}

process deDuplicate {

  publishDir params.alignOutputs, mode: 'copy'
  
  input:
    tuple val(sample_name), path(ch_bam)
  output:
    tuple val(sample_name), path("${sample_name}_deduplicate.bam"), emit: deduplicate
    path "${sample_name}_dMetrics.txt"
  script:
    """
    gatk AddOrReplaceReadGroups I=${ch_bam} O=output.bam RGLB=lib1 RGPL=ILLUMINA RGPU=50 RGSM=TAAGGCGA
    gatk MarkDuplicates -I output.bam -M ${sample_name}_dMetrics.txt -O ${sample_name}_deduplicate.bam
    """
}

process reCalibrate {

  publishDir params.alignOutputs, mode: 'copy'

  input:
    tuple val(sample_name), path(de_bam)
    path ch_fasta
    path ch_vcfs
  output:
    path "${sample_name}_recalibrated.bam"
    
  script:
    """
    samtools faidx ${ch_fasta}
    gatk CreateSequenceDictionary R=${ch_fasta}
    gatk IndexFeatureFile -I ${ch_vcfs}
    gatk BaseRecalibrator -I ${de_bam} -R ${ch_fasta} --known-sites ${ch_vcfs} -O reTabela.table 
    gatk ApplyBQSR -R ${ch_fasta} -I ${de_bam} --bqsr-recal-file reTabela.table -O ${sample_name}_recalibrated.bam
    """
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

workflow alignFiles_wf {
  take:
    ch_fasta
    ch_fastqs
  main:
    alignFiles(ch_fasta, ch_fastqs)
  emit:
    alignedBam = alignFiles.out.alignedBam
}

workflow extractChr20_wf {
  take:
    ch_bam
  main:
    extractChr20(ch_bam)
  emit:
    alignedChr20Bam = extractChr20.out.alignedChr20Bam
}

workflow deDuplicate_wf {
  take:
    ch_bam
  main:
    deDuplicate(ch_bam)
  emit:
    deduplicate = deDuplicate.out.deduplicate

}

workflow reCalibrate_wf {
  take:
    de_bam
    ch_fasta
    ch_vcfs
  main:
    reCalibrate(de_bam, ch_fasta, ch_vcfs)

}

workflow {
  ch_fasta = Channel.fromPath(params.alignInputs+"*.fasta")
  ch_fastqs = Channel.fromFilePairs(params.alignInputs+"*_{1,2}.fastq")
  ch_vcfs = Channel.fromPath(params.alignInputs+"*.vcf.gz")
  alignFiles_wf(ch_fasta, ch_fastqs)
  extractChr20_wf(alignFiles_wf.out)
  deDuplicate_wf(extractChr20_wf.out)
  reCalibrate_wf(deDuplicate_wf.out, ch_fasta, ch_vcfs)
  
}
