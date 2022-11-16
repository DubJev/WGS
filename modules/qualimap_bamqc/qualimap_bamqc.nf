nextflow.enable.dsl=2

process Samtools_sort {
  label "samtools"
  tag "${sample_name}"
  publishDir params.qualimap_bamqc_res , mode: 'copy'

  input:
  tuple val(sample_name), path(recalibrate)

  output:
  tuple val(sample_name), path("${sample_name}_sorted.bam") , emit: bamsorted

  script:
  """
  samtools sort ${recalibrate} > ${sample_name}_sorted.bam
  """
}

process Qualimap_bamqc {
  label "qualimap"
  tag "${sample_name}"
  publishDir params.qualimap_bamqc_res 
  input:
  tuple val(sample_name), path(bamsorted)
  output:
  path("qualimap_outdir_${sample_name}"), emit: outdir
  script:
  """
  qualimap bamqc -bam ${bamsorted} -outdir qualimap_outdir_${sample_name}
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
workflow Samtools_sort_wf {
  take:
    recalibrate
  main:
    Samtools_sort(recalibrate)
  emit:
    bamsorted = Samtools_sort.out.bamsorted
}
workflow Qualimap_bamqc_wf {
  take:
    bamsorted
  main:
    Qualimap_bamqc(bamsorted)
  emit:
    bamqc = Qualimap_bamqc.out.bamqc
  }

  workflow {
    // ch_recalbam = Channel.fromPath(params.recalbam)
    Samtools_sort(ch_recalbam)
    Qualimap_bamqc(Samtools_sort.out)
  }
