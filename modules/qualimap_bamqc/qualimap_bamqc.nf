nextflow.enable.dsl=2

process Samtools_sort {
  publishDir params.qualimap_bamqc_res , mode: 'copy'

  input:
  path(recalbam)

  output:
  path "bamsorted", emit: bamsorted

  script:
  """
  samtools sort ${recalbam} > bamsorted
  """
}

process Qualimap_bamqc {
  publishDir params.qualimap_bamqc_res, mode 'copy'
  input:
  path(bamsorted)

  script:
  """
  qualimap bamqc -bam ${bamsorted} -outdir $params.qualimap_bamqc_res
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
    recalbam
  main:
    Samtools_sort(recalbam)
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
    ch_recalbam = Channel.fromPath(params.recalbam)
    Samtools_sort(ch_recalbam)
    Qualimap_bamqc(Samtools_sort.out)
  }
