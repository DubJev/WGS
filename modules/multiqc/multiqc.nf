nextflow.enable.dsl=2

process MultiQC {
  label "multiqc"
  tag "${sample_name}"
  publishDir "${params.multiqc_res}", mode: 'copy'

  input:
  val qc_val
  path qc_path
  path wait_qc
  // tuple val(sample_name), path("*${sample_name}*.html")
    // tuple val(sample_name), path(fastqc_html)
  // path fastqc_html

  output:
  path "*"

  script:
  """
  multiqc -f ${qc_val}
  """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// workflow MultiQC {
//   take:
//     fastqc_html
//   main:
//     MultiQC(fastqc_html)
// }

workflow {
  ch_fastqc_html = Channel.fromPath(params.fastqc_html)
  MultiQC(ch_fastqc_html)
}
