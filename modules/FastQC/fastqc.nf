#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process checkQ {

  publishDir params.qcOutputs, mode: 'copy'
  
  input:
    path(reads)
  output:
    path "*"
  script:
    """
    fastqc ${reads}
    """
}

process checkVersion {

  publishDir params.qcOutputs, mode: 'copy'
  output:
    path "fastqc_version.yml"
  script:
  """
    FASTQC_VER=\$(fastqc --version 2>&1 |  sed -n -e '1p' | grep -Eo [0-9][.]+[0-9]*[.]+[0-9]*)
    echo FastQC: \$FASTQC_VER > fastqc_version.yml
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

workflow checkQ_wf {
  take:
    ch_reads
  main:
    checkQ(ch_reads)

}

workflow checkVersion_wf {
  main:
    checkVersion()
}

workflow {
  ch_reads = Channel.fromPath(params.qcInputs+"*")
  checkQ_wf(ch_reads)
  checkVersion_wf()
}
