#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process NFVersion {

  publishDir params.verInputs, mode: 'copy'
  
  output:
    path "Nextflow_version.yml", emit: nfversion
  script:
    """
    export NXTFLW_VER=\$(~/nextflow -version |sed -n -e '3p' | grep -Eo [0-9]+[.][0-9]*[.][0-9]*" build "[0-9]*)
    echo Nextflow: \$NXTFLW_VER > Nextflow_version.yml
    """
}

process ToolVersions {

  publishDir params.reportOutputs, mode: 'copy'

  input:
    path tools
  output:
    path "WGS_tool_version.yml"
  script:
    """
    cat WGS_tool_version.yml $tools > WGS_tool_version.yml
    """
}

workflow.onError{
    println "Stopped: ${workflow.errorMessage}"
}

workflow NFVersion_wf {

  main:
    NFVersion()
  emit:
    nfversion = NFVersion.out.nfversion
}


workflow {
  NFVersion_wf()
  tools = Channel.fromPath(params.verInputs+"*.yml").toList()
  ToolVersions(tools)
}
