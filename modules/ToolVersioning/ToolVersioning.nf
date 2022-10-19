#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process ToolVersions {

  publishDir params.reportOutputs, mode: 'copy'
  
  input:
    path tools
  output:
    path "Tool_version.yml", emit: version
  script:
    """
    export NXTFLW_VER=\$(~/nextflow -version |sed -n -e '3p' | grep -Eo [0-9]+[.][0-9]*[.][0-9]*" build "[0-9]*)
    echo Nextflow: \$NXTFLW_VER > Nextflow_version.yml
    
    cat Nextflow_version.yml $tools > Tool_version.yml
    """
}

workflow.onError{
    println "Stopped: ${workflow.errorMessage}"
}

workflow ToolVersions_wf {
  take:
    tools
  main:
    ToolVersions(tools)
  emit:
    version = ToolVersions.out.version
}

workflow {
  tools = Channel.fromPath(params.verInputs+"*.yml").collect()
  ToolVersions_wf(tools)
}
