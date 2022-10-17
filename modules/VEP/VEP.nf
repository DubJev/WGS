#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process VariantPredictor {
  
  label "vep"
  tag "${sample_name}"
  publishDir params.vpOutputs, mode: 'copy'
  
  input:
    tuple val(sample_name), path(vcfs)
  output:
    tuple val(sample_name), path("${sample_name}_variant_effect_output.txt"), emit: vareff
  script:
    """
    vep --input_file ${vcfs} --output_file ${sample_name}_variant_effect_output.txt  --species homo_sapiens --database --force_overwrite

    """
        //vep --offline --cache --dir_cache /opt/vep/.vep --input_file /opt/vep/.vep/input/${vcfs} --output_file /opt/vep/.vep/output/${sample_name}_variant_effect_output.txt  

}

process checkVersion {

  label "vep"
  publishDir params.vpOutputs, mode: 'copy'
  output:
    path "vep_version.yml"
  script:
  """
    VEP_VER=\$(vep | sed -n -e '10p' | grep -Eo [0-9]+[.]+[A-Za-z0-9]*)
    echo VEP: \$VEP_VER > vep_version.yml
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

workflow VariantPredictor_WF {
  take:
    ch_vcf
  main:
    VariantPredictor(ch_vcf)
  emit:
    vareff = VariantPredictor.out.vareff
}

workflow checkVersion_wf {
  main:
    checkVersion()
}

workflow {
  ch_vcf = Channel.fromFilePairs(params.vpInputs+"*.vcf", size:-1)
  VariantPredictor_WF(ch_vcf)
  checkVersion_wf()
  
}
