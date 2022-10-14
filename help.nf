log.info"""
============================================================
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Simple WGS Pipline - Modules: FastQC, BWA MEM and Haplotyper
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Command for running:

    nextflow run /path_to_project/main.nf -c /path_to_project/configs/main.config

------------------------------------------------------------
Input:
    
    Input folder is defined in main.config

    Default input folder /path_to_project/tests/main/Inputs
    Mandatory files: 
        reference FASTA file, 
        reads in FASTQ format (2 files for paired-end reads)
        VCF.GZ file with known polymorphic sites

    Mandatory docker images are defined in environment.config

------------------------------------------------------------
Output:

    Output folders are defined in main.config

    Default output folder /path_to_project/tests/main/Outputs
        Variants in VCF file
    
    Default folder for versions of used software /path_to_project/tests/main/Outputs/Versions
        Versions of Samtools, BWA, GATK and FastQC in yaml files

    Default folder for results of FastQC /path_to_project/tests/main/Outputs/FQCheck
        Results for all input FASTQ files

    Default folder for auxiliary results /path_to_project/tests/main/Outputs/Temp
        Results of alignment, deduplication and recalibration
=============================================================
"""//.stripIndent()