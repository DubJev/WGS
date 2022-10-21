+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Simple WGS Pipline - Modules: FastQC, BWA MEM, Haplotyper and VEP
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Bioinformatics analysis pipeline for WGS (Whole Genome Sequencing) data.
The pipeline is built using Nextflow, a bioinformatics workflow tool to run tasks across multiple compute infrastructures.

It pre-processes raw data from FastQ inputs, aligns the reads, marks duplicates and recalibrates.

FastQC is a program designed to spot potential problems in high-throughput sequencing datasets. 
It runs a set of analyses on one or more raw sequence files in fastq format and produces a report which summarizes the results.

BWA is a software package for mapping DNA sequences against a large reference genome, such as the human genome. 
BWA index is used to index input database sequences in the FASTA format. 
Sequence alignment is done with the BWA-MEM algorithm. The algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW). The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard’s markDuplicates does not work with split alignments. Option -M is used to flag shorter split hits as secondary (for Picard compatibility).

GATK (Genome Analysis Toolkit) is used for variant discovery in high-throughput sequencing data. MarkDuplicate tool is used for second processing step and it is performed per-sample and consists of identifying read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artifactual processes. These are considered to be non-independent observations, so the program tags all but a single read pair within each set of duplicates, causing the marked pairs to be ignored by default during the variant discovery process. At this stage the reads also need to be sorted into coordinate-order for the next step of the processing. 
This third processing step is performed per-sample and consists of applying machine learning to detect and correct for patterns of systematic errors in the base quality scores, which are confidence scores emitted by the sequencer for each base. Base quality scores play an important role in weighing the evidence for or against possible variant alleles during the variant discovery process, so it's important to correct any systematic bias observed in the data. Biases can originate from biochemical processes during library preparation and sequencing, from manufacturing defects in the chips, or instrumentation defects in the sequencer. The recalibration procedure involves collecting covariate measurements from all base calls in the dataset, building a model from those statistics, and applying base quality adjustments to the dataset based on the resulting model. 

The next step we begin by calling variants per sample in order to produce a file in GVCF format. The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. In the GVCF mode used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate file called a GVCF , which can then be used for joint genotyping of multiple samples in a very efficient way. We gather all the per-sample GVCFs and pass them all together to the joint genotyping tool, GenotypeGVCFs. This produces a set of joint-called SNP and indel calls.

VEP (Variant Effect Predictor) predicts the functional effects of genomic variants. VEP needs some input containing variant positions to run. We use cache file to run VEP on the supplied file homo_sapiens_GRCh37.vcf.

![](image.png)

Command for running:

    nextflow run /path_to_project/main.nf -c /path_to_project/configs/main.config

-----------------------------------------------------------------
Input:
    
    Input folder is defined in main.config

    Default input folder /path_to_project/tests/main/Inputs
    Mandatory files: 
        reference FASTA file, 
        reads in FASTQ format (2 files for paired-end reads)
        VCF.GZ file with known polymorphic sites

    Docker images that are used are defined in environment.config

-----------------------------------------------------------------
Output:

    Output folders are defined in main.config

    Default output folder /path_to_project/tests/main/Outputs
        Variants in VCF file
        WGS_tool_version.yml
    
    Default folder for versions of used tools /path_to_project/tests/main/Outputs/Versions
        Versions of Nextflow, Samtools, BWA, GATK, VEP and FastQC in yaml files

    Default folder for results of FastQC /path_to_project/tests/main/Outputs/FQCheck
        Results for all input FASTQ files

    Default folder for auxiliary results /path_to_project/tests/main/Outputs/Temp
        Results of alignment, deduplication, recalibration and variants
==================================================================
