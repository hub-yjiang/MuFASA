[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17531746.svg)](https://doi.org/10.5281/zenodo.17531746)

# MuFASA
Mutation Filtering and Analysis of Somatic Alterations.  
This pipeline will map small sequencing reads to a publicly available reference genome with bwa-mem. After mapping, small somatic variants are called with Strelka and Mutect2. SNV's and INDELs are filtered on quality by several filters.

# Installation Guide

## Nextflow & Conda
Make sure you have [nextflow](https://www.nextflow.io/docs/latest/index.html) available.   
This pipeline is optimized to be run on a cluster with scheduling software like SLURM available.

The Nextflow pipeline primarily uses conda environments. To install the proper tools and dependencies, use the conda env requirements (.yml) files. All Tools and Packages used in the nextflow workflow are specified in these files. 

## workflow configuration

### nextflow specific
Modify the resources.config configuration file with your run information.  
SLURM is set as default in the nextflow.config file. If you don't have a SLURM setup available, modify this to the appropriate executor for your setup.

### general configuration / Input
samples.csv, comma seperated file with experiment information (sampleid, fastq.R1, fastq.R2)   
If you have multiple lanes for one sample, add a second row for the same sampleid.  
tumor_normal_pairs.txt, comma seperated file with information for the somatic variant calling.
 On row 1: (tumor/treated) sampleid, tumor, nr  
 on accompanying second row, same nr: (normal/untreated) sampleid, normal, nr  
(nr increments per sample pairs (two rows))  
filter_parameters.txt. This is the configuration used for [FiNGS](https://pubmed.ncbi.nlm.nih.gov/33602113/)   
This tool is used for filtering SNP's after variant calling

Inputs can also be mapped bam files, the samples.csv need to be a comma seperated file with a header of "sampleID,bam,bai", following the absolute path of the bam and bai files. 

## workflow testing
In the TESTDATA dir is a testset with fastq data to test your installation.  
If your installation is correct, you will be able to run this data set by running:

`./run_nextflow.sh`

### workflow output

| File                                | location        | output                                                        | 
|-------------------------------------|-----------------|---------------------------------------------------------------|
| sampleID.fastqc.html                | FASTQC/         | Quality control checks on raw sequence data                   |
| sampleID_collect_align_metrics.txt  | FASTQC/         | Quality metrics for illumina read alignments                  |
| sampleID_collect_wgs_metrics.txt    | FASTQC/         | Metrics about coverage and performance of (WGS) mapping       |
| sampleID_samtools_markdup.(bam/bai) | MAPPING/        | Read mapping output                                           |
| sampleID.somatic.snvs & indel.vcf   | Strelka/        | Strelka somatic SNV's  and InDel variants                     |
| sampleID_Mutect.snvs & indel.vcf    | Mutect/         | Mutect2 somatic SNV's and InDel variants                      |
| sampleID_PASS.vcf                   | Filtered/       | Strelka and Mutect2 PASS filtered SNV's and InDel variants    |
| sampleID_shared.snvs & indel.vcf    | Filtered/       | Intersected results for Mutect/Strelka                        |
| sampleID_pon.snvs & indel.vcf       | Filtered/       | Customized filtering for Panel of Normal samples              |
| sampleID.filtered.vcf               | snvs_Filtered/  | FiNGS filtered SNV's                                          |
| sampleID.plots.pdf                  | snvs_Filtered/  | FiNGS report                                                  |
| sampleID.vaf.vcf                    | snvs_Filtered/  | SNV's filtered on Variant Allelic Frequency (VAF)             |
| sampleID.AF.vcf                     | indel_Filtered/ | InDels final filtered output file. (filters defined in script)|


## FASTQ processing
All paired-end reads are assesed on read quality and multiple sequencing runs are merged before trimming.   
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) is used to find and remove adapter sequences, primers and poly-A tails. 

## FASTQC 
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Picard CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard) and [Picard CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard) produce reports that can be combined by running:

`multiqc ./FASTQC/`


## MAPPING processing
All reads are mapped against the whole genome reference defined in the resources.config.

## Variant Calling
[Strelka](https://github.com/Illumina/strelka) and [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037269791-Mutect2-BETA) are used for variant calling.  
Somatic SNV's and InDels that are shared between both callers are PASSED.

## Variant Filtering
SNV's and InDels are filtered on quality and additional custimized filters.

## Structural Variant Calling

Structural Variants are called with Manta and Gridss. 
Results are filtered on PASS and encode blacklists and variants called by both Manta and Gridss are presented as output in .bedpe and .vcf format


## Citing MuFASA
For more information and if you use MuFASA in your work, please cite:

 __Jiang, Y., Przybilla, M., Bakker, L. et al. Tissue-specific mutagenesis from endogenous guanine damage is suppressed by Polκ and DNA repair. Nat Commun 17, 436 (2026). https://doi.org/10.1038/s41467-025-67072-1__

## References

_Strelka2: fast and accurate calling of germline and somatic variants. Kim, S., Scheffler, K. et al. (2018), Nature Methods, 15, 591-594._

_FiNGS: high quality somatic mutations using filters for next generation sequencing. Wardell CP, Ashby C, Bauer MA, BMC Bioinformatics. 2021;22(1):77. Published 2021 Feb 18. doi:10.1186/s12859-021-03995-y_
