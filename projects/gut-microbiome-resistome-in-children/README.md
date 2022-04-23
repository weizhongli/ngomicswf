# Antibiotics alter human gut microbiome resistome and colonization resistance against antimicrobial resistant species in children

This folder include the workflow configration and parameters for this study, 
with which the analysis of this study can be reproduced.

The analysis workflow was done using NG-Omics-WF. To understand the
steps in this file, please first read the documents of NG-Omics-WF
at https://github.com/weizhongli/ngomicswf/blob/master/README.md 
and at https://github.com/weizhongli/ngomicswf/wiki.

## System Requirements for this study
This study requires a generic Linux computer or Linux cluster environment. 
Any common Linux systems, such as Ubuntu, CentOS, Fedora, Debian, are good. 
It is recommanded to run the analysis on computers with >=64GB RAM and >=32 cores. 

## Software Requirements for this study
Python 3.6 or higher is required to run NG-Omics-WF. But many other software tools are needed 
for the the analysis in the individual steps in the workflow. 

### Software Requirements
Required for most analysis
* Trimmomatic - A flexible read trimming tool for Illumina NGS data, [link](http://www.usadellab.org/cms/?page=trimmomatic)
* BWA - Burrow-Wheeler Aligner for short-read alignment, [link](https://github.com/lh3/bwa)
* samtools - Tools for manipulating next-generation sequencing data, [link](https://github.com/samtools/samtools)
* Centrifuge - Classifier for metagenomic sequences, [link](https://ccb.jhu.edu/software/centrifuge/)
* SPAdes -  metaSPAdes: a new versatile de novo metagenomics assembler, [link](https://cab.spbu.ru/software/spades/)
* Prodigal  - Gene Prediction Software, [link](https://github.com/hyattpd/Prodigal)
* cd-hit - sequence clustering, [link](https://github.com/weizhongli/cdhit)
* minimap2 - a versatile pairwise aligner for genomic and spliced nucleotide sequences, [link](https://github.com/lh3/minimap2)
* NCBI-blast+, [link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* RGI - Resistance Gene Identifier, [link](https://github.com/arpcard/rgi)

## Steps to run the analysis

1. prepare the sample file and fastq files
The sample file is included [here](https://github.com/weizhongli/ngomicswf/blob/master/projects/gut-microbiome-resistome-in-children/NGS-samples).
For each sample (e.g. 101A), a folder need to be created using the sample ID as the folder name.
Then, the forward reads (R1) and reverse reads (R2), both in gzipped fastq format, need to be
copied into each sample's folder, as below:
<pre>
101A/R1.fq.gz
101A/R2.fq.gz
103A/R1.fq.gz
103A/R2.fq.gz
103B/R1.fq.gz
103B/R2.fq.gz
105A/R1.fq.gz
105A/R2.fq.gz
106A/R1.fq.gz
106A/R2.fq.gz
106B/R1.fq.gz
106B/R2.fq.gz
...
</pre>


