<PRE>
=============================== NG-Omics-WF ==================================
Workflow tools for next generation genomics, metagenomics, RNA-seq 
and other type of omic data analysis

Software originally developed since 2010 by Weizhong Li at UCSD
                                              currently at JCVI

https://github.com/weizhongli/ngomicswf       liwz@sdsc.edu
==============================================================================
</PRE>
  
NG-Omics-WF is a workflow tool to automatically run a bioinformatic pipeline
for multiple datasets under generic Linux computer or Linux cluster environment.
The tool was oritinally implemented since 2010 with supports from the Human Microbiome Project
(HMP) and CAMERA project at UCSD for the analysis of next generation metagenomic 
sequence data. The tool was further improved at UCSD and then at JCVI as the backend workflow 
engine for many other projects for the analysis of metagenomic data, RNA-seq data, genomic data
and other omic data. 

The tool was originally written in Perl, it was re-written in Python 2.7 and in Python 3. 
The Perl and Python 2.7 version were no longer supported (though they still work) since 2019. 
Currently, Python 3.6 or later version is needed to run the tool. 
This program need to run under generic Linux system, either
on a standalone computer or a computer cluster that support queue system, such as 
open grid engine (formally sun grid engine, SGE).

Directory workflow-examples has several workflow examples, which can be directly used
or after some re-configuration.

The detailed documents are available at https://github.com/weizhongli/ngomicswf/wiki. below is 
a brief user's guide: 

## System Requirements
NG-Omics-WF works on any generic Linux computer or Linux cluster environment. 
Any common Linux systems, such as Ubuntu, CentOS, Fedora, Debian, are good. 
The hardware requirements depend on the size and analysis type of input datasets. For example,
for a typical metagenomic datasets with several GB of sequence data / sample, computers 
with >=64GB RAM and >=32 cores are recommanded. 

## Software Requirements
Python 3.6 or higher is required to run NG-Omics-WF. But many other software tools are needed 
for the the analysis in the individual steps in the workflow. The required tools depends on the
type of input data and need of the analysis. Here are some examples:

### Software Requirements for typical metagenomic analysis
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

Optional
* Hmmer3 - Biological sequence analysis using profile hidden Markov models, [link](http://hmmer.org/download.html)
* RGI - Resistance Gene Identifier, [link](https://github.com/arpcard/rgi)
