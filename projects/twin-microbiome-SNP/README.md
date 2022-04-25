# Microbial Species that Initially Colonize the Human Gut at Birth or in Early Childhood Can Stay in Human Body for Lifetime

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
Python 3.6 or higher is required to run NG-Omics-WF. The following software tools are needed: 
for the the analysis in the individual steps in the workflow. 
* Trimmomatic - A flexible read trimming tool for Illumina NGS data, [link](http://www.usadellab.org/cms/?page=trimmomatic)
* BWA - Burrow-Wheeler Aligner for short-read alignment, [link](https://github.com/lh3/bwa)
* samtools - Tools for manipulating next-generation sequencing data, [link](https://github.com/samtools/samtools)
* SPAdes -  metaSPAdes: a new versatile de novo metagenomics assembler, [link](https://cab.spbu.ru/software/spades/)
* Prodigal  - Gene Prediction Software, [link](https://github.com/hyattpd/Prodigal)
* cd-hit - sequence clustering, [link](https://github.com/weizhongli/cdhit)
* NCBI-blast+, [link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* Varscan2, [link](https://github.com/Jeltje/varscan2)

## Databases needed for this study
see https://github.com/weizhongli/ngomicswf/tree/master/refs-db/microbiome-workflow for
steps to prepare reference databases such as
* microbial genome database for mapping the reads
* human genome database for removing human contamination

## Steps to run the analysis
1. prepare the sample file and fastq files  
The sample file is included [here](./NGS-samples).
For each sample (e.g. S-10042), a folder need to be created using the sample ID as the folder name.
Then, the forward reads (R1) and reverse reads (R2), both in gzipped fastq format, need to be
copied into each sample's folder, as below:
<pre>
S-10042/R1.fq.gz
S-10042/R2.fq.gz
S-10057/R1.fq.gz
S-10057/R2.fq.gz
S-10082/R1.fq.gz
S-10082/R2.fq.gz
S-10115/R1.fq.gz
S-10115/R2.fq.gz
S-10185/R1.fq.gz
S-10185/R2.fq.gz
...
</pre>

2. configure NG-Omics-microbiome-twin-SNP.py  
Edit the the ENV block and NGS-executions block in the top (first 44 lines) of this file.
Also, edit each job block (e.g. NGS_batch_jobs['remove-host']{}), 
so that the path of programs (e.g. bwa) used in that job block are correct. 
See https://github.com/weizhongli/ngomicswf/wiki for info about how to edit this file.

3. run the workflow  
run the following command:
<pre>
nohup python3 PATH_to_ngomicswf/NG-Omics-WF.py3 -i NG-Omics-microbiome-twin-SNP.py -s NGS-samples -t NGS-opts \
    -j taxnonomy,SNPb,SNPb1 &

Here, NG-Omics-microbiome-twin-SNP.py, NGS-samples and NGS-opts are all in this folder
</pre>
Depending on the computer hardware, computer cluster size etc. The workflow may take hours to several days.
The workflow will automatically run all the steps for all the samples.

4. parse the results  
run the following command:
<pre>
PATH_to_ngomicswf/workflow-dev/NG-Omics-microbiome-report-xlsx.sh NGS-samples NGS-samples-dir
</pre>
This will produce a list of tsv files under NGS-samples-dir. These are the files of taxonomy abundance and AMR abundance.

## Citation
DOI:10.1007/s00248-020-01636-0, https://link.springer.com/article/10.1007/s00248-020-01636-0,  https://doi.org/10.1007/s00248-020-01636-0
