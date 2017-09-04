=============================== NG-Omics-WF ==================================
Workflow tools for next generation genomics, metagenomics, RNA-seq 
and other type of omics data analyiss, 

Software originally developed since 2010 by Weizhong Li at UCSD
                                              currently at JCVI

https://github.com/weizhongli/ngomicswf       liwz@sdsc.edu
==============================================================================


NG-Omics-WF is a workflow tool to automatically run a bioinformatic pipeline
for multiple datasets under generic Linux computer or Linux cluster environment.
The tool was oritinally implemented with supports from the Human Microbiome Project
(HMP) and CAMERA project at UCSD for the analysis of next generation metagenomic 
sequence data, especially human microbiome data.

The tool was further improved at UCSD and then JCVI to serve as the backend workflow 
engine in the USDA supported RNA-seq portal development project. 

The tool was originally written in Perl, it is now re-written in Python. Python 2.7 
is needed to run the tool. This program need to run under generic Linux system, either
on a standalone computer or a computer cluster that support queue system, such as 
open grid engine (formally sun grid engine, SGE).

Directory workflow-examples has several workflow examples, which can be directly used
or after some re-configuration.


 
