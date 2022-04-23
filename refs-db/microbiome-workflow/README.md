This folder include guide for preparing several reference databases to be used
in the microibome analysis workflow. e.g.
https://github.com/weizhongli/ngomicswf/blob/master/workflow-dev/NG-Omics-microbiome-2019.py

## build microbial reference genome database for mapping
This database is for mapping the reads using Centrifuge to get taxonomy abundance.
It is very complicated and time consuming to build a good and comprehensive
microbial reference genome database for mapping metagenomic reads for taxonomy
profiling.

Only highly experienced bioinformatics expertise with deep knowledge of NCBI
resources can perform the steps below. Unexperienced users will likely fail
during the process.

1. select representative microbial genomes from NCBI
 
The microbial genomes deposited at NCBI are highly redundant. Using all the genomes as
reference database is not only impossible (due to the huge number of genomes), but also 
unnecessary. So it is important to select a subset of representative microbial genomes of
high quality. This process can be done in different ways, but it is very time consuming.
Here, I just give a pre-selected list of 29324 genomes. see file
https://github.com/weizhongli/ngomicswf/blob/master/refs-db/microbiome-workflow/refdb-2020-rep-genomes-compact.tsv.gz.
The process of how these 29324 genomes were selected will be explained elsewhere.


Here are the first few line of this file:
<pre>
accession	taxid	organism_name	superkingdom	species	species_ti	genus	genus_ti
GCF_009729015.1	2283	Acidianus ambivalens	Archaea	Acidianus ambivalens	2283	Acidianus	12914
GCF_003201835.2	41673	Acidianus brierleyi	Archaea	Acidianus brierleyi	41673	Acidianus	12914
GCF_000213215.1	933801	Acidianus hospitalis W1	Archaea	Acidianus hospitalis	563177	Acidianus	12914
GCF_009729545.1	12915	Acidianus infernus	Archaea	Acidianus infernus	12915	Acidianus	12914

Here accession is the NCBI accession number for the genome sequence of the a species or strain
</pre>

2. download the genomic sequences for all these 29324 genomes 

Please write a script to download the genome sequences and reformat the headers of
the fasta file like the following, 
<pre>
>acc|NZ_KN266207.1|taxid|1209072|sptid|39650|getid|10 Cellvibrio mixtus subsp. mixtus J3-8 Scaffold1, whole genome shotgun sequence
GAAGTTAGATTATGAGGCTTTCGGGTGTCTAGTAAAGTCTGGGAAGTCCAGATGAGCTAAAGTTACGAACCCGCTTTACT
AACATGCAATTGCATATTTAATGCTTTTAGCACTTTAAGAATTGTCGCAAAGCTGGGGTTTCCTTCGCCGGATAACGCTT
TATAAAGGCTCTCGCGGCCGAGGTCACAATCTTTGGCTAACTGGGTCATGCCTCTGGCTTTAGCAACATTGCCTAATGCT
TTGGTGATAAACGCTGCATCATCGCCAGCTTCTTCAAAACAAGCATCTAGATAGAGAGCAATATCCTCTTCAGATTTTAA
...
</pre>
Then concatenate the individual 29324 fasta files into a single fasta file, and name
the file as ref-genome-all.fna

You will need https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
to find the download path of genomic sequence based on the accession number. e.g.
ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/985/GCF_000001985.1_JCVI-PMFA1-2.0/GCF_000001985.1_JCVI-PMFA1-2.0_genomic.fna.gz

3. download NCBI taxonomy data 

Download https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz and unzip it

4. prepare an accession to taxid mapping file 

Please write a script to process the headers of ref-genome-all.fna to generate a
text file like below with two columns like below:
<pre>
acc|NZ_KN266207.1	1209072
acc|NZ_KN266208.1	1209072
acc|NZ_KN266209.1	1209072
acc|NZ_ALBT01000032.1	1209072
acc|NZ_ALBT01000033.1	1209072
</pre>

5. generate a taxonomy data file: refseq_genome_taxon.tsv 

Below is a shell script template
<pre>
#!/bin/bash
Path_to_ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl  -o refseq_genome_taxon.tsv.s0 \
  -f ref-genome-all.fna -d /Path_to_NCBI_taxonomy_folder
cat refseq_genome_taxon.tsv.s0 | awk -F '\t' '{ if ($2 == "toprank" || $2=="rank") {print} }' > refseq_genome_taxon.tsv
rm refseq_genome_taxon.tsv.s0
</pre>

6. index the ref-genome-all.fna by program centrifuge 

Below is a shell script template, please run this script on a computer with
at least 256 GB RAM and 32 cores.
<pre>
#!/bin/bash
/Path_to_centrifuge/centrifuge-build -p 16 --bmax 792448239 --dcv 1024 \
  --conversion-table accession_to_taxid_mapping_file \
  --taxonomy-tree /Path_to_NCBI_taxonomy_folder/nodes.dmp \
  --name-table /Path_to_NCBI_taxonomy_folder/names.dmp ref_full.fna ref_full 
/Path_to_centrifuge/centrifuge-inspect --size-table ref_full > ref_full.taxid.len
</pre>

7. index the ref-genome-all.fna by bwa and samtools 

Below is a shell script template, please run this script on a computer with
at least 256 GB RAM and 32 cores.
<pre>
#!/bin/bash
/Path_to_bwa/bwa index -a bwtsw -b 80000000 -p ref_full ref_full.fna
/Path_to_samtools/samtools faidx ref_full.fna
</pre>

After step 6 and 7, the microbial reference genome database is ready. You should have a
list of files like below:
<pre>
-rw-rw-r-- 1 wli 9600-grp 104363722971 Jun 25  2020 ref_full.fna
-rw-rw-r-- 1 wli 9600-grp     51869830 Jun 25  2020 ref_full.3.cf
-rw-rw-r-- 1 wli 9600-grp     20046636 Jun 26  2020 ref_full.4.cf
-rw-rw-r-- 1 wli 9600-grp  34465946584 Jun 26  2020 ref_full.1.cf
-rw-rw-r-- 1 wli 9600-grp  25609902912 Jun 26  2020 ref_full.2.cf
-rw-rw-r-- 1 wli 9600-grp       434616 Jun 26  2020 ref_full.taxid.len
-rw-rw-r-- 1 wli 9600-grp 102839901184 Jun 27  2020 ref_full.bwt
-rw-rw-r-- 1 wli 9600-grp  25709975277 Jun 27  2020 ref_full.pac
-rw-rw-r-- 1 wli 9600-grp    271316745 Jun 27  2020 ref_full.ann
-rw-rw-r-- 1 wli 9600-grp     15920113 Jun 27  2020 ref_full.amb
-rw-rw-r-- 1 wli 9600-grp  51419950600 Jun 28  2020 ref_full.sa
-rw-rw-r-- 1 wli 9600-grp    138162886 Jun 28  2020 ref_full.fna.fai
-rw-rw-r-- 1 wli 9600-grp     14416730 Jun 22  2020 refseq_genome_taxon.tsv
</pre>

## download human genome database
The human genome database, fasta file and index files formated by bwa is
used to map the reads and filter the reads that from human contamination.
Please download GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz from [NCBI](https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/) and unzip it.


