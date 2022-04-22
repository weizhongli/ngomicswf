#!/bin/bash
##$ -q default.q
#$ -q centos7.q
#$ -P 0804
#$ -v PATH
#$ -V


#$ -pe threaded 16

my_host=`hostname`
my_pid=$$
my_core=16
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir minimap-binning-2020
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/minimap-binning-2020/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/minimap-binning-2020/WF.start.date;  fi


#### print taxids with depth > cutoff
grep -v "^#" centrifuge-mapping/taxon.toprank.txt | gawk -F '\t' '($12+0.0) > 0.01 {print $1}' > minimap-binning-2020/high_cov.taxids
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_idx_fetch_by_ids.pl -i minimap-binning-2020/high_cov.taxids -s /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/bwa/ref_full.fna \
  -a taxid -h /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/bwa/ref_full.fna.header -o minimap-binning-2020/high_cov_ref.fna

/local/ifs2_projdata/9600/projects/BIOINFO/apps/minimap2/minimap2 minimap-binning-2020/high_cov_ref.fna assembly/assembly/scaffold.fa -N 20 -t 8 > minimap-binning-2020/hits.paf
rm -f minimap-binning-2020/high_cov_ref.fna

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/JCVI/minimap-binning-taxon.pl -i minimap-binning-2020/hits.paf -s assembly/assembly/scaffold.fa \
  -t /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/refseq_genome_taxon.tsv -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/bwa/ref_full.taxid.len -o minimap-binning-2020/assembly-bin

if ! [ -s minimap-binning-2020/assembly-bin ]; then echo "zero size minimap-binning-2020/assembly-bin"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/minimap-binning-2020/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=minimap-binning-2020 host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/minimap-binning-2020/WF.cpu

