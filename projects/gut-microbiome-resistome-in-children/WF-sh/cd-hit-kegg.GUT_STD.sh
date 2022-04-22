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
mkdir cd-hit-kegg
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/cd-hit-kegg/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/cd-hit-kegg/WF.start.date;  fi

/local/ifs2_projdata/9600/projects/BIOINFO/apps/cd-hit/cd-hit-2d -i /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/keggf -i2 ORF-prediction/ORF.faa -o cd-hit-kegg/out \
  -c 0.75 -n 5 -d 0 -g 1 -G 0 -aS 0.9 -A 60 -aL 0.25 -T 16 -M 32000 > cd-hit-kegg/out.log
/local/ifs2_projdata/9600/projects/BIOINFO/apps/cd-hit/clstr_select.pl 2 99999999 < cd-hit-kegg/out.clstr > cd-hit-kegg/out.clstr.1
mv -f cd-hit-kegg/out.clstr.1 cd-hit-kegg/out.clstr
/local/ifs2_projdata/9600/projects/BIOINFO/apps/cd-hit/cd-hit-clstr_2_blm8.pl < cd-hit-kegg/out.clstr > cd-hit-kegg/out.bl

mkdir cd-hit-kegg/orf-split
/local/ifs2_projdata/9600/projects/BIOINFO/apps/cd-hit/cd-hit-div.pl cd-hit-kegg/out cd-hit-kegg/orf-split/split 64


date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/cd-hit-kegg/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=cd-hit-kegg host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/cd-hit-kegg/WF.cpu

