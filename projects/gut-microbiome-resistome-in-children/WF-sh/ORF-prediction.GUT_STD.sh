#!/bin/bash
##$ -q default.q
#$ -q centos7.q
#$ -P 0804
#$ -v PATH
#$ -V


#$ -pe threaded 8

my_host=`hostname`
my_pid=$$
my_core=8
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir ORF-prediction
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/ORF-prediction/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/ORF-prediction/WF.start.date;  fi

if [ "prodigal" = "metagene" ]
then
  /local/ifs2_projdata/9600/projects/BIOINFO/apps/metagene/metagene_run.pl assembly/assembly/scaffold.fa ORF-prediction/ORF.faa

elif [ "prodigal" = "prodigal" ]
then
  /local/ifs2_projdata/9600/projects/BIOINFO/apps/Prodigal/prodigal -i assembly/assembly/scaffold.fa -p meta -o ORF-prediction/ORF.gff -f gff -a ORF-prediction/ORF.faa -d ORF-prediction/ORF.fna

else
  echo "undefined ORF-prediction method"
  exit 1  
fi

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_filter_short_seq.pl -i ORF-prediction/ORF.faa -c 20 -o ORF-prediction/ORF-new.faa
mv  ORF-prediction/ORF-new.faa ORF-prediction/ORF.faa

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/assembly-cov-pass-to-orf.pl -i ORF-prediction/ORF.faa -d assembly/assembly/scaffold-cov -o ORF-prediction/ORF-cov

if ! [ -s ORF-prediction/ORF.faa ]; then echo "zero size ORF-prediction/ORF.faa"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/ORF-prediction/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=ORF-prediction host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/ORF-prediction/WF.cpu

