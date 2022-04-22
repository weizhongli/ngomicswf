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
mkdir blast-kegg
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg/WF.start.date;  fi


#### skip this - only use cd-hit-kegg's result
if [ "Skip" = "Skip" ] 
then
  mkdir blast-kegg/blast
  echo "Skip" >> blast-kegg/skip.txt
else
  for i in `seq 1 4`
    do /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_batch_run_dir.pl --INDIR1=cd-hit-kegg/orf-split --OUTDIR1=blast-kegg/blast --CPU=blast-kegg/WF.cpu /local/ifs2_projdata/9600/projects/BIOINFO/apps/blast+/bin/blastp  -query {INDIR1} -out {OUTDIR1} \
    -db /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/keggf -evalue 1e-6 -num_threads 4 -num_alignments 5 -outfmt 6 -seg yes &
  done
  wait
fi



date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=blast-kegg host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg/WF.cpu

