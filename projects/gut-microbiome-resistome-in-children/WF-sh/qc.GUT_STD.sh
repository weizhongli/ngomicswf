#!/bin/bash
#$ -q default.q
#$ -P 8550
#$ -v PATH
#$ -V


#$ -pe threaded 16

my_host=`hostname`
my_pid=$$
my_core=16
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir qc
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/qc/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/qc/WF.start.date;  fi


#### Trim adapter options
if [ "None" = "NexteraPE" ] 
then
  java -jar /local/ifs2_projdata/9600/projects/BIOINFO/apps/Trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 R1.fq.gz R2.fq.gz qc/R1.fq qc/R1-s.fq qc/R2.fq qc/R2-s.fq \
      ILLUMINACLIP:/local/ifs2_projdata/9600/projects/BIOINFO/apps/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 \
      SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:100 MAXINFO:80:0.5 1>qc/qc.stdout 2>qc/qc.stderr
else
  java -jar /local/ifs2_projdata/9600/projects/BIOINFO/apps/Trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 R1.fq.gz R2.fq.gz qc/R1.fq qc/R1-s.fq qc/R2.fq qc/R2-s.fq \
      SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:100 MAXINFO:80:0.5 1>qc/qc.stdout 2>qc/qc.stderr
fi

perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|GUT_STD|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < qc/R1.fq > qc/R1.fa &
perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|GUT_STD|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < qc/R2.fq > qc/R2.fa &
wait

GZIP=gzip
PIGZ=$(command -v pigz)
if [ $PIGZ ]; then
  GZIP="$PIGZ -p 8"
fi

$GZIP qc/R1.fa &
$GZIP qc/R2.fa &
wait

rm -f qc/R1.fq qc/R2.fq qc/R1-s.fq qc/R2-s.fq

#### qc.txt
NUM_reads_total=$(grep "Input Read Pairs"  qc/qc.stderr | cut -f 4 -d " ")
NUM_reads_pass=$(grep "Input Read Pairs"  qc/qc.stderr | cut -f 7 -d " ")
echo -e "#Reads\tNumber" > qc/qc.txt
echo -e "Total_reads\t$NUM_reads_total" >> qc/qc.txt
echo -e "QC_reads\t$NUM_reads_pass" >> qc/qc.txt


if ! [ -s qc/R1.fa.gz ]; then echo "zero size qc/R1.fa.gz"; exit; fi
if ! [ -s qc/R2.fa.gz ]; then echo "zero size qc/R2.fa.gz"; exit; fi
if ! [ -s qc/qc.txt ]; then echo "zero size qc/qc.txt"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/qc/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=qc host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/qc/WF.cpu

