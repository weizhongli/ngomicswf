#!/bin/bash
#$ -q default.q
#$ -P 8550
#$ -v PATH
#$ -V


#$ -pe threaded 24

my_host=`hostname`
my_pid=$$
my_core=24
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir reads-filtering
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/reads-filtering/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/reads-filtering/WF.start.date;  fi


#### db1 for host filtering
if [ ! "host/GRCh38.fa" = "NONE" ]; then
#### -F 1 for sam-filter-top-pair-or-single.pl means only accept the reads that are mapped in proper pairs
  /local/ifs2_projdata/9600/projects/BIOINFO/apps/bin/bwa mem -t 16 -T 60 /local/ifs2_projdata/9600/projects/BIOINFO/refs/host/GRCh38.fa \
    qc/R1.fa.gz qc/R2.fa.gz | /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/sam-filter-top-pair-or-single.pl -T 60 -F 0 -O reads-filtering/host-hit.ids | \
    /local/ifs2_projdata/9600/projects/BIOINFO/apps/bin/samtools view -b -S - > reads-filtering/host.top.bam
else
  touch reads-filtering/host-hit.ids
fi

#### db2 for rRNA filtering
if [ ! "NONE" = "NONE" ]; then
  /local/ifs2_projdata/9600/projects/BIOINFO/apps/bin/bwa mem -t 16 -T 60  /local/ifs2_projdata/9600/projects/BIOINFO/refs/NONE \
    qc/R1.fa.gz qc/R2.fa.gz | /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/sam-filter-top-pair-or-single.pl -T 60 -F 0 -O reads-filtering/rRNA-hit.ids | \
    /local/ifs2_projdata/9600/projects/BIOINFO/apps/bin/samtools view -b -S - > reads-filtering/rRNA.top.bam
else
  touch reads-filtering/rRNA-hit.ids
fi

cat reads-filtering/host-hit.ids reads-filtering/rRNA-hit.ids | sort | uniq > reads-filtering/filter-hit.ids

GZIP=gzip
PIGZ=$(command -v pigz)
if [ $PIGZ ]; then
  GZIP="$PIGZ -p 8"
fi

if [ -s "reads-filtering/filter-hit.ids" ]; then
  /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_fetch_exclude_ids.pl -i reads-filtering/filter-hit.ids -s  qc/R1.fa.gz -o reads-filtering/filtered-R1.fa &
  /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_fetch_exclude_ids.pl -i reads-filtering/filter-hit.ids -s  qc/R2.fa.gz -o reads-filtering/filtered-R2.fa &
  wait

  $GZIP reads-filtering/filtered-R1.fa &
  $GZIP reads-filtering/filtered-R2.fa &
  wait
else
  #### do nothing, simply link 
  ln -s ../qc/R1.fa.gz reads-filtering/filtered-R1.fa.gz
  ln -s ../qc/R2.fa.gz reads-filtering/filtered-R2.fa.gz
fi

#### filter.txt
NUM_HOST=$(grep -c "." reads-filtering/host-hit.ids)
NUM_RRNA=$(grep -c "." reads-filtering/rRNA-hit.ids)
NUM_HOST_RNA=$(grep -c "." reads-filtering/filter-hit.ids)
echo -e "#Filter\tdb\tNumber"                          > reads-filtering/filter.txt
echo -e "Host\tHost\t$NUM_HOST"                       >> reads-filtering/filter.txt
echo -e "rRNA_tRNA\trRNA_tRNA\t$NUM_RRNA"             >> reads-filtering/filter.txt
echo -e "Both_Host_RNA\tBoth_Host_RNA\t$NUM_HOST_RNA" >> reads-filtering/filter.txt


if ! [ -s reads-filtering/filtered-R1.fa.gz ]; then echo "zero size reads-filtering/filtered-R1.fa.gz"; exit; fi
if ! [ -s reads-filtering/filtered-R2.fa.gz ]; then echo "zero size reads-filtering/filtered-R2.fa.gz"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/reads-filtering/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=reads-filtering host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/reads-filtering/WF.cpu

