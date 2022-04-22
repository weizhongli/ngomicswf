#!/bin/bash
##$ -q default.q
#$ -q centos7.q
#$ -P 0804
#$ -v PATH
#$ -V


#$ -pe threaded 4

my_host=`hostname`
my_pid=$$
my_core=4
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir RGI
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/RGI/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/RGI/WF.start.date;  fi


#### #source /local/ifs2_projdata/9600/projects/BIOINFO/local-py3/bin/activate
#### #export LD_LIBRARY_PATH=/usr/local/packages/gcc-4.7.1/lib64
#### source /local/ifs2_projdata/0804/projects/amr_children/IFX/anaconda3.source.me
#### export PATH=/local/ifs2_projdata/9600/projects/BIOINFO/apps/blast+/bin/:$PATH

#### rgi main --input_sequence ORF-prediction/ORF.faa --output_file RGI/rgi.out --input_type protein
#### rm -f RGI/ORF.faa.temp.*

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/JCVI/ann_parse_RGI.pl -i RGI/rgi.out.txt -o RGI/rgi.out -a ORF-prediction/ORF-cov -b minimap-binning-2020/assembly-bin
#### delete below no use
#### /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/JCVI/ann_parse_RGI.pl -i RGI/rgi.out.txt -o RGI/count-rgi.out -b minimap-binning-2020/assembly-bin


if ! [ -s RGI/rgi.out.txt ]; then echo "zero size RGI/rgi.out.txt"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/RGI/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=RGI host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/RGI/WF.cpu

