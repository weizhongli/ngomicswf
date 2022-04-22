#!/bin/bash
##$ -q default.q
#$ -q centos7.q
#$ -P 0804
#$ -v PATH
#$ -V


#$ -pe threaded 24

my_host=`hostname`
my_pid=$$
my_core=24
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir assembly
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/assembly/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/assembly/WF.start.date;  fi

if [ "spade" = "idba-ud" ]
then
  /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/PE-2file-to-1file.pl -i reads-filtering/filtered-R1.fa.gz,reads-filtering/filtered-R2.fa.gz -r 0 > assembly/input.fa
  /local/ifs2_projdata/9600/projects/BIOINFO/apps/bin/idba_ud -r assembly/input.fa -o assembly/assembly --mink=50 --maxk=80 --step=10 --num_threads=16 --min_contig=250 
  rm -f assembly/input.fa assembly/assembly/kmer assembly/assembly/local* assembly/assembly/contig* assembly/assembly/graph*
  
elif [ "spade" = "spade" ]
then
  python /local/ifs2_projdata/9600/projects/BIOINFO/apps/SPAdes/bin/spades.py -1 reads-filtering/filtered-R1.fa.gz -2 reads-filtering/filtered-R2.fa.gz --meta --only-assembler -o assembly/assembly -t 16
  mv assembly/assembly/scaffolds.fasta assembly/assembly/scaffold.fa
  rm -rf assembly/assembly/K*
  rm -rf assembly/assembly/tmp

  gzip assembly/assembly/*fasta &
  gzip assembly/assembly/*fastg &
  gzip assembly/assembly/*gfa   &
  gzip assembly/assembly/*paths &
  wait

else
  echo "not defined assembly method"
  exit 1  
fi


## ensure to have >SAMPLE_name|scaffold|n header format
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_filter_short_seq.pl -i assembly/assembly/scaffold.fa -c 250 -o - -a len | \
  /local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/fasta_rename.pl -i - -s "GUT_STD|scaffold|" -b 0 -o assembly/assembly/scaffold-new.fa
mv -f assembly/assembly/scaffold-new.fa assembly/assembly/scaffold.fa

## depth of coverage

perl -e 'while(<>){ if ($_ =~ /^>(\S+)/) { $id=$1;  if ($_ =~ /_cov_([\d\.]+)/) { print "$id\t$1\n";} } }' < assembly/assembly/scaffold.fa > assembly/assembly/scaffold-cov

## $NGS_bin_dir/bwa index -a bwtsw -p assembly/scaffold assembly/scaffold.fa
## insert coverage 

if ! [ -s assembly/assembly/scaffold.fa ]; then echo "zero size assembly/assembly/scaffold.fa"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/assembly/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=assembly host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/assembly/WF.cpu

