#!/bin/bash
#$ -q default.q
#$ -P 8550
#$ -v PATH
#$ -V


#$ -pe threaded 32

my_host=`hostname`
my_pid=$$
my_core=32
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir centrifuge-mapping
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/centrifuge-mapping/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/centrifuge-mapping/WF.start.date;  fi


/local/ifs2_projdata/9600/projects/BIOINFO/apps/centrifuge/centrifuge -1 reads-filtering/filtered-R1.fa.gz -2 reads-filtering/filtered-R2.fa.gz \
  -x /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/centrifuge/ref_full \
  -S centrifuge-mapping/centrifuge-out --report-file centrifuge-mapping/centrifuge-taxon -p 16 -f --out-fmt tab 
/local/ifs2_projdata/9600/projects/BIOINFO/apps/centrifuge/centrifuge-kreport -x /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/centrifuge/ref_full centrifuge-mapping/centrifuge-out > centrifuge-mapping/centrifuge-report

grep -vP "0.0$" centrifuge-mapping/centrifuge-taxon > centrifuge-mapping/centrifuge-taxon.1
(head -n 1 centrifuge-mapping/centrifuge-taxon.1 && tail -n +2 centrifuge-mapping/centrifuge-taxon.1 | sort -grt $'\t' -k 7,7 ) > centrifuge-mapping/centrifuge-taxon.2

# count total number of input reads
NUM_reads=$(zcat reads-filtering/filtered-R1.fa.gz | grep -c "^>")
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/centrifuge-taxon.pl -i centrifuge-mapping/centrifuge-out -j centrifuge-mapping/centrifuge-report \
  -t /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/refseq_genome_taxon.tsv -a /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/centrifuge/ref_full.ann \
  -o centrifuge-mapping/taxon -c 1e-7 -N $NUM_reads -l 60

GZIP=gzip
PIGZ=$(command -v pigz)
if [ $PIGZ ]; then
  GZIP="$PIGZ -p 8"
fi
$GZIP -f centrifuge-mapping/centrifuge-out

cp -p centrifuge-mapping/taxon.superkingdom.txt centrifuge-mapping/taxon.superkingdom-whost.txt
grep "^Host" reads-filtering/filter.txt >> centrifuge-mapping/taxon.superkingdom-whost.txt
grep "^rRNA" reads-filtering/filter.txt >> centrifuge-mapping/taxon.superkingdom-whost.txt


if ! [ -s centrifuge-mapping/centrifuge-out.gz ]; then echo "zero size centrifuge-mapping/centrifuge-out.gz"; exit; fi
if ! [ -s centrifuge-mapping/centrifuge-taxon ]; then echo "zero size centrifuge-mapping/centrifuge-taxon"; exit; fi
if ! [ -s centrifuge-mapping/centrifuge-report ]; then echo "zero size centrifuge-mapping/centrifuge-report"; exit; fi

date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/centrifuge-mapping/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=centrifuge-mapping host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/centrifuge-mapping/WF.cpu

