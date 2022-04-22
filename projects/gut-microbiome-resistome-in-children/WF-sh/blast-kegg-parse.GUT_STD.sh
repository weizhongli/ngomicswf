#!/bin/bash
##$ -q default.q
#$ -q centos7.q
#$ -P 0804
#$ -v PATH
#$ -V


#$ -pe threaded 2

my_host=`hostname`
my_pid=$$
my_core=2
my_queue=qsub_1
my_time_start=`date +%s`

cd /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD
mkdir blast-kegg-parse
if ! [ -f /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg-parse/WF.start.date ]; then date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg-parse/WF.start.date;  fi

ln -s ../../cd-hit-kegg/out.bl blast-kegg/blast/out.bl

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_prebin.pl -i blast-kegg/blast -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/keggf.clstr.ann \
  -a ORF-prediction/ORF.faa -o blast-kegg-parse/ORF -t /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/keggf_taxon.txt -s minimap-binning-2020/assembly-bin \
  -x /local/ifs2_projdata/9600/projects/BIOINFO/refs/ref-genomes-2020-0622/refseq_genome_taxon.tsv -p blast-kegg-parse/scaffold-ann.txt -d assembly/assembly/scaffold-cov

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00001-biome.keg -o blast-kegg-parse/ORF-ann-kegg-pathway     -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-edit.keg  -o blast-kegg-parse/ORF-ann-kegg-module      -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01000.keg       -o blast-kegg-parse/ORF-ann-kegg-EC          -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko02000.keg       -o blast-kegg-parse/ORF-ann-kegg-transporter -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01504.keg       -o blast-kegg-parse/ORF-ann-kegg-AMR         -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00001-biome.keg -o blast-kegg-parse/sp-pathway     -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-edit.keg  -o blast-kegg-parse/sp-module      -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01000.keg       -o blast-kegg-parse/sp-EC          -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko02000.keg       -o blast-kegg-parse/sp-transporter -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01504.keg       -o blast-kegg-parse/sp-AMR         -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg

/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00001-biome.keg -o blast-kegg-parse/ge-pathway     -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg -K 8
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-edit.keg  -o blast-kegg-parse/ge-module      -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg -K 8
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01000.keg       -o blast-kegg-parse/ge-EC          -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg -K 8
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko02000.keg       -o blast-kegg-parse/ge-transporter -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg -K 8
/local/ifs2_projdata/9600/projects/BIOINFO/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i blast-kegg-parse/ORF-ann.txt -d ORF-prediction/ORF-cov -k /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko01504.keg       -o blast-kegg-parse/ge-AMR         -r /local/ifs2_projdata/9600/projects/BIOINFO/refs/kegg/ko00002-M00178.keg -K 8



date +%s > /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg-parse/WF.complete.date
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample=GUT_STD job=blast-kegg-parse host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> /local/ifs2_projdata/0804/projects/amr_children/Novaseq/GUT_STD/blast-kegg-parse/WF.cpu

