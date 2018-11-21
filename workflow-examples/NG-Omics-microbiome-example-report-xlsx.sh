#!/bin/sh

# $1 input NGS-sample file
# $2 output directory

if [ $# -lt 2 ]
  then
    echo "Usage:"; echo
    echo "$0 NGS-sample-file output-dir"; echo
    exit
fi

mkdir -p $2

s_dir=`dirname $0`
s_dir=`dirname $s_dir`
p_pl="$s_dir/NGS-tools/NGS-table-merge.pl -s $1 -t 0 -p 0"
xlsx_py="$s_dir/NGS-tools/NGS-tsv_2_xlsx.py"


  $p_pl -f cdd-parse/cog.txt        -o $2/cog_abundance.tsv               -c 0 -i 0 -v 4 -a 1,5
  $p_pl -f cdd-parse/cog-class.txt  -o $2/cog_class.tsv                   -c 0 -i 0 -v 3 -a 4
  $p_pl -f cdd-parse/kog.txt        -o $2/kog_abundance.tsv               -c 0 -i 0 -v 4 -a 1,5
  $p_pl -f cdd-parse/kog-class.txt  -o $2/kog_class.tsv                   -c 0 -i 0 -v 3 -a 4
  
  $p_pl -f pfam-parse/pfam-ann.txt  -o $2/pfam_abundance.tsv               -c 0 -i 0 -v 3 -a 4,5
  
  $p_pl -f qc/qc.txt -o $2/qc.tsv -c 0 -i 0 -v 1
  
cd $2

  $xlsx_py -i qc.tsv,cog_abundance.tsv,cog_class.tsv,kog_abundance.tsv,kog_class.tsv,pfam_abundance.tsv -o annotation.xlsx

