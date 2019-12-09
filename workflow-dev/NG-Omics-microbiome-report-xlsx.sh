#!/bin/sh

# $1 input NGS-sample file
# $2 output directory

if [ $# -lt 2 ]
  then
    echo "Usage:"; echo
    echo "$0 NGS-sample-file output-dir [taxonomy_dir function_dir optional]"; echo
    exit
fi

mkdir -p $2

s_dir=`dirname $0`
s_dir=`dirname $s_dir`
p_pl="$s_dir/NGS-tools/NGS-table-merge.pl -s $1 -t 0 -p 0"
xlsx_py="$s_dir/NGS-tools/NGS-tsv_2_xlsx.py"

if [ -n "$3" ]; then
  TAXDIR=$3
else
  TAXDIR=taxonomy
fi

if [ -n "$4" ]; then
  FUNDIR=$4
else
  FUNDIR=blast-kegg-parse
fi

#### first sample
SAM1=$(head -n 1 $1 | cut -f 1)

if [ -s $SAM1/$FUNDIR/ORF-ann-kegg-module-full-des ]
then
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-full-des      -o $2/module-KO-depth.tsv               -c 0 -i 0 -v 8  -a 1,2,3,4,5,6,7
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-full-des      -o $2/module-KO-depth-adjusted.tsv      -c 0 -i 0 -v 9  -a 1,2,3,4,5,6,7
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-full-des      -o $2/module-KO-R.depth.tsv             -c 0 -i 0 -v 12 -a 1,2,3,4,5,6,7
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-full-des      -o $2/module-KO-R.depth-adjusted.tsv    -c 0 -i 0 -v 13 -a 1,2,3,4,5,6,7
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-D             -o $2/module-depth.tsv                  -c 0 -i 3 -v 8  -a 0,1,2,4
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-D             -o $2/module-depth-adjusted.tsv         -c 0 -i 3 -v 9  -a 0,1,2,4
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-D             -o $2/module-R.depth.tsv                -c 0 -i 3 -v 10 -a 0,1,2,4
  $p_pl -f $FUNDIR/ORF-ann-kegg-module-D             -o $2/module-R.depth-adjusted.tsv       -c 0 -i 3 -v 11 -a 0,1,2,4
  
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-full-des     -o $2/pathway-KO-depth.tsv                -c 0 -i 0 -v 7  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-full-des     -o $2/pathway-KO-depth-adjusted.tsv       -c 0 -i 0 -v 8  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-full-des     -o $2/pathway-KO-R.depth.tsv              -c 0 -i 0 -v 11 -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-full-des     -o $2/pathway-KO-R.depth-adjusted.tsv     -c 0 -i 0 -v 12 -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-C            -o $2/pathway-depth.tsv                   -c 0 -i 2 -v 7  -a 0,1,3
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-C            -o $2/pathway-depth-adjusted.tsv          -c 0 -i 2 -v 8  -a 0,1,3
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-C            -o $2/pathway-R.depth.tsv                 -c 0 -i 2 -v 9  -a 0,1,3
  $p_pl -f $FUNDIR/ORF-ann-kegg-pathway-C            -o $2/pathway-R.depth-adjusted.tsv        -c 0 -i 2 -v 10 -a 0,1,3
  
  $p_pl -f $FUNDIR/ORF-ann-kegg-transporter-full-des -o $2/transporter-KO-depth.tsv            -c 0 -i 0 -v 7  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-transporter-full-des -o $2/transporter-KO-depth-adjusted.tsv   -c 0 -i 0 -v 8  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-transporter-full-des -o $2/transporter-KO-R.depth.tsv          -c 0 -i 0 -v 11 -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-transporter-full-des -o $2/transporter-KO-R.depth-adjusted.tsv -c 0 -i 0 -v 12 -a 1,2,3,4,5,6
  
  $p_pl -f $FUNDIR/ORF-ann-kegg-EC-full-des          -o $2/EC-KO-depth.tsv                     -c 0 -i 0 -v 7  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-EC-full-des          -o $2/EC-KO-depth-adjusted.tsv            -c 0 -i 0 -v 8  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-EC-full-des          -o $2/EC-KO-R.depth.tsv                   -c 0 -i 0 -v 11 -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-EC-full-des          -o $2/EC-KO-R.depth-adjusted.tsv          -c 0 -i 0 -v 12 -a 1,2,3,4,5,6
  
  $p_pl -f $FUNDIR/ORF-ann-kegg-AMR-full-des         -o $2/AMR-KO-depth.tsv                    -c 0 -i 0 -v 7  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-AMR-full-des         -o $2/AMR-KO-depth-adjusted.tsv           -c 0 -i 0 -v 8  -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-AMR-full-des         -o $2/AMR-KO-R.depth.tsv                  -c 0 -i 0 -v 11 -a 1,2,3,4,5,6
  $p_pl -f $FUNDIR/ORF-ann-kegg-AMR-full-des         -o $2/AMR-KO-R.depth-adjusted.tsv         -c 0 -i 0 -v 12 -a 1,2,3,4,5,6
fi
  
  
if [ -s $SAM1/ann-summary/assembly-orf-summary.txt ]
then
  $p_pl -f ann-summary/assembly-orf-summary.txt -o $2/assembly-num_scaffolds.txt  -c 0 -i 0 -v 1
  $p_pl -f ann-summary/assembly-orf-summary.txt -o $2/assembly-len_scaffolds.txt  -c 0 -i 0 -v 2
  $p_pl -f ann-summary/assembly-orf-summary.txt -o $2/assembly-num_ORFs.txt       -c 0 -i 0 -v 3
  
  $p_pl -f ann-summary/assembly-taxon-genome  -o $2/assembly-genome-scaffolds.tsv -c 0 -i 0 -v 10 -a 1,2,3,4,5,6,7,8,9
  $p_pl -f ann-summary/assembly-taxon-genome  -o $2/assembly-genome-length.tsv    -c 0 -i 0 -v 11 -a 1,2,3,4,5,6,7,8,9
  $p_pl -f ann-summary/assembly-taxon-genome  -o $2/assembly-genome-coverage.tsv  -c 0 -i 0 -v 12 -a 1,2,3,4,5,6,7,8,9
  $p_pl -f ann-summary/assembly-taxon-genome  -o $2/assembly-genome-depth.tsv     -c 0 -i 0 -v 13 -a 1,2,3,4,5,6,7,8,9
  $p_pl -f ann-summary/assembly-taxon-genome  -o $2/assembly-genome-ORFs.tsv      -c 0 -i 0 -v 14 -a 1,2,3,4,5,6,7,8,9
  
  $p_pl -f ann-summary/assembly-taxon-species -o $2/assembly-species-scaffolds.tsv -c 0 -i 0 -v 9  -a 1,2,3,4,5,6,7,8
  $p_pl -f ann-summary/assembly-taxon-species -o $2/assembly-species-length.tsv    -c 0 -i 0 -v 10 -a 1,2,3,4,5,6,7,8
  $p_pl -f ann-summary/assembly-taxon-species -o $2/assembly-species-coverage.tsv  -c 0 -i 0 -v 11 -a 1,2,3,4,5,6,7,8
  $p_pl -f ann-summary/assembly-taxon-species -o $2/assembly-species-depth.tsv     -c 0 -i 0 -v 12 -a 1,2,3,4,5,6,7,8
  $p_pl -f ann-summary/assembly-taxon-species -o $2/assembly-species-ORFs.tsv      -c 0 -i 0 -v 13 -a 1,2,3,4,5,6,7,8
fi  

if [ -s $SAM1/$TAXDIR/taxon.superkingdom.txt ]
then
  $p_pl -f qc/qc.txt -o $2/qc.tsv -c 0 -i 0 -v 1
  $p_pl -f $TAXDIR/taxon.superkingdom-whost.txt    -o $2/taxon-superkingdom-whost.tsv -c 0      -i 0 -a 1 -v 2

  $p_pl -f $TAXDIR/taxon.superkingdom.txt    -o $2/taxon-superkingdom.tsv       -c 0      -i 0 -a 1 -v 2
  $p_pl -f $TAXDIR/taxon.toprank.txt         -o $2/taxon-toprank.tsv            -c 1e-4   -i 0 -a 1,2,3,4,5,6,7,8,9 -v 10
  $p_pl -f $TAXDIR/taxon.species.txt         -o $2/taxon-species.tsv            -c 1e-4   -i 0 -a 1,2,3,4,5,6,7,8 -v 9
  $p_pl -f $TAXDIR/taxon.genus.txt           -o $2/taxon-genus.tsv              -c 1e-4   -i 0 -a 1,2,3,4,5,6,7   -v 8
  $p_pl -f $TAXDIR/taxon.family.txt          -o $2/taxon-family.tsv             -c 1e-4   -i 0 -a 1,2,3,4,5,6     -v 7
  $p_pl -f $TAXDIR/taxon.order.txt           -o $2/taxon-order.tsv              -c 1e-4   -i 0 -a 1,2,3,4,5       -v 6
  $p_pl -f $TAXDIR/taxon.class.txt           -o $2/taxon-class.tsv              -c 1e-4   -i 0 -a 1,2,3,4         -v 5
  $p_pl -f $TAXDIR/taxon.phylum.txt          -o $2/taxon-phylum.tsv             -c 1e-4   -i 0 -a 1,2,3           -v 4
  
  $p_pl -f $TAXDIR/taxon.toprank.txt         -o $2/taxon-depth-toprank.tsv      -c 1e-4   -i 0 -a 1,2,3,4,5,6,7,8,9 -v 11
  $p_pl -f $TAXDIR/taxon.species.txt         -o $2/taxon-depth-species.tsv      -c 1e-4   -i 0 -a 1,2,3,4,5,6,7,8 -v 10
  $p_pl -f $TAXDIR/taxon.genus.txt           -o $2/taxon-depth-genus.tsv        -c 1e-4   -i 0 -a 1,2,3,4,5,6,7   -v 9
  $p_pl -f $TAXDIR/taxon.family.txt          -o $2/taxon-depth-family.tsv       -c 1e-4   -i 0 -a 1,2,3,4,5,6     -v 8
  $p_pl -f $TAXDIR/taxon.order.txt           -o $2/taxon-depth-order.tsv        -c 1e-4   -i 0 -a 1,2,3,4,5       -v 7
  $p_pl -f $TAXDIR/taxon.class.txt           -o $2/taxon-depth-class.tsv        -c 1e-4   -i 0 -a 1,2,3,4         -v 6
  $p_pl -f $TAXDIR/taxon.phylum.txt          -o $2/taxon-depth-phylum.tsv       -c 1e-4   -i 0 -a 1,2,3           -v 5
  
  $p_pl -f $TAXDIR/taxon.toprank.txt         -o $2/taxon-reads-toprank.tsv      -c 1      -i 0 -a 1,2,3,4,5,6,7,8,9 -v 12
  $p_pl -f $TAXDIR/taxon.species.txt         -o $2/taxon-reads-species.tsv      -c 1      -i 0 -a 1,2,3,4,5,6,7,8 -v 11
  $p_pl -f $TAXDIR/taxon.genus.txt           -o $2/taxon-reads-genus.tsv        -c 1      -i 0 -a 1,2,3,4,5,6,7   -v 10
  $p_pl -f $TAXDIR/taxon.family.txt          -o $2/taxon-reads-family.tsv       -c 1      -i 0 -a 1,2,3,4,5,6     -v 9
  $p_pl -f $TAXDIR/taxon.order.txt           -o $2/taxon-reads-order.tsv        -c 1      -i 0 -a 1,2,3,4,5       -v 8
  $p_pl -f $TAXDIR/taxon.class.txt           -o $2/taxon-reads-class.tsv        -c 1      -i 0 -a 1,2,3,4         -v 7
  $p_pl -f $TAXDIR/taxon.phylum.txt          -o $2/taxon-reads-phylum.tsv       -c 1      -i 0 -a 1,2,3           -v 6
fi  
  
#### chdir
cd $2

if [ -s module-KO-depth.tsv ]
then
  $xlsx_py -i module-KO-depth.tsv,module-KO-depth-adjusted.tsv,module-KO-R.depth.tsv,module-KO-R.depth-adjusted.tsv,module-depth.tsv,module-depth-adjusted.tsv,module-R.depth.tsv,module-R.depth-adjusted.tsv -o module.xlsx
  $xlsx_py -i pathway-KO-depth.tsv,pathway-KO-depth-adjusted.tsv,pathway-KO-R.depth.tsv,pathway-KO-R.depth-adjusted.tsv,pathway-depth.tsv,pathway-depth-adjusted.tsv,pathway-R.depth.tsv,pathway-R.depth-adjusted.tsv -o pathway.xlsx
  $xlsx_py -i transporter-KO-depth.tsv,transporter-KO-depth-adjusted.tsv,transporter-KO-R.depth.tsv,transporter-KO-R.depth-adjusted.tsv -o transporter.xlsx
  $xlsx_py -i EC-KO-depth.tsv,EC-KO-depth-adjusted.tsv,EC-KO-R.depth.tsv,EC-KO-R.depth-adjusted.tsv -o EC.xlsx
  $xlsx_py -i AMR-KO-depth.tsv,AMR-KO-depth-adjusted.tsv,AMR-KO-R.depth.tsv,AMR-KO-R.depth-adjusted.tsv -o AMR.xlsx
fi

if [ -s assembly-num_scaffolds.txt ]
then
  $xlsx_py -i assembly-num_scaffolds.txt,assembly-len_scaffolds.txt,assembly-num_ORFs.txt -o assembly-summary.xlsx
  $xlsx_py -i assembly-genome-scaffolds.tsv,assembly-genome-length.tsv,assembly-genome-coverage.tsv,assembly-genome-depth.tsv,assembly-genome-ORFs.tsv -o assembly-genome.xlsx
  $xlsx_py -i assembly-species-scaffolds.tsv,assembly-species-length.tsv,assembly-species-coverage.tsv,assembly-species-depth.tsv,assembly-species-ORFs.tsv -o assembly-species.xlsx
fi

if [ -s taxon-superkingdom.tsv ]
then
  $xlsx_py -i taxon-depth-phylum.tsv,taxon-depth-class.tsv,taxon-depth-order.tsv,taxon-depth-family.tsv,taxon-depth-genus.tsv,taxon-depth-species.tsv,taxon-depth-toprank.tsv -o taxon-depth.xlsx
  $xlsx_py -i taxon-reads-phylum.tsv,taxon-reads-class.tsv,taxon-reads-order.tsv,taxon-reads-family.tsv,taxon-reads-genus.tsv,taxon-reads-species.tsv,taxon-reads-toprank.tsv -o taxon-reads.xlsx
  $xlsx_py -i qc.tsv,taxon-superkingdom-whost.tsv,taxon-superkingdom.tsv,taxon-phylum.tsv,taxon-class.tsv,taxon-order.tsv,taxon-family.tsv,taxon-genus.tsv,taxon-species.tsv,taxon-toprank.tsv -o taxon.xlsx
fi

