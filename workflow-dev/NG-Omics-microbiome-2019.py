#!/usr/bin/python
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

queue_system = 'SGE'

########## local variables etc. Please edit
ENV={
  'NGS_root' : '/home/oasis/gordon-data/NGS-ann-project-new',
}

########## computation resources for execution of jobs
NGS_executions = {}
NGS_executions['qsub_1'] = {
  'type'                : 'qsub-pe',
  'pe_para'             : '-pe orte', #### '-pe orte' work with orte environment with $pe_slot allocation rule, '-pe threaded' work with multiple threading env
  'qsub_exe'            : 'qsub',
  'cores_per_node'      : 32,
  'number_nodes'        : 64,
  'user'                : 'weizhong', #### I will use command such as qstat -u weizhong to query submitted jobs
  'command'             : 'qsub',
  'command_name_opt'    : '-N',
  'command_err_opt'     : '-e',
  'command_out_opt'     : '-o',
  'template'            : '''#!/bin/bash
#$ -q RNA.q
#$ -v PATH
#$ -V

'''
}

NGS_executions['sh_1'] = {
  'type'                : 'sh',
  'cores_per_node'      : 8,
  'number_nodes'        : 1,
  'template'            : '''#!/bin/bash

'''
}


NGS_batch_jobs = {}
NGS_batch_jobs['qc'] = {
  'CMD_opts'         : ['100','NexteraPE'],
  'non_zero_files' : ['R1.fa.gz','R2.fa.gz','qc.txt'],
  'execution'        : 'qsub_1',               # where to execute
  'cores_per_cmd'    : 4,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

#### Trim adapter options
if [ "$CMDOPTS.1" = "NexteraPE" ] 
then
  java -jar $ENV.NGS_root/apps/Trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DATA.0 $DATA.1 $SELF/R1.fq $SELF/R1-s.fq $SELF/R2.fq $SELF/R2-s.fq \\
      ILLUMINACLIP:$ENV.NGS_root/apps/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 \\
      SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:$CMDOPTS.0 MAXINFO:80:0.5 1>$SELF/qc.stdout 2>$SELF/qc.stderr
else
  java -jar $ENV.NGS_root/apps/Trimmomatic/trimmomatic-0.36.jar PE -threads 4 -phred33 $DATA.0 $DATA.1 $SELF/R1.fq $SELF/R1-s.fq $SELF/R2.fq $SELF/R2-s.fq \\
      SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:$CMDOPTS.0 MAXINFO:80:0.5 1>$SELF/qc.stdout 2>$SELF/qc.stderr
fi

perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R1.fq > $SELF/R1.fa &
perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R2.fq > $SELF/R2.fa &
wait

gzip $SELF/R1.fa &
gzip $SELF/R2.fa &
wait

rm -f $SELF/R1.fq $SELF/R2.fq $SELF/R1-s.fq $SELF/R2-s.fq

#### qc.txt
NUM_reads_total=$(grep "Input Read Pairs"  $SELF/qc.stderr | cut -f 4 -d " ")
NUM_reads_pass=$(grep "Input Read Pairs"  $SELF/qc.stderr | cut -f 7 -d " ")
echo -e "#Reads\\tNumber" > $SELF/qc.txt
echo -e "Total_reads\\t$NUM_reads_total" >> $SELF/qc.txt
echo -e "QC_reads\\t$NUM_reads_pass" >> $SELF/qc.txt

'''
}


## remove reads from host (e.g. human for human microbiome)
## three methods:
##   bwa:      bwa based mapping, human genome formatted with bwa need to be at $ENV.NGS_root/refs/host
##   bowtie2:  bowtie2 based mapping, human genome formatted with bowtie2 need to be at $ENV.NGS_root/refs/host
##   skip:     do not run, for non-host related samples
NGS_batch_jobs['reads-filtering'] = {
  'injobs'         : ['qc'],          # start with high quality reads
  'CMD_opts'       : ['host/GRCh38.fa','NULL'],  # 1st option db1 for filtering against host, metaG and metaT, NONE to skip
#  'CMD_opts'       : ['host/GRCh38.fa','total_RNA/total_RNA'],  # 1st option db1 for filtering against host, metaG and metaT, NONE to skip
                                                                # 2nd option db2 for filtering against rRNA etc, metaT, NONE to skip
  'non_zero_files' : ['filtered-R1.fa.gz','filtered-R2.fa.gz'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

#### db1 for host filtering
if [ ! "$CMDOPTS.0" = "NONE" ]; then
#### -F 1 for sam-filter-top-pair-or-single.pl means only accept the reads that are mapped in proper pairs
  $ENV.NGS_root/apps/bin/bwa mem -t 16 -T 60 $ENV.NGS_root/refs/$CMDOPTS.0 \\
    $INJOBS.0/R1.fa.gz $INJOBS.0/R2.fa.gz | $ENV.NGS_root/NGS-tools/sam-filter-top-pair-or-single.pl -T 60 -F 0 -O $SELF/host-hit.ids | \\
    $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/host.top.bam
else
  touch $SELF/host-hit.ids
fi

#### db2 for rRNA filtering
if [ ! "$CMDOPTS.1" = "NONE" ]; then
  $ENV.NGS_root/apps/bin/bwa mem -t 16 -T 60  $ENV.NGS_root/refs/$CMDOPTS.1 \\
    $INJOBS.0/R1.fa.gz $INJOBS.0/R2.fa.gz | $ENV.NGS_root/NGS-tools/sam-filter-top-pair-or-single.pl -T 60 -F 0 -O $SELF/rRNA-hit.ids | \\
    $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/rRNA.top.bam
else
  touch $SELF/rRNA-hit.ids
fi

cat $SELF/host-hit.ids $SELF/rRNA-hit.ids | sort | uniq > $SELF/filter-hit.ids

if [ -s "$SELF/filter-hit.ids" ]; then
  $ENV.NGS_root/NGS-tools/fasta_fetch_exclude_ids.pl -i $SELF/filter-hit.ids -s  $INJOBS.0/R1.fa.gz -o $SELF/filtered-R1.fa &
  $ENV.NGS_root/NGS-tools/fasta_fetch_exclude_ids.pl -i $SELF/filter-hit.ids -s  $INJOBS.0/R2.fa.gz -o $SELF/filtered-R2.fa &
  wait

  gzip $SELF/filtered-R1.fa &
  gzip $SELF/filtered-R2.fa &
  wait
else
  #### do nothing, simply link 
  ln -s ../$INJOBS.0/R1.fa.gz $SELF/filtered-R1.fa.gz
  ln -s ../$INJOBS.0/R2.fa.gz $SELF/filtered-R2.fa.gz
fi

#### filter.txt
NUM_HOST=$(grep -c "." $SELF/host-hit.ids)
NUM_RRNA=$(grep -c "." $SELF/rRNA-hit.ids)
NUM_HOST_RNA=$(grep -c "." $SELF/filter-hit.ids)
echo -e "#Filter\\tdb\\tNumber"                         > $SELF/filter.txt
echo -e "Host\\tHost\\t$NUM_HOST"                      >> $SELF/filter.txt
echo -e "rRNA_tRNA\\trRNA_tRNA\\t$NUM_TRNA_RRNA"       >> $SELF/filter.txt
echo -e "Both_Host_RNA\\tBoth_Host_RNA\\tNUM_HOST_RNA" >> $SELF/filter.txt

'''
}

#### mitochondron analysis
#### ~100 mitochondria in a mammalian cell, and each mitochondrion has 2-10 copies of mtDNA
NGS_batch_jobs['mito-ana'] = {
  'injobs'         : ['qc'],          # start with high quality reads
  'CMD_opts'       : ['mito/HS_mito'],         # can be bwa, bowtie2 or skip (do nothing for non-host related sample)
  'non_zero_files' : ['mito.vcf'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

$ENV.NGS_root/apps/bin/bwa mem -t 16 -T 60 $ENV.NGS_root/refs/$CMDOPTS.0 \\
  $INJOBS.0/R1.fa.gz $INJOBS.0/R2.fa.gz | $ENV.NGS_root/NGS-tools/sam-filter-top-pair-or-single.pl -T 60 | \\
  $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/mito.top.bam

$ENV.NGS_root/apps/bin/samtools view $SELF/mito.top.bam -F 0x004 | cut -f 1 > $SELF/mito-hit.ids
$ENV.NGS_root/apps/bin/samtools view $SELF/mito.top.bam | \\
  $ENV.NGS_root/NGS-tools/sam-HS-depth.pl -i - -a $ENV.NGS_root/refs/$CMDOPTS.0.ann -o $SELF/mito-depth

$ENV.NGS_root/apps/bin/samtools sort $SELF/mito.top.bam -o $SELF/mito.sorted.bam
$ENV.NGS_root/apps/bin/samtools mpileup -o $SELF/mito.mpileup -f /local/ifs2_projdata/8460/projects/WOUND/wli/refs/mito/HS_mito    $SELF/mito.sorted.bam
$ENV.NGS_root/apps/bin/samtools mpileup -o $SELF/mito.raw.vcf -f /local/ifs2_projdata/8460/projects/WOUND/wli/refs/mito/HS_mito -v $SELF/mito.sorted.bam

java -jar $ENV.NGS_root/apps/VarScan-2.4.2/VarScan.v2.4.2.jar  pileup2snp $SELF/mito.mpileup > $SELF/mito.vcf

#$ENV.NGS_root/NGS-tools/fasta_fetch_by_ids.pl -i $SELF/mito-hit.ids -s $INJOBS.0/R1.fa.gz -o $SELF/mito-R1.fa
#$ENV.NGS_root/NGS-tools/fasta_fetch_by_ids.pl -i $SELF/mito-hit.ids -s $INJOBS.0/R2.fa.gz -o $SELF/mito-R2.fa

'''
}

NGS_batch_jobs['centrifuge-mapping'] = {
  'injobs'         : ['reads-filtering'],
  'CMD_opts'       : ['ref-genomes-2018-1109/centrifuge/ref_full','ref-genomes-2018-1109/centrifuge/map',
                      'ref-genomes-2018-1109/refseq_genome_taxon.tsv', 'ref-genomes-2018-1109/centrifuge/ref_full.ann'], ## 2018 version
  'non_zero_files' : ['centrifuge-out.gz','centrifuge-taxon','centrifuge-report'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

$ENV.NGS_root/apps/centrifuge/centrifuge -1 $INJOBS.0/filtered-R1.fa.gz -2 $INJOBS.0/filtered-R2.fa.gz \\
  -x $ENV.NGS_root/refs/$CMDOPTS.0 \\
  -S $SELF/centrifuge-out --report-file $SELF/centrifuge-taxon -p 16 -f --out-fmt tab 
$ENV.NGS_root/apps/centrifuge/centrifuge-kreport -x $ENV.NGS_root/refs/$CMDOPTS.0 $SELF/centrifuge-out > $SELF/centrifuge-report

grep -vP "0.0$" $SELF/centrifuge-taxon > $SELF/centrifuge-taxon.1
(head -n 1 $SELF/centrifuge-taxon.1 && tail -n +2 $SELF/centrifuge-taxon.1 | sort -grt $'\\t' -k 7,7 ) > $SELF/centrifuge-taxon.2

# count total number of input reads
NUM_reads=$(zcat $INJOBS.0/filtered-R1.fa.gz | grep -c "^>")
$ENV.NGS_root/NGS-tools/centrifuge-taxon.pl -i $SELF/centrifuge-out -j $SELF/centrifuge-report \\
  -t $ENV.NGS_root/refs/$CMDOPTS.2 -a $ENV.NGS_root/refs/$CMDOPTS.3 \\
  -o $SELF/taxon -c 1e-7 -N $NUM_reads -l 60

gzip -f $SELF/centrifuge-out

cp -p $SELF/taxon.superkingdom.txt $SELF/taxon.superkingdom-whost.txt
grep "^Host" $INJOBS.0/filter.txt >> $SELF/taxon.superkingdom-whost.txt
grep "^rRNA" $INJOBS.0/filter.txt >> $SELF/taxon.superkingdom-whost.txt

'''
}

NGS_batch_jobs['assembly'] = {
  'injobs'         : ['reads-filtering'],
  'CMD_opts'       : ['spade', '250'],       # can be idba-ud or spade , ## 250 confit length cutoff
  'non_zero_files' : ['assembly/scaffold.fa'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
if [ "$CMDOPTS.0" = "idba-ud" ]
then
  $ENV.NGS_root/NGS-tools/PE-2file-to-1file.pl -i $INJOBS.0/filtered-R1.fa.gz,$INJOBS.0/filtered-R2.fa.gz -r 0 > $SELF/input.fa
  $ENV.NGS_root/apps/bin/idba_ud -r $SELF/input.fa -o $SELF/assembly --mink=50 --maxk=80 --step=10 --num_threads=16 --min_contig=$CMDOPTS.1 
  rm -f $SELF/input.fa $SELF/assembly/kmer $SELF/assembly/local* $SELF/assembly/contig* $SELF/assembly/graph*
  
elif [ "$CMDOPTS.0" = "spade" ]
then
  python $ENV.NGS_root/apps/SPAdes/bin/spades.py -1 $INJOBS.0/filtered-R1.fa.gz -2 $INJOBS.0/filtered-R2.fa.gz --meta --only-assembler -o $SELF/assembly -t 16
  mv $SELF/assembly/scaffolds.fasta $SELF/assembly/scaffold.fa
  rm -rf $SELF/assembly/K*
  rm -rf $SELF/assembly/tmp

else
  echo "not defined assembly method"
  exit 1  
fi


## ensure to have >SAMPLE_name|scaffold|n header format
$ENV.NGS_root/NGS-tools/fasta_filter_short_seq.pl -i $SELF/assembly/scaffold.fa -c $CMDOPTS.1 -o - -a len | \\
  $ENV.NGS_root/NGS-tools/fasta_rename.pl -i - -s "$SAMPLE|scaffold|" -b 0 -o $SELF/assembly/scaffold-new.fa
mv -f $SELF/assembly/scaffold-new.fa $SELF/assembly/scaffold.fa

## depth of coverage

perl -e 'while(<>){ if ($_ =~ /^>(\S+)/) { $id=$1;  if ($_ =~ /_cov_([\d\.]+)/) { print "$id\\t$1\\n";} } }' < $SELF/assembly/scaffold.fa > $SELF/assembly/scaffold-cov

## $NGS_bin_dir/bwa index -a bwtsw -p $SELF/scaffold $SELF/scaffold.fa
## insert coverage 
'''
}

NGS_batch_jobs['ORF-prediction'] = {
  'injobs'         : ['assembly'],
  'CMD_opts'       : ['prodigal', '20'],    # can be metagene or prodigal 
  'non_zero_files' : ['ORF.faa'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 4,               # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
if [ "$CMDOPTS.0" = "metagene" ]
then
  $ENV.NGS_root/apps/metagene/metagene_run.pl $INJOBS.0/assembly/scaffold.fa $SELF/ORF.faa

elif [ "$CMDOPTS.0" = "prodigal" ]
then
  $ENV.NGS_root/apps/Prodigal/prodigal -i $INJOBS.0/assembly/scaffold.fa -p meta -o $SELF/ORF.gff -f gff -a $SELF/ORF.faa -d $SELF/ORF.fna

else
  echo "undefined ORF-prediction method"
  exit 1  
fi

$ENV.NGS_root/NGS-tools/fasta_filter_short_seq.pl -i $SELF/ORF.faa -c $CMDOPTS.1 -o $SELF/ORF-new.faa
mv  $SELF/ORF-new.faa $SELF/ORF.faa

$ENV.NGS_root/NGS-tools/assembly-cov-pass-to-orf.pl -i $SELF/ORF.faa -d $INJOBS.0/assembly/scaffold-cov -o $SELF/ORF-cov
'''
}


NGS_batch_jobs['minimap-binning'] = {
  'injobs'         : ['assembly','centrifuge-mapping','ORF-prediction'],
  'CMD_opts'       : ['0.01','ref-genomes-2018-1109/bwa/ref_full.fna','ref-genomes-2018-1109/bwa/ref_full.fna.header',
                      'ref-genomes-2018-1109/refseq_genome_taxon.tsv'],
  'non_zero_files' : ['assembly-bin'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

#### print taxids with depth > cutoff
grep -v "^#" $INJOBS.1/taxon.toprank.txt | gawk -F '\\t' '($12+0.0) > $CMDOPTS.0 {print $1}' > $SELF/high_cov.taxids
$ENV.NGS_root/NGS-tools/fasta_idx_fetch_by_ids.pl -i $SELF/high_cov.taxids -s $ENV.NGS_root/refs/$CMDOPTS.1 \\
  -a taxid -h $ENV.NGS_root/refs/$CMDOPTS.2 -o $SELF/high_cov_ref.fna

$ENV.NGS_root/apps/minimap2/minimap2 $SELF/high_cov_ref.fna $INJOBS.0/assembly/scaffold.fa -N 20 -t 8 > $SELF/hits.paf
rm -f $SELF/high_cov_ref.fna

$ENV.NGS_root/NGS-tools/JCVI/minimap-binning-taxon.pl -i $SELF/hits.paf -s $INJOBS.0/assembly/scaffold.fa \\
  -t $ENV.NGS_root/refs/$CMDOPTS.3 -o $SELF/assembly-bin
'''
}


NGS_batch_jobs['cd-hit-kegg'] = {
  'injobs'         : ['ORF-prediction'],
  'CMD_opts'       : ['kegg/keggf'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
$ENV.NGS_root/apps/cd-hit/cd-hit-2d -i $ENV.NGS_root/refs/$CMDOPTS.0 -i2 $INJOBS.0/ORF.faa -o $SELF/out \\
  -c 0.75 -n 5 -d 0 -g 1 -G 0 -aS 0.9 -A 60 -aL 0.25 -T 16 -M 32000 > $SELF/out.log
$ENV.NGS_root/apps/cd-hit/clstr_select.pl 2 99999999 < $SELF/out.clstr > $SELF/out.clstr.1
mv -f $SELF/out.clstr.1 $SELF/out.clstr
$ENV.NGS_root/apps/cd-hit/cd-hit-clstr_2_blm8.pl < $SELF/out.clstr > $SELF/out.bl

mkdir $SELF/orf-split
$ENV.NGS_root/apps/cd-hit/cd-hit-div.pl $SELF/out $SELF/orf-split/split 64
'''
}

NGS_batch_jobs['orf-split'] = {
  'injobs'         : ['ORF-prediction'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
mkdir $SELF/orf-split
$ENV.NGS_root/apps/cd-hit/cd-hit-div.pl $INJOBS.0/ORF.faa               $SELF/orf-split/split        256
'''
}

NGS_batch_jobs['pfam'] = {
  'injobs'         : ['orf-split'],
  'CMD_opts'       : ['pfam/Pfam-A.hmm'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 4,              # number of threads used by command below
  'no_parallel'    : 4,               # number of total jobs to run using command below
  'command'        : '''
for i in `seq 1 4`
  do $ENV.NGS_root/NGS-tools/ann_batch_run_dir.pl --INDIR1=$INJOBS.0/orf-split --OUTDIR1=$SELF/pfam --OUTDIR2=$SELF/pfam.2 --OUTDIR3=$SELF/pfam.3 \\
    --CPU=$SELF/WF.cpu $ENV.NGS_root/apps/hmmer/binaries/hmmscan -E 0.001 -o {OUTDIR1} --notextw --noali --cpu 1 --tblout {OUTDIR2} \\
    --domtblout {OUTDIR3} $ENV.NGS_root/refs/$CMDOPTS.0 {INDIR1} &
done
wait

'''
}

NGS_batch_jobs['pfam-parse'] = {
  'injobs'         : ['pfam', 'ORF-prediction'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 2,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
$ENV.NGS_root/NGS-tools/ann_parse_hmm.pl -i $INJOBS.0/pfam.3 -o $SELF/pfam-ann.txt -a $INJOBS.1/ORF-cov
'''
}

NGS_batch_jobs['blast-kegg'] = {
  'injobs'         : ['cd-hit-kegg'],
  'CMD_opts'       : ['kegg/keggf','Run'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 2,               # number of total jobs to run using command below
  'command'        : '''

#### skip this - only use cd-hit-kegg's result
if [ "$CMDOPTS.1" = "Skip" ] 
then
  mkdir $SELF/blast
  echo "Skip" >> $SELF/skip.txt
else
  for i in `seq 1 4`
    do $ENV.NGS_root/NGS-tools/ann_batch_run_dir.pl --INDIR1=$INJOBS.0/orf-split --OUTDIR1=$SELF/blast --CPU=$SELF/WF.cpu $ENV.NGS_root/apps/blast+/bin/blastp  -query {INDIR1} -out {OUTDIR1} \\
    -db $ENV.NGS_root/refs/$CMDOPTS.0 -evalue 1e-6 -num_threads 4 -num_alignments 5 -outfmt 6 -seg yes &
  done
  wait
fi

'''
}

NGS_batch_jobs['blast-kegg-parse'] = {
  'injobs'         : ['blast-kegg','cd-hit-kegg','ORF-prediction','minimap-binning','assembly'],
  'CMD_opts'       : ['kegg/keggf_taxon.txt', 'ref-genomes-2018-1109/refseq_genome_taxon.tsv', 'kegg/keggf.clstr.ann', 'kegg/ko00002-M00178.keg'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 2,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
ln -s ../../$INJOBS.1/out.bl $INJOBS.0/blast/out.bl

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_prebin.pl -i $INJOBS.0/blast -r $ENV.NGS_root/refs/$CMDOPTS.2 \\
  -a $INJOBS.2/ORF.faa -o $SELF/ORF -t $ENV.NGS_root/refs/$CMDOPTS.0 -s $INJOBS.3/assembly-bin \\
  -x $ENV.NGS_root/refs/$CMDOPTS.1 -p $SELF/scaffold-ann.txt -d $INJOBS.4/assembly/scaffold-cov

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00001-biome.keg -o $SELF/ORF-ann-kegg-pathway     -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00002-edit.keg  -o $SELF/ORF-ann-kegg-module      -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01000.keg       -o $SELF/ORF-ann-kegg-EC          -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko02000.keg       -o $SELF/ORF-ann-kegg-transporter -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01504.keg       -o $SELF/ORF-ann-kegg-AMR         -r $ENV.NGS_root/refs/$CMDOPTS.3

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00001-biome.keg -o $SELF/sp-pathway     -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00002-edit.keg  -o $SELF/sp-module      -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01000.keg       -o $SELF/sp-EC          -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko02000.keg       -o $SELF/sp-transporter -r $ENV.NGS_root/refs/$CMDOPTS.3
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01504.keg       -o $SELF/sp-AMR         -r $ENV.NGS_root/refs/$CMDOPTS.3

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00001-biome.keg -o $SELF/ge-pathway     -r $ENV.NGS_root/refs/$CMDOPTS.3 -K 8
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko00002-edit.keg  -o $SELF/ge-module      -r $ENV.NGS_root/refs/$CMDOPTS.3 -K 8
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01000.keg       -o $SELF/ge-EC          -r $ENV.NGS_root/refs/$CMDOPTS.3 -K 8
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko02000.keg       -o $SELF/ge-transporter -r $ENV.NGS_root/refs/$CMDOPTS.3 -K 8
$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func_kegg_per_sp.pl -i $SELF/ORF-ann.txt -d $INJOBS.2/ORF-cov -k $ENV.NGS_root/refs/kegg/ko01504.keg       -o $SELF/ge-AMR         -r $ENV.NGS_root/refs/$CMDOPTS.3 -K 8

'''
}


NGS_batch_jobs['ann-summary'] = {
  'injobs'         : ['minimap-binning','blast-kegg-parse','assembly'],
  'CMD_opts'       : ['kegg/keggf_taxon.txt','ref-genomes-2018-1109/refseq_genome_taxon.tsv','ref-genomes-2018-1109/centrifuge/ref_full.ann'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 4,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
$ENV.NGS_root/NGS-tools/JCVI/assembly-orf-summary.pl -t $ENV.NGS_root/refs/$CMDOPTS.0 -x $ENV.NGS_root/refs/$CMDOPTS.1 \\
  -b $INJOBS.0/assembly-bin -i $INJOBS.1/scaffold-ann.txt -o $SELF/assembly-orf-summary.txt \\
  -a $ENV.NGS_root/refs/$CMDOPTS.2 -p $SELF/assembly-taxon -N 10 -d $INJOBS.2/assembly/scaffold-cov

for i in `ls -1 $SELF/assembly-taxon*sids`;
do
    $ENV.NGS_root/NGS-tools/fasta_fetch_by_ids.pl -i $i -s $INJOBS.2/assembly/scaffold.fa  -o $i.faa
    rm -f $i
done;

'''
}


