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
  'CMD_opts'         : ['100'],
  'execution'        : 'qsub_1',               # where to execute
  'cores_per_cmd'    : 4,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''
java -jar $ENV.NGS_root/apps/Trimmomatic/trimmomatic-0.32.jar PE -threads 4 -phred33 $DATA.0 $DATA.1 $SELF/R1.fq $SELF/R1-s.fq $SELF/R2.fq $SELF/R2-s.fq \\
    SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:$CMDOPTS.0 MAXINFO:80:0.5 1>$SELF/qc.stdout 2>$SELF/qc.stderr

perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R1.fq > $SELF/R1.fa &
perl -e '$i=0; while(<>){ if (/^@/) {$i++;  print ">Sample|$SAMPLE|$i ", substr($_,1); $a=<>; print $a; $a=<>; $a=<>;}}' < $SELF/R2.fq > $SELF/R2.fa &

wait
rm -f $SELF/R1.fq $SELF/R2.fq $SELF/R1-s.fq $SELF/R2-s.fq
'''
}


## remove reads from host (e.g. human for human microbiome)
## three methods:
##   bwa:      bwa based mapping, human genome formatted with bwa need to be at $ENV.NGS_root/refs/host
##   bowtie2:  bowtie2 based mapping, human genome formatted with bowtie2 need to be at $ENV.NGS_root/refs/host
##   skip:     do not run, for non-host related samples
NGS_batch_jobs['remove-host'] = {
  'injobs'         : ['qc'],          # start with high quality reads
  'CMD_opts'       : ['bwa'],         # can be bwa, bowtie2 or skip (do nothing for non-host related sample)
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

if [ "$CMDOPTS.0" = "bwa" ]
then
  $ENV.NGS_root/apps/bin/bwa mem -t 16 -a $ENV.NGS_root/refs/host $INJOBS.0/R1.fa | $ENV.NGS_root/NGS-tools/NGS-sam-raw-reduce-to-tophits.pl > $SELF/host-R1-top.sam &
  $ENV.NGS_root/apps/bin/bwa mem -t 16 -a $ENV.NGS_root/refs/host $INJOBS.0/R2.fa | $ENV.NGS_root/NGS-tools/NGS-sam-raw-reduce-to-tophits.pl > $SELF/host-R2-top.sam &
  wait
  cat $SELF/host-R1-top.sam | $ENV.NGS_root/apps/bin/samtools view -S - -F 0x004 | cut -f 1 | uniq > $SELF/host-hit-R1.ids
  cat $SELF/host-R2-top.sam | $ENV.NGS_root/apps/bin/samtools view -S - -F 0x004 | cut -f 1 | uniq > $SELF/host-hit-R2.ids
  cat $SELF/host-hit-R1.ids $SELF/host-hit-R2.ids | sort | uniq > $SELF/host-hit.ids
  rm -f $SELF/host-hit-R1.ids $SELF/host-hit-R2.ids
  $ENV.NGS_root/NGS-tools/NGS-fasta-fetch-exclude-ids.pl -i $SELF/host-hit.ids -s  $INJOBS.0/R1.fa -o $SELF/non-host-R1.fa
  $ENV.NGS_root/NGS-tools/NGS-fasta-fetch-exclude-ids.pl -i $SELF/host-hit.ids -s  $INJOBS.0/R2.fa -o $SELF/non-host-R2.fa
  
elif [ "$CMDOPTS.0" = "bowtie2" ]
then
  $ENV.NGS_root/apps/bin/bowtie -f -k 1 -v 2 -p 16 $ENV.NGS_root/refs/host $INJOBS.0/R1.fa $SELF/host-hit-1 &
  $ENV.NGS_root/apps/bin/bowtie -f -k 1 -v 2 -p 16 $ENV.NGS_root/refs/host $INJOBS.0/R2.fa $SELF/host-hit-2 &
  wait
  cut -f 1 $SELF/host-hit-1 > $SELF/host-hit-R1.ids
  cut -f 1 $SELF/host-hit-2 > $SELF/host-hit-R2.ids
  rm -f $SELF/host-hit-R1.ids $SELF/host-hit-R2.ids
  $ENV.NGS_root/NGS-tools/NGS-fasta-fetch-exclude-ids.pl -i $SELF/host-hit.ids -s  $INJOBS.0/R1.fa -o $SELF/non-host-R1.fa
  $ENV.NGS_root/NGS-tools/NGS-fasta-fetch-exclude-ids.pl -i $SELF/host-hit.ids -s  $INJOBS.0/R2.fa -o $SELF/non-host-R2.fa

elif [ "$CMDOPTS.0" = "skip" ]
then
  #### do nothing, simply link 
  ln -s ../$INJOBS.0/R1.fa $SELF/non-host-R1.fa
  ln -s ../$INJOBS.0/R2.fa $SELF/non-host-R2.fa

else
  echo "not defined filter-host method"
  exit 1  
fi
'''
}

NGS_batch_jobs['reads-mapping'] = {
  'injobs'         : ['remove-host'],
  'CMD_opts'       : ['75'],          # significant score cutoff
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

$ENV.NGS_root/apps/bin/bwa mem -t 16 -T $CMDOPTS.0 -a $ENV.NGS_root/refs/ref-genomes $INJOBS.0/non-host-R1.fa | $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/R1-raw.bam
$ENV.NGS_root/apps/bin/bwa mem -t 16 -T $CMDOPTS.0 -a $ENV.NGS_root/refs/ref-genomes $INJOBS.0/non-host-R2.fa | $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/R2-raw.bam
$ENV.NGS_root/NGS-tools/NGS-sam-raw-reduce-to-tophits-2SE.pl -i $SELF/R1-raw.bam -j $SELF/R2-raw.bam -T $CMDOPTS.0 \\
  -l $SELF/bam-filter-log -p $ENV.NGS_root/apps/bin/samtools -b 1 | $ENV.NGS_root/apps/bin/samtools view -b -S - > $SELF/R12-top.bam

$ENV.NGS_root/NGS-tools/NGS-sam-genome-cov.pl -i $SELF/R12-top.bam -o $SELF/R12--genome-cov -b 1 -t $ENV.NGS_root/apps/bin/samtools
$ENV.NGS_root/NGS-tools/NGS-sam-genome-cov-filter-sam.pl -i $SELF/R12-top.bam -j $SELF/R12--genome-cov -c 0.1 -b 1 -t $ENV.NGS_root/apps/bin/samtools -o $SELF/R12-topf.sam 
mv $SELF/R12-topf.sam $SELF/R12-top.sam
$NGS_bin_dir/samtools view -S $SELF/R12-T60-top.sam         -F 0x004 | cut -f 1 | uniq > $SELF/mapped-genome.ids

'''
}


NGS_batch_jobs['assembly'] = {
  'injobs'         : ['remove-host'],
  'CMD_opts'       : ['spade', '250'],       # can be idba-ud or spade , ## 250 confit length cutoff
  'non_zero_files' : ['assembly/scaffold.fa'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
if [ "$CMDOPTS.0" = "idba-ud" ]
then
  $ENV.NGS_root/NGS-tools/PE-2file-to-1file.pl -i $INJOBS.0/non-host-R1.fa,$INJOBS.0/non-host-R2.fa -r 0 > $SELF/input.fa
  $ENV.NGS_root/apps/bin/idba_ud -r $SELF/input.fa -o $SELF/assembly --mink=50 --maxk=80 --step=10 --num_threads=16 --min_contig=$CMDOPTS.1 
  rm -f $SELF/input.fa $SELF/assembly/kmer $SELF/assembly/local* $SELF/assembly/contig* $SELF/assembly/graph*
  
elif [ "$CMDOPTS.0" = "spade" ]
then
  python $ENV.NGS_root/apps/SPAdes/bin/spades.py -1 $INJOBS.0/non-host-R1.fa -2 $INJOBS.0/non-host-R2.fa --meta --only-assembler -o $SELF/assembly -t 16
  mv $SELF/assembly/scaffolds.fasta $SELF/assembly/scaffold.fa
  rm -rf $SELF/assembly/K*
  rm -rf $SELF/assembly/tmp

else
  echo "not defined assembly method"
  exit 1  
fi


$ENV.NGS_root/NGS-tools/fasta_filter_short_seq.pl -i $SELF/assembly/scaffold.fa -c $CMDOPTS.1 -o - -a len | \\
  $ENV.NGS_root/NGS-tools/fasta_rename.pl -i - -s "$SAMPLE|scaffold|" -b 0 -o $SELF/assembly/scaffold-new.fa
mv -f $SELF/assembly/scaffold-new.fa $SELF/assembly/scaffold.fa

## $NGS_bin_dir/bwa index -a bwtsw -p $SELF/scaffold $SELF/scaffold.fa
## insert coverage 
'''
}


NGS_batch_jobs['ORF-prediction'] = {
  'injobs'         : ['assembly'],
  'CMD_opts'       : ['metagene', '20'],    # can be metagene or prodigal 
  'non_zero_files' : ['ORF.faa'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 1,               # number of threads used by command below
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

$ENV.NGS_root/NGS-tools/fasta_filter_short_seq.pl -i $SELF/ORF.faa -c $CMDOPTS.1 -a len -o $SELF/ORF-new.faa
mv  $SELF/ORF-new.faa $SELF/ORF.faa
'''
}

NGS_batch_jobs['cd-hit-kegg'] = {
  'injobs'         : ['ORF-prediction'],
  'CMD_opts'       : ['kegg/kegg'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
$ENV.NGS_root/apps/cd-hit/cd-hit-2d -i $ENV.NGS_root/refs/$CMDOPTS.0 -i2 $INJOBS.0/ORF.faa -o $SELF/out \\
  -c 0.75 -n 5 -d 0 -g 1 -G 0 -aS 0.9 -T 16 -M 32000 > $SELF/out.log
$ENV.NGS_root/apps/cd-hit/cd-hit-clstr_2_blm8.pl < $SELF/out.clstr > $SELF/out.bl

mkdir $SELF/orf-split
$ENV.NGS_root/apps/cd-hit/cd-hit-div.pl $SELF/out $SELF/orf-split/split 64
'''
}


#NGS_batch_jobs['orf-split'] = {
#  'injobs'         : ['ORF-prediction'],
#  'execution'      : 'qsub_1',        # where to execute
#  'cores_per_cmd'  : 16,              # number of threads used by command below
#  'no_parallel'    : 1,               # number of total jobs to run using command below
#  'command'        : '''
#mkdir $SELF/orf-split
#$ENV.NGS_root/apps/cd-hit/cd-hit-div.pl $INJOBS.0/ORF.faa               $SELF/orf-split/split        256
#'''
#}

NGS_batch_jobs['blast-kegg'] = {
  'injobs'         : ['cd-hit-kegg'],
  'CMD_opts'       : ['kegg/kegg'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 16,              # number of threads used by command below
  'no_parallel'    : 2,               # number of total jobs to run using command below
  'command'        : '''

for i in `seq 1 4`
  do $ENV.NGS_root/NGS-tools/ann_batch_run_dir.pl --INDIR1=$INJOBS.0/orf-split --OUTDIR1=$SELF/blast --CPU=$SELF/WF.cpu $ENV.NGS_root/apps/blast+/bin/blastp  -query {INDIR1} -out {OUTDIR1} \\
  -db $ENV.NGS_root/refs/$CMDOPTS.0 -evalue 0.001 -num_threads 4 -num_alignments 5 -outfmt 6 -seg yes &
done
wait

'''
}

NGS_batch_jobs['blast-kegg-parse'] = {
  'injobs'         : ['blast-kegg','cd-hit-kegg','ORF-prediction'],
  'CMD_opts'       : ['kegg/kegg_all.faa'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 2,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''
ln -s ../../$INJOBS.1/out.bl $INJOBS.0/blast/out.bl

#$ENV.NGS_root/NGS-tools/ann_parse_blm8.pl     -i $INJOBS.0/blast -o $SELF/protein-ann.txt  -d $ENV.NGS_root/refs/$CMDOPTS.0
#$ENV.NGS_root/NGS-tools/ann_parse_blm8-raw.pl -i $INJOBS.0/blast -o $SELF/protein-full.txt -d $ENV.NGS_root/refs/$CMDOPTS.0

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func.pl -i $INJOBS.0/blast -r $ENV.NGS_root/refs/kegg/kegg.clstr.ann \\
  -a $INJOBS.2/ORF.faa -o $SELF/ORF -t $ENV.NGS_root/refs/kegg/kegg_taxon.txt

'''
}

NGS_batch_jobs['blastx-viral'] = {
  'injobs'         : ['assembly'],
  'CMD_opts'       : ['ref-genomes/ref_viral_prot_95'],
  'execution'      : 'qsub_1',        # where to execute
  'cores_per_cmd'  : 8,              # number of threads used by command below
  'no_parallel'    : 1,               # number of total jobs to run using command below
  'command'        : '''

$ENV.NGS_root/apps/blast+/bin/blastx -query  $INJOBS.0/assembly/scaffold.fa -out $SELF/viral-bl \\
  -db $ENV.NGS_root/refs/ref-genomes/ref_viral_prot_95 -evalue 0.001 -num_threads 8 -num_alignments 5 -outfmt 6 -seg yes

$ENV.NGS_root/NGS-tools/ann_blastx_orf_call_format.pl -i $SELF/viral-bl -o $SELF/viral-blnew
perl -e 'while(<>){ chop($_); $_ =~ s/\s.+$//; print ">$_\\nXXXXXXXXXX\\n"; }' < $SELF/viral-blnew > $SELF/blastx-ORF

$ENV.NGS_root/NGS-tools/ann_ORF_taxon_func.pl -i $SELF/viral-blnew -r $ENV.NGS_root/refs/ref-genomes/ref_viral_prot_95.clstr.ann \\
  -a $SELF/blastx-ORF -o $SELF/ORF -t $ENV.NGS_root/refs/ref-genomes/ref_viral_prot_taxon.txt
'''
}

