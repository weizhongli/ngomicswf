#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

## input BAM / SAM must already be filtered to remove non-top hits 

use Getopt::Std;
getopts("i:o:a:t:r:N:c:P:",\%opts);
die usage() unless ($opts{o} and $opts{a});


my $sam_file    = $opts{i}; ### input sam file, usually pipe in from samtools view input.bam
   $sam_file    = "-" unless ($sam_file);
my $bwa_db_ann  = $opts{a}; ### bwa db .ann file
my $output      = $opts{o};
my $read_length = $opts{r}; #### read length
   $read_length = 150*2 unless ($read_length);

my ($i, $j, $k, $ll, $cmd);

#### read genome length
# /local/ifs2_projdata/8460/projects/WOUND/wli/refs/host/GRCh38.fa.ann
# 3099750718 194 11
# 0 NC_000001.11 Homo sapiens chromosome 1, GRCh38.p7 Primary Assembly
# 0 248956422 168
# 0 NC_000002.12 Homo sapiens chromosome 2, GRCh38.p7 Primary Assembly

my %refid_2_len = ();
my %refid_2_des = ();
open(TMP, $bwa_db_ann) || die "can not open $bwa_db_ann";
$ll=<TMP>;
while($ll=<TMP>){
  chop($ll);
  my ($t1, $refid, $des) = split(/\s+/, $ll, 3);
  $refid_2_des{$refid} = $des;

  $ll=<TMP>;
  chop($ll);
  $refid_2_len{$refid} = (split(/\s+/, $ll))[1];
}
close(TMP);


my $fh;
if ($sam_file eq "-") { $fh = "STDIN";}
else {
  open(TMP, $sam_file) || die "can not open $sam_file";
  $fh = "TMP";
}

my %refid_reads_count = ();
my $last_id = "";
my %hit_refids = ();
my $num_mapped_reads = 0;
while($ll=<$fh>){
  if ($ll =~ /^\@/) { #### headers
    next;
  }
  else { #### alignment
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $id  = $lls[0];
    my $refid = $lls[2];    if ($refid eq "*") {  next; }

    my $flag_info = <<EOD;
        Chr     Flag    Description
        p       0x0001  the read is paired in sequencing
        P       0x0002  the read is mapped in a proper pair
        u       0x0004  the query sequence itself is unmapped
        U       0x0008  the mate is unmapped
        r       0x0010  strand of the query (1 for reverse)
        R       0x0020  strand of the mate
        1       0x0040  the read is the first read in a pair
        2       0x0080  the read is the second read in a pair
        s       0x0100  the alignment is not primary
        f       0x0200  QC failure
        d       0x0400  optical or PCR duplicate
EOD

    my $FLAG = $lls[1];
    next unless ($FLAG & 0x0001 );
    next unless ($FLAG & 0x0002 );

    if    (($id    ne $last_id) and $last_id ) {
      #### for this read pair, if both R1 and R2 hit a refid (t1) and only R1 or R2 hit another refid (t2), skip t2 
      my @t = sort {$b <=> $a} values %hit_refids;
      my $m = $t[0]; #### the refid with most 
      foreach $i (keys %hit_refids) {
        delete($hit_refids{$i}) unless ($m == $hit_refids{$i});
      }
      my $n = scalar keys %hit_refids;
      foreach $i (keys %hit_refids) {
        $refid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_refids = ();
    }
    $hit_refids{$refid}++;
    $last_id = $id;

  } #### alignment section
}

   #### last reads
    if    ($last_id) {
      #### for this read pair, if both R1 and R2 hit a refid (t1) and only R1 or R2 hit another refid (t2), skip t2 
      my @t = sort {$b <=> $a} values %hit_refids;
      my $m = $t[0]; #### the refid with most 
      foreach $i (keys %hit_refids) {
        delete($hit_refids{$i}) unless ($m == $hit_refids{$i});
      }
      my $n = scalar keys %hit_refids;
      foreach $i (keys %hit_refids) {
        $refid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_refids = ();
    }

if ($sam_file ne "-") { close(TMP); }

#####


my %refid_depth_of_cov = ();
foreach $refid (keys %refid_reads_count) {
  $refid_depth_of_cov{$refid} = $read_length * $refid_reads_count{$refid} / $refid_2_len{$refid};
}

open(OUT, "> $output") || die "can not write to $output";
print OUT "RefID\tReads\tDepth\tRef_length\tDescription\n";
my @refids = keys %refid_depth_of_cov;
   @refids = sort { $refid_depth_of_cov{$b} <=> $refid_depth_of_cov{$a} } @refids;
foreach $refid (@refids) {
  print OUT "$refid\t$refid_reads_count{$refid}\t$refid_depth_of_cov{$refid}\t$refid_2_len{$refid}\t$refid_2_des{$refid}\n";
}
close(OUT);


sub usage {
<<EOD
Given a reads to human genome mapping file in SAM format,
this script calculate depth of each reference sequence.
The input SAM may need to be filtered and pre-procssed.

usage:
  $0  -i input_SAM -a human_genome.ann -o output

  options
    -i input sam file, default STDIN
    -a path to one of the bwa index files for reference genome, the one with .ann suffix
    -o output file
    -r read length, default 300, which is 150x2
EOD
}
########## END usage

