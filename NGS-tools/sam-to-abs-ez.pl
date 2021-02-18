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
getopts("i:o:a:t:r:N:c:P:f:",\%opts);
die usage() unless ($opts{o} and $opts{a});


my $sam_file    = $opts{i}; ### input sam file, usually pipe in from samtools view input.bam
   $sam_file    = "-" unless ($sam_file);
my $bwa_db_ann  = $opts{a}; ### bwa db .ann file
my $output      = $opts{o};
my $read_length = $opts{r}; #### read length
   $read_length = 150*2 unless ($read_length);
my $cutoff      = $opts{c}; #### relative abundance cutoff
   $cutoff      = 1e-4 unless ($cutoff);
my $frac_bp     = $opts{f};

my $num_total_reads = $opts{N};
my $num_unmapped_reads = 0;

my ($i, $j, $k, $ll, $cmd);

my %rid_2_frac = ();
if ($frac_bp) {
  open(TMP, "$frac_bp") || die "can not open $frac_bp";
  while($ll=<TMP>){
    my ($rid, $frac) = split(/\s+/, $ll);
    $rid_2_frac{$rid} = $frac;
  }
  close(TMP);
}

#### read genome length
my %rid_2_len = ();
my %rid_2_des = ();
open(TMP, $bwa_db_ann) || die "can not open $bwa_db_ann";
$ll=<TMP>;

#2557576839 28773 11
#0 NZ_M10917.1 Bacillus thuringiensis plasmid unnamed, complete sequence
#0 2955 0

while($lla=<TMP>){
  $llb=<TMP>;
  chop($lla);
  @lls = split(/\s+/, $lla, 3);
  my $rid = $lls[1];
  my $des = $lls[2];
  my $len = (split(/\s+/, $llb))[1];
  $rid_2_len{$rid} = $len;
  $rid_2_des{$rid} = $des;
}
close(TMP);


my $fh;
if ($sam_file eq "-") { $fh = "STDIN";}
else {
  open(TMP, $sam_file) || die "can not open $sam_file";
  $fh = "TMP";
}

my %rid_reads_count = ();
my $last_id = "";
my %hit_rids = ();
my $num_mapped_reads = 0;
while($ll=<$fh>){
  if ($ll =~ /^\@/) { #### headers
    next;
  }
  else { #### alignment
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $id  = $lls[0];
    my $rid = $lls[2];    if ($rid eq "*") {  next; }

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
      #### for this read pair, if both R1 and R2 hit a tid (t1) and only R1 or R2 hit another tid (t2), skip t2 
      my @t = sort {$b <=> $a} values %hit_rids;
      my $m = $t[0]; #### the tid with most 
      foreach $i (keys %hit_rids) {
        delete($hit_rids{$i}) unless ($m == $hit_rids{$i});
      }
      my $n = scalar keys %hit_rids;
      foreach $i (keys %hit_rids) {
        $rid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_rids = ();
    }
    $hit_rids{$rid}++;
    $last_id = $id;

  } #### alignment section
}

   #### last reads
    if    ($last_id) {
      #### for this read pair, if both R1 and R2 hit a tid (t1) and only R1 or R2 hit another tid (t2), skip t2 
      my @t = sort {$b <=> $a} values %hit_rids;
      my $m = $t[0]; #### the tid with most 
      foreach $i (keys %hit_rids) {
        delete($hit_rids{$i}) unless ($m == $hit_rids{$i});
      }
      my $n = scalar keys %hit_rids;
      foreach $i (keys %hit_rids) {
        $rid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_rids = ();
    }

if ($sam_file ne "-") { close(TMP); }

#####
if (not defined($num_total_reads)) {
  $num_total_reads = $num_mapped_reads;
  $num_unmapped_reads = 0; 
}
elsif ($num_total_reads <  $num_mapped_reads) {
  $num_total_reads = $num_mapped_reads;
  $num_unmapped_reads = 0; 
}
else {
  $num_unmapped_reads = $num_total_reads - $num_mapped_reads;  
}




#### apply relative abundance cutoff 
my %rid_depth_of_cov = ();
my $total_depth = 0;
foreach $rid (keys %rid_reads_count) {
  $rid_depth_of_cov{$rid} = $read_length * $rid_reads_count{$rid} / $rid_2_len{$rid};
  $total_depth += $rid_depth_of_cov{$rid};
}

my %rid_relative_abs = ();
foreach $rid (keys %rid_reads_count) {
  $rid_relative_abs{$rid} = $rid_depth_of_cov{$rid} /  $total_depth;
}


my @rids = %rid_reads_count;
   @rids = sort {$rid_relative_abs{$b} <=> $rid_relative_abs{$a}} @rids;

open(OUT, "> $output") || die "can not write to $output";
print OUT "Ref_id\tRef_length\tDescription\tNum_reads\tDepth\tAbundance\tFraction\n";
foreach $rid (@rids) {
  my $frac = 0;
  $frac = int($rid_2_frac{$rid} / $rid_2_len{$rid}*100000)/100000 if ($frac_bp and $rid_2_len{$rid});
  next unless ( $rid_relative_abs{$rid}>= $cutoff);
  print OUT "$rid\t$rid_2_len{$rid}\t$rid_2_des{$rid}\t$rid_reads_count{$rid}\t$rid_depth_of_cov{$rid}\t$rid_relative_abs{$rid}\t$frac\n";
}

close(OUT);


sub usage {
<<EOD
Given a reads to reference genome mapping file in SAM format,
this script generate reference abundance profile.
The input SAM need to be filtered and pre-procssed.

usage:
  $0  -i input_SAM -a bwa_ref_genome_db.ann -t taxon_info_file -o output

  options
    -i input sam file, default STDIN
    -a path to one of the bwa index files for reference genome, the one with .ann suffix
    -o output file
    -r read length, default 150
    -N total number of reads, used for calculating unmapped reads
    -c relative abundance cutoff, default 1e-9
    -f fraction of coverage file
EOD
}
########## END usage

