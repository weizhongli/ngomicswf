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
die usage() unless ($opts{o} and $opts{a} and $opts{t});


my $sam_file    = $opts{i}; ### input sam file, usually pipe in from samtools view input.bam
   $sam_file    = "-" unless ($sam_file);
my $bwa_db_ann  = $opts{a}; ### bwa db .ann file
my $output      = $opts{o};
my $taxon_file  = $opts{t}; #### taxon file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl
my $read_length = $opts{r}; #### read length
   $read_length = 150*2 unless ($read_length);
my $cutoff      = $opts{c}; #### relative abundance cutoff
   $cutoff      = 1e-9 unless ($cutoff);
my $num_total_reads = $opts{N};
my $num_unmapped_reads = 0;

my ($i, $j, $k, $ll, $cmd);

#### read genome length
my %tid_2_len = ();
my %tid_2_seqs = ();
open(TMP, $bwa_db_ann) || die "can not open $bwa_db_ann";
$ll=<TMP>;
while($ll=<TMP>){
  if ($ll =~ /taxid\|(\d+)/) {
    my $tid = $1;
    $ll=<TMP>;
    $tid_2_len{$tid} += (split(/\s+/, $ll))[1];
  }
}
close(TMP);



my $fh;
if ($sam_file eq "-") { $fh = "STDIN";}
else {
  open(TMP, $sam_file) || die "can not open $sam_file";
  $fh = "TMP";
}

my %tid_reads_count = ();
my $last_id = "";
my %hit_tids = ();
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
    my $tid = "";
    if ($rid =~ /taxid\|(\d+)/) { $tid = $1; }
    else { next; }

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
      my @t = sort {$b <=> $a} values %hit_tids;
      my $m = $t[0]; #### the tid with most 
      foreach $i (keys %hit_tids) {
        delete($hit_tids{$i}) unless ($m == $hit_tids{$i});
      }
      my $n = scalar keys %hit_tids;
      foreach $i (keys %hit_tids) {
        $tid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_tids = ();
    }
    $hit_tids{$tid}++;
    $last_id = $id;

  } #### alignment section
}

   #### last reads
    if    ($last_id) {
      #### for this read pair, if both R1 and R2 hit a tid (t1) and only R1 or R2 hit another tid (t2), skip t2 
      my @t = sort {$b <=> $a} values %hit_tids;
      my $m = $t[0]; #### the tid with most 
      foreach $i (keys %hit_tids) {
        delete($hit_tids{$i}) unless ($m == $hit_tids{$i});
      }
      my $n = scalar keys %hit_tids;
      foreach $i (keys %hit_tids) {
        $tid_reads_count{$i} += 1/$n;
      }
      $num_mapped_reads++;
      %hit_tids = ();
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

#### read taxon info
my $taxon_format = <<EOD;
Col.0   #taxid  511145
Col.1   rank    toprank
Col.2   name    Escherichia coli str. K-12 substr. MG1655
Col.3   superkingdom    Bacteria
Col.4   superkingdom_ti 2
Col.5   kingdom \\N
Col.6   kingdom_ti      \\N
Col.7   phylum  Proteobacteria
Col.8   phylum_ti       1224
Col.9   class   Gammaproteobacteria
Col.10  class_ti        1236
Col.11  order   Enterobacterales
Col.12  order_ti        91347
Col.13  family  Enterobacteriaceae
Col.14  family_ti       543
Col.15  genus   Escherichia
Col.16  genus_ti        561
Col.17  species Escherichia coli
Col.18  species_ti      562
Col.19  toprank Escherichia coli str. K-12 substr. MG1655
Col.20  toprank_ti      511145
EOD
open(TMP, $taxon_file) || die "can not open $taxon_file";
my %taxon_info = ();
while($ll=<TMP>) {
  chop($ll);
  next if ($ll =~ /^#/);
  $ll =~ s/\\N/NA/g;
  my ($tid, $rank, @lls) = split(/\t/, $ll);
  next unless ($rank eq "toprank");
  $taxon_info{$tid} = [$tid, $rank, @lls];
}
close(TMP);

my $rank = "superkingdom";
my %superkingdom_reads = ();
foreach $tid (keys %tid_reads_count) {
  my @tid_info = @{ $taxon_info{$tid} };
  $superkingdom_ti = $tid_info[4];
  $superkingdom_reads{$superkingdom_ti} += $tid_reads_count{$tid};
}

my $out = "$output.$rank.txt";
#TaxID   superkingdom   
#2       Bacteria       
#2157    Archaea        
#2759    Eukaryota      
#10239   Viruses        
#12884   Viroids 
my %tid_2_des = qw/2 Bacteria 2157 Archaea 2759 Eukaryota 10239 Viruses/;       
open(OUT, "> $out") || die "can not write to $out";
print OUT "#tid\tsuperkingdom\tnumber_reads\n";
foreach $tid (qw/2 2157 2759 10239/) {
  $j = $superkingdom_reads{$tid}; $j=0 unless ($j);
  $j = int($j);
  print OUT "$tid\t$tid_2_des{$tid}\t$j\n";
}
if ($num_unmapped_reads) {
  print OUT "-1\tunmaped\t$num_unmapped_reads\n";
}
close(OUT);


#### apply relative abundance cutoff 
my %tid_depth_of_cov = ();
my %superkingdom_total_cov = ();
foreach $tid (keys %tid_reads_count) {
  $tid_depth_of_cov{$tid} = $read_length * $tid_reads_count{$tid} / $tid_2_len{$tid};
  my @tid_info = @{ $taxon_info{$tid} };
  $superkingdom_ti = $tid_info[4];
  $superkingdom_total_cov{$superkingdom_ti} += $tid_depth_of_cov{$tid};
}

my %tid_relative_abs = ();
my $total_adj_abs = 0;
foreach $tid (keys %tid_reads_count) {
  my @tid_info = @{ $taxon_info{$tid} };
  $superkingdom_ti = $tid_info[4]; 

  my $abs = $tid_depth_of_cov{$tid} / $superkingdom_total_cov{$superkingdom_ti} * 
           ( $superkingdom_reads{$superkingdom_ti} / $num_mapped_reads);
  $abs  = int($abs * 100000000)/100000000;
  if ($abs >= $cutoff) {
    $total_adj_abs +=$abs;
  }
  $tid_relative_abs{$tid} = $abs;
}


#taxid	rank	name	superkingdom	superkingdom_ti	kingdom	kingdom_ti	phylum	phylum_ti	class	class_ti	order	order_ti	family	family_ti	genus	genus_ti	species	species_ti	toprank	toprank_ti
#1000373	toprank	Rosellinia necatrix quadrivirus 1	Viruses	10239	\N	\N	dsRNA viruses	35325	\N	\N	\N	\N	Quadriviridae	1299296	Quadrivirus	1299297	Rosellinia necatrix quadrivirus 1	1000373	Rosellinia necatrix quadrivirus 1	1000373

my @ranks = qw/phylum class order family genus species toprank/;
my %rank_col = qw/phylum 7 class 9 order 11 family 13 genus 15 species 17 toprank 19/;
foreach $rank (@ranks) {

  my %rank_ti_info = ();
  my %rank_ti_abs = ();
  my %rank_ti_cov = ();
  my $rank_col = $rank_col{$rank};

  foreach $tid (keys %tid_reads_count) {
    my $abs = $tid_relative_abs{$tid};
    next unless ($abs >= $cutoff);
    my @tid_info = @{ $taxon_info{$tid} };
    my $abs_adj = $abs / $total_adj_abs;

    my $rank_tid = $tid_info[ $rank_col + 1];
    if ($rank_tid =~ /\d+/) { #### this is good
      ;
    }
    else { #### missing tid at this rank
      $j = $rank_col + 2; 
      while(1) {
        if ($tid_info[$j+1] =~ /\d+/) {
          $rank_tid = $tid_info[$j+1];
          last;
        }
        $j +=2; 
      }
    }
    if (not defined($rank_ti_info{$rank_tid})) {
      my @j = ();
      for ($j=3; $j<=$rank_col{$rank}; $j+=2) {
         push(@j, $j);
      }
      $rank_ti_info{$rank_tid} = [ @tid_info[@j] ];
    }
    $rank_ti_abs{$rank_tid} += $abs_adj;
    $rank_ti_cov{$rank_tid} += $tid_depth_of_cov{$tid};
  }
  my @output_rank_tids = keys %rank_ti_abs;
     @output_rank_tids = sort { $rank_ti_abs{$b} <=> $rank_ti_abs{$a} } @output_rank_tids;

  my $out = "$output.$rank.txt";
  open(OUT, "> $out") || die "can not write to $out";
  print OUT "#tid\tsuperkingdom\tkingdom";
  foreach $i (@ranks) {
    print OUT "\t$i";
    last if ($i eq $rank);
  }
  print OUT "\trelative_abundance\tdepth_coverage\n";
  foreach $rank_tid (@output_rank_tids) {
    print OUT $rank_tid, "\t";
    print OUT join("\t", @{$rank_ti_info{$rank_tid}});
    print OUT "\t$rank_ti_abs{$rank_tid}\t$rank_ti_cov{$rank_tid}\n";
  }
  close(OUT);
}


sub usage {
<<EOD
Given a reads to reference genome mapping file in SAM format,
this script generate taxonomy abundance profile.
The input SAM need to be filtered and pre-procssed.

usage:
  $0  -i input_SAM -a bwa_ref_genome_db.ann -t taxon_info_file -o output

  options
    -i input sam file, default STDIN
    -a path to one of the bwa index files for reference genome, the one with .ann suffix
    -o output file
    -t taxon info file for reference genomes, 
       created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl based on blast ref db
    -r read length, default 150
    -N total number of reads, used for calculating unmapped reads
    -c relative abundance cutoff, default 1e-9
EOD
}
########## END usage

