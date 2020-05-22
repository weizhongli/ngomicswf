#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================


## copied from  assembly-binning.pl 
## added taxon parsing function
## parse from taxid and get species taxid from taxon file
use Getopt::Std;
getopts("i:j:c:s:o:p:d:r:n:a:t:q:m:N:l:r:x:X:f:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{o} and $opts{t} and $opts{a});

my $centrifuge_in    = $opts{i};
my $centrifuge_report= $opts{j};
my $bwa_db_ann       = $opts{a}; ### bwa db .ann file
my $hit_len_cutoff   = $opts{l};
   $hit_len_cutoff   = 60 unless ($hit_len_cutoff > 0);
my $taxon_file       = $opts{t}; #### taxon file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl
my $read_length = $opts{r}; #### read length
   $read_length = 150*2 unless (defined($read_length));
my $cutoff           = $opts{c}; #### relative abundance cutoff
   $cutoff           = 1e-6 unless ($cutoff);
my $output           = $opts{o};
my $exclude_f        = $opts{x};
my $exclude_s        = $opts{X};
my $cutoff_evenness  = $opts{f};
   $cutoff_evenness  = 0.1 unless ($cutoff_evenness);
my $out_abs          = "$output-abundance";
my $out_genome       = "$output-genome";
my $out_species      = "$output-species";
my $num_total_reads  = $opts{N};
my $num_unmapped_reads = 0;
my $log_f            = "$output-log";
my $output_mapped    = "$output-mapped-reads.ids";

my ($i, $j, $k, $ll, $cmd);

open(LOG, ">$log_f") || die "can not write to $log_f";

#### read genome length
#0 acc|NZ_DS480350.1|taxid|428125|sptid|1535|getid|N [Clostridium] leptum DSM 753 Scfld_02_19, whole genome shotgun sequence
my %tid_2_len = ();
my %acc_2_len = ();
my %tid_2_accs = ();
open(TMP, $bwa_db_ann) || die "can not open $bwa_db_ann";
$ll=<TMP>;
while($ll=<TMP>){
  if ($ll =~ /taxid\|(\d+)/) {
    my $tid = $1;
    if (not defined($tid_2_accs{$tid})) {
      $tid_2_accs{$tid} = [];
    }
    my $acc = "";
    if ($ll =~ /(acc\|[^\|]+)/) {
      $acc = $1;
    }
    $ll=<TMP>;
    my $len1 = (split(/\s+/, $ll))[1];
    $tid_2_len{$tid} += $len1;
    if ($acc) {
      $acc_2_len{$acc} = $len1;
      push(@{ $tid_2_accs{$tid} }, $acc);
    }
  }
}
close(TMP);


my %exclude_tids = ();
if ($exclude_f and (-e $exclude_f)) {
  open(TMP, $exclude_f) || die "can not open $exclude_f";
  while($ll=<TMP>) {
    chop($ll);
    next if ($ll =~ /^#/);
    my $tid = (split(/\s+/, $ll))[0];
    next unless ($tid =~ /^\d+$/);
    $exclude_tids{$tid} = 1;
  }
  close(TMP);
}
if ($exclude_s) {
  my @lls = split(/,/, $exclude_s);
  foreach my $tid (@lls) {
    next unless ($tid =~ /^\d+$/);
    $exclude_tids{$tid} = 1;
  }
}

#### open centrifuge report first time
#### to get list of valid taxids, these taxids is after EM calc, so excluding spurious hits
#### that maybe in $centrifuge_in
my %em_valid_taxids = ();
open(TMP, $centrifuge_report) || die "can not open $centrifuge_report";
while($ll=<TMP>) {
  $ll =~ s/^\s+//;
  my ($abs, $t1, $t2, $rank, $tid, $name) = split(/\s+/, $ll);
  next unless ($tid =~ /^\d+$/); #### is a valid taxid
  $em_valid_taxids{$tid} = 1;
}
close(TMP);


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

#### based on sam-to-taxon-abs-ez.pl

open(TMP, $centrifuge_in) || die "can not open $centrifuge_in";
my @reads_mapped = ();
my %tid_reads_count = ();
my %tid_bases_count = ();
my %acc_reads_count = ();
my %acc_bases_count = ();
my $last_id = "";
my %hit_tids = ();
my $num_mapped_reads = 0;
my ($id,  $rid, $tid, $score, $score2, $hit_len, $len, $nhit);
$ll = <TMP>;
while($ll=<TMP>) {
  chop($ll);
  ($id,  $rid, $tid, $score, $score2, $hit_len, $len, $nhit) = split(/\t/, $ll);
  next if ($acc eq "unclassified");
  next unless ($hit_len >= $hit_len_cutoff);
  next unless ($em_valid_taxids{$tid}); #### avoid non EM passed taxon, for hits with ties
  next unless (defined( $taxon_info{$tid}) );
  next unless (defined( $tid_2_len{$tid}));
  next if     (defined( $exclude_tids{$tid}));

  #### $rid can be "species" "genus"
  $acc_reads_count{$rid}++     if( $acc_2_len{$rid} );
  $acc_bases_count{$rid}+=$len if( $acc_2_len{$rid} );
  if    (($id    ne $last_id) and $last_id ) {
    my $n = scalar keys %hit_tids;
    foreach $i (keys %hit_tids) {
      $tid_reads_count{$i} += 1/$n;
      $tid_bases_count{$i} += $len/$n;
    }
    push(@reads_mapped, [$last_id, keys %hit_tids]);
    $num_mapped_reads++;
    %hit_tids = ();
  }

  $hit_tids{$tid}++;
  $last_id = $id;
}
  if    ($last_id) {
    my $n = scalar keys %hit_tids;
    foreach $i (keys %hit_tids) {
      $tid_reads_count{$i} += 1/$n;
      $tid_bases_count{$i} += $len/$n;
    }
    push(@reads_mapped, [$last_id, keys %hit_tids]);
    $num_mapped_reads++;
    %hit_tids = ();
  }
close(TMP);
#### finish reading centrifuge input


#### evenness filtering
my %failed_evenness_tids = ();
foreach $tid (keys %tid_reads_count) {
  if (not defined($tid_2_len{$tid})) {
    print LOG "not defind $tid\n"; next;
  }
  my $dc = $read_length * $tid_reads_count{$tid} / $tid_2_len{$tid};
  $dc    =                $tid_bases_count{$tid} / $tid_2_len{$tid} if ($read_length == 0);
  next unless ($dc > 0);

  my $expected_fc = 1 - 1 / exp($dc);
  my $observed_fc = 0;
  my $max_d = 0;
  my $max_acc = "";
  foreach $acc (@{ $tid_2_accs{$tid} }) {
    next unless (defined($acc_reads_count{$acc}));
    next unless (defined($acc_2_len{$acc}));
    ## I don't have mapping coordinate, assume the max coverage on this acc =
    my $total1 = $read_length * $acc_reads_count{$acc}; 
       $total1 =                $acc_bases_count{$acc} if ($read_length == 0);
    if ($total1 / $acc_2_len{$acc} > $max_d) {
      $max_d = $total1 / $acc_2_len{$acc};
      $max_acc = $acc;
    }
    ## but max acc should < length of this acc
    $total1 = $acc_2_len{$acc} if ($acc_2_len{$acc} < $total1 );
    $observed_fc += $total1;
  }
  $observed_fc /= $tid_2_len{$tid};

  ## f can be zero if tid got hits from only "species" "genus", without exact accessions
  if (($observed_fc / $expected_fc) < $cutoff_evenness) {
    my $n = $taxon_info{$tid}->[2];
    print LOG "$tid ($n, $tid_2_len{$tid} bp) evenness error! reads $tid_reads_count{$tid}, depth $dc, observed/expected fraction $observed_fc / $expected_fc\n";
    print LOG "\tf= ",  ($observed_fc/$expected_fc), ", max contig depth $max_d at $max_acc ($acc_2_len{$max_acc}bp, $acc_reads_count{$max_acc} reads)\n";
    delete($tid_reads_count{$tid});
    delete($tid_bases_count{$tid});
    $failed_evenness_tids{$tid} = 1;
  }
}

open(OUTR, "> $output_mapped") || die "can not write to $output_mapped";
foreach $i (@reads_mapped) {
  my ($id, @tids) = @{ $i };
  my @mapped_tids = ();
  foreach $j (@tids) {
    push(@mapped_tids, $j) unless (defined($failed_evenness_tids{$j} ));
  }
  if (@mapped_tids) {
    print OUTR "$i\t", join("|", @mapped_tids), "\n";
  }
}
close(OUTR);

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
  if (not defined($tid_2_len{$tid})) {
    print LOG "not defind $tid\n"; next;
  }
  $tid_depth_of_cov{$tid} = $read_length * $tid_reads_count{$tid} / $tid_2_len{$tid};
  $tid_depth_of_cov{$tid} =                $tid_bases_count{$tid} / $tid_2_len{$tid} if ($read_length == 0);

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

#taxid  rank    name    superkingdom    superkingdom_ti kingdom kingdom_ti  phylum  phylum_ti   class   class_ti    order   order_ti    family  family_ti   genus   genus_ti    species species_ti  toprank toprank_ti
#1000373    toprank Rosellinia necatrix quadrivirus 1   Viruses 10239   \N  \N  dsRNA viruses   35325   \N  \N  \N  \N  Quadriviridae   1299296 Quadrivirus 1299297 Rosellinia necatrix quadrivirus 1   1000373 Rosellinia necatrix quadrivirus 1   1000373
my @ranks = qw/phylum class order family genus species toprank/;
my %rank_col = qw/phylum 7 class 9 order 11 family 13 genus 15 species 17 toprank 19/;
foreach $rank (@ranks) {

  my %rank_ti_info = ();
  my %rank_ti_abs = ();
  my %rank_ti_cov = ();
  my %rank_ti_reads = ();
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
    $rank_ti_reads{$rank_tid} += $tid_reads_count{$tid};
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
  print OUT "\trelative_abundance\tdepth_coverage\tnum_reads\n";
  foreach $rank_tid (@output_rank_tids) {
    print OUT $rank_tid, "\t";
    print OUT join("\t", @{$rank_ti_info{$rank_tid}});
    print OUT "\t$rank_ti_abs{$rank_tid}\t$rank_ti_cov{$rank_tid}\t$rank_ti_reads{$rank_tid}\n";
  }
  close(OUT);
}
close(LOG);

sub usage {
<<EOD
Given centrifuge output, 
this script generate taxonomy abundance profile.

usage:
  $0  -i centrifuge_output -j centrifuge_report -m acc_to_taxid_file -a taxid_to_length_file -t taxon_info_file -o output

  options
    -i input, centrifuge output
    -j input, centrifuge report
    -a reference taxid to length file: .ann file made by bwa index
    -r read length, default 300 bp
       use 0 if your read length vary a lot, e.g. Nanopore reads 
    -l centrifuge hit length cutoff, default 60
    -c abundance cutoff, default 1e-6
    -t taxon info file for reference genomes, 
       created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl based on blast ref db
    -o output prefix
    -x filename of list of taxids to exclude, this can be used to filter out known contaminant species
    -f coverage evenness cutoff, default 0.1
       given depth of coverage of a genome (dc), the expected fraction of coverage (fc) =  1 - 1 / exp(dc)  (Poisson distribution)
       if (observed fc / expected fc) < this cutoff, it suggests aligned reads are piled up at a small fraction, very likely due to
       artifacts.  Ref: Lander et al. doi:10.1016/0888-7543(88)90007-9

EOD
}
########## END usage

