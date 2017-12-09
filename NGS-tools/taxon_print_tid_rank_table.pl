#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);

use Getopt::Std;

getopts("i:o:d:t:s:r:n:l:f:a:",\%opts);
die usage() unless ($opts{o} and $opts{d} and $opts{f});

my $output              = $opts{o};
my $tax_dir             = $opts{d};
my $fasta_file          = $opts{f}; #### fasta file
my $tid_str             = $opts{t}; #### fasta header with $tid_str|taxid
   $tid_str             = "taxid" unless ($tid_str);
my $des_str             = $opts{a}; #### fasta header with $des_str|taxid_description
   $des_str             = "taxdes" unless ($des_str);

my ($i, $j, $k, $ll, $cmd);
my @ranks = qw/superkingdom kingdom phylum class order family genus species toprank/;
my %ranks = map {$_, 1} @ranks;

my %toprank_taxids = ();
my %taxid_2_name=();

open(TMP, $fasta_file) || die "Can not open $fasta_file";
while($ll=<TMP>){
  next unless ($ll =~ /^>/);
  if ($ll =~ /$tid_str\|(\d+)/) {
    my $taxid = $1;
    $toprank_taxids{$taxid} = 1;
    if ($ll =~ /$des_str\|(\w+)/) {
      my $t1 = $1; $t1 =~ s/_/ /g;
      $taxid_2_name{$taxid} = $t1;
    }
  }
}
close(TMP);


$taxid_2_name{"-1"} = "unknown";

#### read tax name
#### re-read tax name
open(TMP, "$tax_dir/names.dmp") || die "Can not open $tax_dir/names.dmp\n";
while($ll=<TMP>){
  chop($ll);
  next unless ($ll =~ /scientific name/);
  my ($t_id, $name, @t) = split(/\|/, $ll);
  $t_id  =~ s/\s+//g;
  $name  =~ s/^\s+//;
  $name  =~ s/\s+$//;
  #$name  =~ s/\W+/_/g;  ## no white space
  $taxid_2_name{$t_id} = $name;
}
close(TMP);

#### read tax id tree
my %taxid_2_parent = ();
my %taxid_2_rank = ();
my %taxid_2_children = ();
open(TMP, "$tax_dir/nodes.dmp") || die "Can not open $tax_dir/nodes.dmp\n";
while($ll=<TMP>){
  chop($ll);
  my ($t_id, $parent_id, $rank, @t) = split(/\|/, $ll);
  $t_id      =~ s/\s+//g;
  $parent_id =~ s/\s+//g;
  $rank      =~ s/\s+//g;

  #### for viruses, since it don't have phylum, treat no rank under 10239 as phylum
  $rank = "phylum" if (($parent_id eq "10239") and ($rank eq "norank"));
  $taxid_2_parent{$t_id} = $parent_id;
  $taxid_2_rank{$t_id}   = $rank;
  if (not defined($taxid_2_children{$parent_id})) {
    $taxid_2_children{$parent_id} = [];
  }
  push(@{$taxid_2_children{$parent_id}}, $t_id);
}
close(TMP);

####parse tax tree
my %taxid_involved = ();
my %toptax_rank_2_taxid = ();
foreach $ttax (keys %toprank_taxids) {
  $taxid = $ttax;
  $taxid_involved{$taxid} = 1;
  $toptax_rank_2_taxid{$ttax} = {};
  $toptax_rank_2_taxid{$ttax}{"toprank"} = $ttax;
  while(1) {
    $tax_rank = $taxid_2_rank{$taxid};
    foreach $r (@ranks) {
      if ($r eq $tax_rank) {
        $toptax_rank_2_taxid{$ttax}{$r} = $taxid;
        $taxid_involved{$taxid} = 1;
      }
    }
    if ($taxid_2_parent{$taxid}>1) {
      $taxid = $taxid_2_parent{$taxid};
    } #### > root id
    else {
      last;
    }
  } #### while
}

####also get toptax_rank_2_taxid for other ranks
foreach $ttax (keys %taxid_involved) {
  next if (defined($toptax_rank_2_taxid{$ttax}));
  $taxid = $ttax;
  $toptax_rank_2_taxid{$ttax} = {};
  while(1) {
    $tax_rank = $taxid_2_rank{$taxid};
    foreach $r (@ranks) {
      if ($r eq $tax_rank) {
        $toptax_rank_2_taxid{$ttax}{$r} = $taxid;
      }
    }
    if ($taxid_2_parent{$taxid}>1) {
      $taxid = $taxid_2_parent{$taxid};
    } #### > root id
    else {
      last;
    }
  } #### while
}

open(OUT, "> $output") || die "can not write to $output\n";
print OUT "#taxid\trank\tname";
foreach $i (@ranks) { print "\t$i\t$i\_ti"; }
print OUT "\n";
if (1) {
  my @top_tax_list = keys %toprank_taxids;
  for ($i=0; $i<@ranks; $i++) {
    my $rank = $ranks[$i];
    my %current_taxids = ();

    foreach $ttax (@top_tax_list) {
      my $taxid = $toptax_rank_2_taxid{$ttax}{$rank};
      if (not $taxid) {
         $taxid = -1;
         for (my $i0 = $i; $i0>=0; $i0--) {
           my $j0 = $toptax_rank_2_taxid{$ttax}{$ranks[$i0]};
           if ($j0) {
             $taxid = -$j0;
             $taxid_2_name{$taxid} = "Unnamed $rank in $taxid_2_name{$j0} $ranks[$i0]";
             last;
           }
         }
      }
      $current_taxids{$taxid} = 1;
    }
    my @ref_tax_list = sort keys %current_taxids;
    foreach my $taxid (@ref_tax_list) {
      print OUT "$taxid\t$rank\t$taxid_2_name{$taxid}";
      for ($j=0; $j<= $#ranks; $j++) { 
        my $rank_tax =  ($taxid > 0) ? $toptax_rank_2_taxid{$taxid}{$ranks[$j]} :
                                       $toptax_rank_2_taxid{-$taxid}{$ranks[$j]}; 
        my $rank_name = $taxid_2_name{$rank_tax}; 
           $rank_name = "\\N" unless (defined($rank_name));
        print OUT "\t$rank_name";
        $rank_tax = "\\N" unless (defined($rank_tax));
        print OUT "\t$rank_tax";
      }
      print OUT "\n";
    }
  }
}


sub usage {
<<EOD
Given a fasta file and path to NCBI taxonomy dir, this script produce a table of 
taxid, rank and tree
usage:
  $script_name  -o output -f path_to_reference_fasta_file -d path_to_ncbi_taxonomy_dir

  options
    -t tid_str, taxid need to be coded in fasta header in format like taxid|12345
       or ti|12345. here 12345 is the taxid and 'taxid' is the tid_str
       default: taxid
    -a des_str, similar to tid_str, tax name may be included in fasta header in format
       like taxdes|species_xyz
       default: taxdes

EOD
}
########## END usage





