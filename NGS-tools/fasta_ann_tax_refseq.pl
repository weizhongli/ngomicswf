#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================
## 

use Getopt::Std;

getopts("i:o:d:f:",\%opts);
die usage() unless ($opts{i} and $opts{d});


my $fasta           = $opts{i};  #
my $refseq_catalog  = $opts{d};  # RefSeq-release66.catalog
my $output          = $opts{o};  $output = "$fasta.t" unless (defined($output));

my %accs          = ();
open(TMP, $fasta) || die "Can not open $fasta\n";
while($ll=<TMP>){
  if ($ll =~ /^>(\S+)/) {
    $accs{$1} = 1;
  }
}
close(TMP);

#### read tax name
my %taxid_2_name = ();
my %acc_2_taxid = ();
open(TMP, "$refseq_catalog") || die "Can not open $refseq_catalog\n";
while($ll=<TMP>){
  chop($ll);
  my ($t_id, $name, $acc, $div, $status, $len) = split(/\t/, $ll);

  next unless ($accs{$acc});
  $acc_2_taxid{$acc} = $t_id;

  next if (defined($taxid_2_name{$t_id}));
  $name  =~ s/\W+/_/g;  ## no white space
  $taxid_2_name{$t_id} = $name;
}
close(TMP);


my $fh = "OUT";
open($fh, "> $output") || die "can not write to $output";
open(TMP, $fasta) || die "Can not open $fasta\n";
while($ll=<TMP>){
  if ($ll =~ />(\S+)/) {
    $taxid = $acc_2_taxid{$1};
    $taxdes = $taxid_2_name{$taxid}; 
    $taxdes = "none" unless ($taxdes);
    if ($taxid) {
      $ll =~ s/>(\S+)/>$1|taxid|$taxid|taxdes|$taxdes/;
    }
  }
  print $fh $ll;
}
close(TMP);
close($fh);

sub usage {
<<EOD;

This script add tax IDs and tax name to the headers of Refseq fasta file
so that the header will looks like ">NC_002483.1|taxid|83333|taxdes|Escherichia_coli_K-12 ..."

    options:
    -i fasta iutput
    -d refseq release catalog file
    -o output file name, default input.t

EOD
}
########## end usage
