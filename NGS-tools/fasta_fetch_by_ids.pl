#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================


use Getopt::Std;
getopts("i:s:o:a:",\%opts);
die usage() unless ($opts{i} and $opts{s});

my $taxid_file = $opts{i};
my $seq_file   = $opts{s};
my $id_str     = $opts{a};
my $output     = $opts{o};

my %taxid1 = ();
my $taxid;

open(TMP, $taxid_file) || die;
while($ll = <TMP>) {
  chop($ll);
  $ll =~ s/\s.+$//;
  $ll =~ s/^>//;
  $taxid1{$ll} = 1;
}
close(TMP);

my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $OUT";
  $fh = "OUT";
}
my $flag = 0;
open(TMP, $seq_file) || die;
while($ll = <TMP>) {
  if ($ll =~ /^>/) {
    if ($id_str) {
      if ($ll =~ /$id_str\|(\d+)/) {
        $taxid = $1;
      }
    }
    else {
      $taxid = substr($ll,1);
      chop($taxid);
      $taxid =~ s/\s.+$//;
    }
    $flag = ( $taxid1{$taxid} ) ? 1 : 0;
  }
  print $fh $ll  if ($flag);
}
close(TMP);

close(OUT) if ($output);
