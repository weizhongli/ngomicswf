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

my $seqid_file = $opts{i};
my $seq_file   = $opts{s};
my $output     = $opts{o};

my %seqid1 = ();
my $seqid;

open(TMP, $seqid_file) || die;
while($ll = <TMP>) {
  chop($ll);
  $ll =~ s/\s.+$//;
  if ($ll =~ /(PF\d+)/) {
    $seqid1{$1} = 1;
  }
}
close(TMP);

my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $OUT";
  $fh = "OUT";
}
my $flag = 0;
if ($seq_file =~ /.gz$/) { open(TMP, "gunzip < $seq_file | ") || die "can not gunzip $seq_file"; }
else                     { open(TMP, $seq_file              ) || die; }

my $model = "";
my $acc = "";
while($ll = <TMP>) {
  $model .= $ll;
  if ($ll =~ /^\/\//) {
    if (defined($seqid1{$acc})) {
      print $fh $model;
    }
    $model = "";
    $acc = "";
  }
  if ($ll =~ /^ACC\s+(PF\d+)/) {
    $acc = $1;
  }
}
close(TMP);
close(OUT) if ($output);

sub usage {
<<EOD;
Fetch subset of HMMs by Accessions from HMM database
Version number of accessions will be ignored, for example PF01234.x will match PF01234.y

$script_name -i HMM_accession_file -s original_HMM file -o output

  -i input file of HMM accessions,  one accession per line, accession must be the first column, 
     - for STDIN
  -s original HMM file
  -o output, default STDOUT

EOD
}

