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
getopts("i:s:o:a:b:c:",\%opts);
die usage() unless ($opts{i} and $opts{c} and $opts{o});

my $seq_file   = $opts{i};
my $len_cutoff = $opts{c};
my $output     = $opts{o};

die usage () unless ($len_cutoff > 0); #### at least = 1

my $fh;
if ($seq_file eq "-") { $fh = "STDIN";}
else {
  open(TMP, $seq_file) || die "can not open $seq_file";
  $fh = "TMP";
}

my $fho;
if ($output eq "-") { $fho = "STDOUT";}
else {
  open(OUT, ">$output") || die "Can not write to $output";
  $fho = "OUT";
}

my $len = 0;
my $seq = "";
my $des = "";
while($ll = <$fh>){
  if ($ll =~ /^>/) {
    if ($len <= $len_cutoff) {
      print $fho $des, $seq;
    }
    $len = 0;
    $seq = "";
    $des = $ll;
  }
  else {
    $seq .= $ll;
    chop($ll);
    $ll =~ s/\s//g;
    $len += length($ll);
  }
}
    if ($len <= $len_cutoff) {
      print $fho $des, $seq;
    }

close(TMP) unless ($seq_file eq "-");
close(OUT) unless ($output   eq "-");

