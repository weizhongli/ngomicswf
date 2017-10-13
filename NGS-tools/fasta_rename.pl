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
getopts("i:s:o:a:b:",\%opts);
die usage() unless ($opts{i} and $opts{s} and $opts{o});

my $seq_file   = $opts{i};
my $name_str   = $opts{s};
my $no_start   = $opts{b}; $no_start = 1 unless (defined($no_start));
my $output     = $opts{o};

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

while($ll = <$fh>) {
  if ($ll =~ /^>/) {
    $i = substr($ll,1);
    $ll = ">$name_str$no_start $i";
    $no_start++;
  }
  print $fho $ll;
}

close(TMP) unless ($seq_file eq "-");
close(OUT) unless ($output   eq "-");

