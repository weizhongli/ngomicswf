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
die usage() unless ($opts{i} and $opts{o});

my $seq_file   = $opts{i};
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

my $seqid;
my $len = 0;
my $seq = "";
while($ll = <$fh>){
  if ($ll =~ /^>/) {
    print $fho "$seqid\t$len\n" if ($len>0);
    $len = 0;
    $seq = "";
    $seqid = substr($ll,1);
    chop($seqid);
    $seqid =~ s/\s.+$//;
  }
  else {
    $seq .= $ll;
    chop($ll);
    $ll =~ s/\s//g;
    $len += length($ll);
  }
}
    print $fho "$seqid\t$len\n" if ($len>0);

close(TMP) unless ($seq_file eq "-");
close(OUT) unless ($output   eq "-");

