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
getopts("i:j:o:p:",\%opts);
die usage() unless ($opts{o} and $opts{p});

my $input      = $opts{i};
my $output_R1  = $opts{o};
my $output_R2  = $opts{p};

my $fh = "STDIN";
if ($input) {
  open(TMP, $input) || die "Can not open $input";
  $fh = "TMP";
}

open(OUT1, "> $output_R1") || die "can not write to $output_R1";
open(OUT2, "> $output_R2") || die "can not write to $output_R2";

$line = 0;
while($ll = <$fh>) {
  $i = $line % 8;
  if ($i<=3) {
    print OUT1 $ll;
  }
  else {
    print OUT1 $ll;
  }
  $line++;
}
close(TMP) if ($input);

close(OUT1);
close(OUT2);

sub usage {
<<EOD;
split a interlaced fastq file into seperate fq file

$0 -o output_R1 -q output_R2 < input_interlaced_fq
  -i input interlaced fastq, default STDIN
  -o output_R1
  -p output_R2
EOD

}
