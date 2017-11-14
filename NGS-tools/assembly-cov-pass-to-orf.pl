#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

## given coverage of assembled contig
## pass the cov to ORFs predicted from contigs

use Getopt::Std;

getopts("i:o:d:t:s:",\%opts);
die usage() unless ($opts{i} and $opts{o} and $opts{d});

my $orf_fasta    = $opts{i};
my $contig_cov   = $opts{d};
my $output_file  = $opts{o};

my $contig_2_cov = ();
open(TMP, $contig_cov) || die "can not open $contig_cov\n";
while($ll=<TMP>){
  chop($ll);
  next if ($ll =~ /^#/); #### skip comments
  my ($id, $cov) = split(/\t/, $ll);
  $contig_2_cov{$id} = $cov;
}
close(TMP);


open(TMP, $orf_fasta) || die "can not open $orf_fasta\n";
open(OUT, "> $output_file") || die "can not write $output_file\n";

while($ll = <TMP>){
  if ($ll =~ /^>(\S+)/) {
    $ORF = $1;
    $contig = $ORF; 
    $contig =~ s/_\d+$//;
    if ($contig_2_cov{$contig}) {
      print OUT "$ORF\t$contig_2_cov{$contig}\n";
    }
  }
}
close(OUT);
close(TMP);

