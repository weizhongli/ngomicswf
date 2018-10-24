#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================
## sort one or more fasta file by length

my @faa = ();
my $rm = 0;
my $suf = "sorted";
my $decreasing = 1;
while($arg = shift) {
  if    ($arg eq "-r") { $rm = 1;}
  elsif ($arg eq "-s") { $suf = shift; }
  elsif ($arg eq "-d") { $decreasing = shift; }
  else { push(@faa, $arg); }
}

if (not @faa) {
  print <<EOD;

This script sorts one or more fasta files by length. sorted file will be saved to input.$suf files

  options:
    -r replace the old file
    -s suffix, default "sorted"
    -d decreasing order, default 1

  examples:
    $0 input.fa
    $0 input1.fa input2.fa input3.fa ...

EOD
  exit();
}

foreach $faa (@faa) {
  my @records = ();

  open(TMP, $faa) || die "can not open $faa";
  
  my $len = 0;
  my $seq = "";
  my $des = "";
  while($ll = <TMP>){
    if ($ll =~ /^>/) {
      if ($len) {
        push(@records, [$len, $des, $seq]);
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
  close(TMP);
      if ($len) {
        push(@records, [$len, $des, $seq]);
      }
  if ($decreasing) {
    @records = sort {$b->[0] <=> $a->[0]} @records;
  }
  else {
    @records = sort {$a->[0] <=> $b->[0]} @records;
  }

  my $out = "$faa.$suf";
  open(OUT, "> $out") || die "can not write to $out";
  foreach $a (@records) {
    print OUT $a->[1], $a->[2];
  }
  close(OUT);

  if ($rm) {
    my $cmd = `mv -f $out $faa`;
  }
  print STDERR "$faa sorted\n";
}
