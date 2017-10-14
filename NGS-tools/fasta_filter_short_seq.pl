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

use Getopt::Std;
getopts("i:s:o:a:b:c:",\%opts);
die usage() unless ($opts{i} and defined($opts{c}) and $opts{o});

my $seq_file   = $opts{i};
my $len_cutoff = $opts{c};
my $output     = $opts{o};
my $ann        = $opts{a}; #### annotate fasta header with |$ann|length_of_seq

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
    if ($len >= $len_cutoff) {
      $des = ann_str($des, $ann, $len) if $ann;
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
    if ($len >= $len_cutoff) {
      $des = ann_str($des, $ann, $len) if $ann;
      print $fho $des, $seq;
    }

close(TMP) unless ($seq_file eq "-");
close(OUT) unless ($output   eq "-");

sub ann_str {
  my ($des, $ann, $len) = @_;
  if ($des =~ /\s/) {
    $des =~ s/(\s)/\|$ann|$len$1/;
  }
  else {
    $des = "$des|$ann|$len";
  }
  return $des;
}

sub usage {
<<EOD;
filter out short sequences
$script_name -i input -o output -c length_cutoff
  -i input, - for stdin
  -o output, - for stdout
  -c cutoff
  -a annotation string, annotate fasta header
     with -a LEN, then "|LEN|sequence_length" will be added before the first white space
EOD
}
