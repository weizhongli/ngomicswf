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
getopts("i:s:o:a:r:L:",\%opts);
die usage() unless ($opts{i});

my $seq_file   = $opts{i};
my $ref_file   = $opts{r};
my $ref_len    = $opts{L};

if ($ref_file) {
  die usage() if (defined($ref_len));
  if ($ref_file =~ /.gz$/) { open(TMP, "gunzip < $ref_file| ") || die "can not gunzip $ref_file"; }
  else                     { open(TMP,           $ref_file   ) || die; }
  $ref_len = 0;
  while($ll=<TMP>){
    next if ($ll =~ /^>/);
    chop($ll);
    $ll =~ s/\s//g;
    $ref_len += length($ll);
  }
  close(TMP);
}

if ($seq_file =~ /.gz$/) { open(TMP, "gunzip < $seq_file | ") || die "can not gunzip $seq_file"; }
else                     { open(TMP, $seq_file              ) || die; }
if (1) {
  my @lens = ();
  my $sum = 0;
  my $ave = 0;
  my $N50 = 0;
  my $N80 - 0;
  my $N = 0;

  my $len = 0;
  while($ll = <TMP>){
    if ($ll =~ /^>/) {
      if ($len) {
        $N++;
        $sum += $len;
        push(@lens, $len);
      }
      $len = 0;
    }
    else {
      chop($ll);
      $ll =~ s/\s//g;
      $len += length($ll);
    }
  }
      if ($len) {
        $N++;
        $sum += $len;
        push(@lens, $len);
      }

  my $ave = int($sum/$N);
  @lens = sort {$b <=> $a} @lens;
  $ref_len = $sum if (not defined($ref_len));

  my $longest = $lens[1];
  my $sum1 = 0;
  foreach $len (@lens) {
    $sum1 += $len;
    if (not $N50) { $N50 = $len if ($sum1 >= $ref_len * 0.5); }
    if (not $N80) { $N80 = $len if ($sum1 >= $ref_len * 0.8); }
    last if ($N50 and $N80);
  }
  print "#Contig_stat\tNumber\n";
  print "Number_contig\t$N\n";
  print "Total_length\t$sum\n";
  print "Ave_length\t$ave\n";
  print "Longest\t$longest\n";
  print "N50\t$N50\n";
  print "N80\t$N80\n";
}

close(TMP);

sub usage {
<<EOD;
Calculate N50 and several parameters

$script_name -i sequence_id_file -s original_fastq_sequence file -o output

  -i input fasta file 
  -r reference fasta file, to calculate the total genome length
     the program will use total contig length as total genome length, if -r and -L parameter are missing
     use only -r or only -L parameter
  -L total genome length
     the program will use total contig length as total genome length, if -r and -r parameter are missing
     use only -r or only -L parameter
EOD
}

