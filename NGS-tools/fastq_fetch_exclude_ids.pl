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
my $id_str     = $opts{a};
my $output     = $opts{o};

my %seqid1 = ();
my $seqid;

open(TMP, $seqid_file) || die;
while($ll = <TMP>) {
  chop($ll);
  $ll =~ s/\s.+$//;
  $ll =~ s/^>//;
  $seqid1{$ll} = 1;
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
while($ll = <TMP>) {
  if ($ll =~ /^\@/) {
    if ($id_str) {
      if ($ll =~ /$id_str\|(\d+)/) {
        $seqid = $1;
      }
    }
    else {
      $seqid = substr($ll,1);
      chop($seqid);
      $seqid =~ s/\s.+$//;
    }
    $flag = ( $seqid1{$seqid} ) ? 0 : 1;
    print $fh $ll  if ($flag);
    $a=<TMP>; print $fh $a if ($flag);
    $a=<TMP>; print $fh $a if ($flag);
    $a=<TMP>; print $fh $a if ($flag);
  }
}
close(TMP);
close(OUT) if ($output);

sub usage {
<<EOD;
fetch subset of sequences by excluding IDs

$script_name -i sequence_id_file -s original_fastq_sequence file -o output

  -i input file of sequence IDs,  one ID per line, ID must be the first column, - for stdin
  -s original fastq file
  -o output, default STDOUT
  -a ID string, optional, if it is supplied, it will match sequence id with ID_string|ID
     for example, if with -a taxid, the input file can be a list of taxids and the sequence name
     need to be have taxid|12345 pattern to match
EOD
}

