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
getopts("i:s:o:a:d:",\%opts);
die usage() unless ($opts{i} and $opts{s});

my $seqid_file = $opts{i};
my $seq_file   = $opts{s};
my $idx_file   = $opts{d};
   $idx_file   = "$seq_file.fai" unless ($idx_file);
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

#.fai template
#acc|NZ_UGSS01000002.1|taxid|757|sptid|757|getid|N  2596673 120 80  81
#acc|NZ_UGSS01000001.1|taxid|757|sptid|757|getid|N  16175   2629372 80  81
#acc|NZ_DS480351.1|taxid|428125|sptid|1535|getid|N  3450    2645873 80  81
#acc|NZ_DS480350.1|taxid|428125|sptid|1535|getid|N  523705  2649490 80  81
#...
#acc|NZ_KB911070.1|taxid|1121334|sptid|1549|getid|N 205920  7891764 80  81

my @seq_idx = ();
open(TMP, $idx_file) || die "can not open $idx_file";
while($ll=<TMP>){
  chop($ll);
  my($seqid, $len, $pos, $ll1, $ll2) = split(/\t/,$ll);
  next unless  ($seqid1{$seqid});
  my $num_lines = int($len / $ll1);
  $num_lines ++ if ($len % $ll1);
  my $new_line_letters = $num_lines * ($ll2 - $ll1);
  my $total_c = $len + $new_line_letters;

  push(@seq_idx, [$seqid, $pos, $total_c]);
}
close(TMP);

@seq_idx = sort { $a->[1] <=> $b->[1] } @seq_idx;

my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $OUT";
  $fh = "OUT";
}
my $flag = 0;

open(TMP, $seq_file) || "can not open $seq_file";

my $data;
foreach $i (@seq_idx) {
  my ($seq_id, $pos, $total_c) = @{$i};
  seek(TMP, $pos, 0);
  $data_c = read(TMP, $data, $total_c);
  if ($data_c) {
    print $fh ">$seq_id\n$data";
  }
}

close(TMP);
close(OUT) if ($output);

sub usage {
<<EOD;
fetch subset of sequences by IDs
This works for HUGE databases and the database has to be pre-indexed by
samtools faidx

for small database, use fasta_fetch_by_ids.pl

$script_name -i sequence_id_file -s original_fastq_sequence file -o output

  -i input file of sequence IDs,  one ID per line, ID must be the first column, - for stdin
  -s original fasta file, must be a physical unzipped file, can not be stdin
  -d original fasta index file, default is to add .fai to the original fasta file
  -o output, default STDOUT
EOD
}

