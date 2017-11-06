#!/usr/bin/perl -w
## ==============================================================================
## Automated annotation tools
##
## program written by
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
##                                      http://weizhong-lab.ucsd.edu
## ==============================================================================

## given blastx output querying scaffold against protein database in tabluar format 
## change the format, add cols 6,7 to column 0
## so that column 0 be new ORF id
    my $output_looks_like = <<EOD; 
#query                          subject         %       alnln   mis     gap     q_b     q_e     s_b     s_e     expect  bits
#0                              1               2       3       4       5       6       7       8       9       10      11
==> input <==
sample01|scaffold|0	ti|40521|sptax|40521|NP_463477.1	86.61	1792	213	3	162960	157636	1	1782	0.0	 2682
sample01|scaffold|0	ti|40521|sptax|40521|NP_463464.1	94.99	499	22	1	170594	169098	1	496	0.0	  981
==> output  <==
sample01|scaffold|0.162960.157636	ti|40521|sptax|40521|NP_463477.1	86.61	1792	213	3	1	1775	1	1782	0.0	 2682
sample01|scaffold|0.170594.169098	ti|40521|sptax|40521|NP_463464.1	94.99	499	22	1	1	499	1	496	0.0	  981
EOD

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/ann_local.pl";

use Getopt::Std;
getopts("i:a:o:e:d:s:c:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $input        = $opts{i};
my $output       = $opts{o};
my $overlap_cutoff = $opts{c}; 
   $overlap_cutoff = 0.5 unless defined($overlap_cutoff);

my ($i, $j, $k, $ll, $cmd);

open(TMP, $input) || die "Can not open $input";
open(OUT, "> $output") || die "Can not write to $output";

my $last_seq_id = "";
my @last_e;
my @last_b;

while($ll=<TMP>){
  chop($ll);
  my @lls = split(/\t/, $ll);
  my $seq_id = $lls[0];

  my $new_id = "$lls[0].$lls[6].$lls[7]";
  $locb = ($lls[6] < $lls[7]) ? $lls[6] : $lls[7];
  $loce = ($lls[6] > $lls[7]) ? $lls[6] : $lls[7];
  $len  = int( ($loce - $locb + 1) / 3);

  if ($seq_id ne $last_seq_id) { @last_e = (); @last_b = (); } #### start a new sequence

  my $overlap_with_before = 0;
  for ($j=0; $j<@last_b; $j++) {
    my $lb=$last_b[$j];
    my $le=$last_e[$j];
    if ( overlap1($lb,$le,$locb,$loce) > ($loce-$locb+1) * $overlap_cutoff) {
     $overlap_with_before=1; last;
    }
  }
  next if ($overlap_with_before);
  print OUT "$new_id\t$lls[1]\t$lls[2]\t$lls[3]\t$lls[4]\t$lls[5]\t1\t$len\t$lls[8]\t$lls[9]\t$lls[10]\t$lls[11]\n";

  push(@last_e, $loce);
  push(@last_b, $locb);
  $last_seq_id = $seq_id;
}
close(TMP);
close(OUT);



sub overlap1 {
  my ($b1, $e1, $b2, $e2) = @_;
  return 0 if ($e2 < $b1);
  return 0 if ($b2 > $e1);
  return ( ($e1<$e2)? $e1:$e2 )-( ($b1>$b2)? $b1:$b2);
}


sub usage {
<<EOD;
given blastx output querying scaffold against protein database in tabluar format 
change the format, add cols 6,7 to column 0 so that column 0 be new ORF id.
It also remove redundant alignments.

$script_name -i input_blastx_alignment -o output_blastx_alignment

  options:
    -i blast alignment file in tab format
    -o output blastt alignment
    -c overlap cutoff, default 0.5 
       delete redundant alignments that overlap with previous alignments
       overlap is on query sequences (i.e. scaffolds) 

EOD
}

