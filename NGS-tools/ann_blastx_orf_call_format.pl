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
mHE-SRS012902|scaffold|86.16    gnl|CDD|226997  47.62   42      17      2       164     201     210     250     5e-04   37.6
mHE-SRS012902|scaffold|109.23   gnl|CDD|225183  47.46   236     122     1       1       236     475     708     1e-92    284
mHE-SRS012902|scaffold|109.23   gnl|CDD|224055  44.35   239     130     2       1       239     332     567     2e-84    259
mHE-SRS012902|scaffold|109.23   gnl|CDD|227321  39.50   238     140     3       1       238     324     557     9e-69    218
EOD

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/ann_local.pl";

use Getopt::Std;
getopts("i:a:o:e:d:s:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $input        = $opts{i};
my $output       = $opts{o};

my ($i, $j, $k, $ll, $cmd);

open(TMP, $input) || die "Can not open $input";
open(OUT, "> $output") || die "Can not write to $output";

while($ll=<TMP>){
  chop($ll);
  my @lls = split(/\t/, $ll);

  my $new_id = "$lls[0].$lls[6].$lls[7]";
  $locb = ($lls[6] < $lls[7]) ? $lls[6] : $lls[7];
  $loce = ($lls[6] < $lls[7]) ? $lls[6] : $lls[7];
  $len  = int( ($loce - $locb + 1) / 3);

  print OUT "$new_id\t$lls[1]\t$lls[2]\t$lls[3]\t$lls[4]\t$lls[5]\t1\t$len\t$lls[8]\t$lls[9]\t$lls[10]\t$lls[11]\n";
}
close(TMP);
close(OUT);

