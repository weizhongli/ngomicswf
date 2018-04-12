#!/usr/bin/perl -w
## ==============================================================================
## Automated annotation tools
##
## program written by
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
##                                      http://weizhong-lab.ucsd.edu
## ==============================================================================

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/ann_local.pl";

use Getopt::Std;
getopts("i:a:o:e:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $results_dir  = $opts{i};
my $parse_output = $opts{o};
my $orf_abs      = $opts{a}; #### abundance file for ORFs
my $e_cutoff     = $opts{e}; $e_cutoff = 0.001 unless (defined($e_cutoff));
my $overlap_cutoff = 0.5;
die "hmmer3 output results dir $results_dir not found" unless (-e $results_dir);
my ($i, $j, $k, $ll, $cmd);

my %orf_abs = ();
my $abs_flag = 0;
if ($orf_abs) {
  open(TMP, $orf_abs) || die "Can not open $orf_abs";
  while($ll=<TMP>){
    next if ($ll =~ /^#/);
    chop($ll);
    my @lls = split(/\t/,$ll);     #### ID	length	cov	relative abundance
    $orf_abs{$lls[0]} = $lls[2]; 
  }
  close(TMP);
  $abs_flag = 1;
}

my $hmm_acc_2_name = ();
my $hmm_acc_2_des = ();
my $hmm_acc_2_len = ();

my %hmm_orf_adjusted = ();
my %hmm_cov_adjusted = ();

my @os;
if (-d $results_dir) { ##### a directory of output files
  @os = LL_get_active_ids($results_dir);
  @os = sort @os;
  foreach $i (@os) { parse_it($i);}
}
else { ##### a single output file
  parse_it();
}


my $total_cov_adjusted=0;
foreach $hmm_acc (keys %hmm_orf_adjusted) {
  $total_cov_adjusted += $hmm_cov_adjusted{$hmm_acc};
}

open(OUT, "> $parse_output") || die "Can not write to $parse_output";
print OUT "#", join("\t", qw/Accession No_orfs Coverage Abundance Name Description/), "\n";
my @hmm_accs = sort keys %hmm_orf_adjusted;
foreach $hmm_acc (@hmm_accs) {
  my $countd = int($hmm_orf_adjusted{$hmm_acc});
  my $covd   = int($hmm_cov_adjusted{$hmm_acc} * 1000000) / 1000000;
  my $absd   = int($covd / $total_cov_adjusted * 1000000) / 1000000;
  print OUT "$hmm_acc\t$countd\t$covd\t$absd\t$hmm_acc_2_name{$hmm_acc}\t$hmm_acc_2_des{$hmm_acc}\n";
}
close(OUT);

#### using output from hmmscan option --domtblout  -E 0.001 --notextw --cut_tc --noali 
sub parse_it {
  my $id = shift;
  my ($i, $j, $k, $ll);
  my $tout = (defined($id)) ? "$results_dir/$id" : $results_dir;

  my $last_seq_id = "";
  my @last_e;
  my @last_b;
  open(TMP, $tout) || next;
  while($ll=<TMP>) {
    chop($ll);
    next if ($ll =~ /^#/);
    my @lls = split(/\s+/, $ll, 23); 
    my $output_looks_like = <<EOD; 
ABC2_membrane        PF01061.19   210 mHE-SRS012902|scaffold|2.77 -            278   2.8e-21   75.6  24.6   1   2   1.9e-25   2.8e-21   75.6  17.1     5   205    27   234    23   239 0.84 ABC-2 type transporter
ABC2_membrane        PF01061.19   210 mHE-SRS012902|scaffold|2.77 -            278   2.8e-21   75.6  24.6   2   2     0.091   1.4e+03   -1.8   1.1   122   149   248   264   234   274 0.46 ABC-2 type transporter
EOD
    my $hmm_name = $lls[0];
    my $hmm_acc  = $lls[1];
    my $hmm_len  = $lls[2];
    my $seq_id   = $lls[3];
    my $seq_len  = $lls[5];
    my $e_value  = $lls[6];
    my $hmm_b    = $lls[15];
    my $hmm_e    = $lls[16];
    my $seq_b    = $lls[17];
    my $seq_e    = $lls[18];
    my $hmm_des  = $lls[22];
    my $hmm_aln_ln = $hmm_e - $hmm_b + 1; 
    next unless ($e_value <= $e_cutoff);

    if (not defined($hmm_acc_2_len{$hmm_acc})) {
      $hmm_acc_2_name{$hmm_acc} = $hmm_name;
      $hmm_acc_2_len{$hmm_acc}  = $hmm_len;
      $hmm_acc_2_des{$hmm_acc}  = $hmm_des;
    }
    if ($seq_id ne $last_seq_id) { @last_e = (); @last_b = (); } #### start a new sequence
    my $overlap_with_before = 0;
    for ($j=0; $j<@last_b; $j++) {
      my $lb=$last_b[$j];
      my $le=$last_e[$j];
      if ( overlap1($lb,$le,$seq_b,$seq_e) > ($seq_e-$seq_b+1) * $overlap_cutoff) {
       $overlap_with_before=1; last;
      }
    }
    #print $overlap_with_before ? "$ll***\n": "$ll\n";

    next if ($overlap_with_before);
    my $cov1 = $hmm_aln_ln / $hmm_len;
    my $weight = ($abs_flag and defined($orf_abs{$seq_id})) ?  $orf_abs{$seq_id} : 1;
    $hmm_orf_adjusted{$hmm_acc} +=         $weight;
    $hmm_cov_adjusted{$hmm_acc} += $cov1 * $weight;
    push(@last_e, $seq_e);
    push(@last_b, $seq_b);
    $last_seq_id = $seq_id;
  }
  close(TMP);
}


sub overlap1 {
  my ($b1, $e1, $b2, $e2) = @_;
  return 0 if ($e2 < $b1);
  return 0 if ($b2 > $e1);
  return ( ($e1<$e2)? $e1:$e2 )-( ($b1>$b2)? $b1:$b2);
}

