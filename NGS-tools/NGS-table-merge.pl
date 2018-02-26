#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

## make sample ~ feature table. Feature can be taxonomy, functional etc
## row is feature, with possible annotation on the first few following columns
## other columns are samples

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);

use Getopt::Std;

getopts("i:o:d:f:c:t:s:r:n:l:s:a:v:h:p:",\%opts);
die usage() unless (defined($opts{i}) and $opts{s} and $opts{o} and $opts{f});

my $sample_file         = $opts{s}; #### NGS-samples
my $file                = $opts{f}; #### file name eg blast-kegg-parse-bin/pathway.txt
my $id_col              = $opts{i}; #### 0-based column for ID, default 0;
   $id_col              = 0  unless (defined $id_col);
my $ann_cols            = $opts{a}; #### 0-based columns for annotations delimited by ",", defaults none
my $value_col           = $opts{v}; #### 0-based value column, default -1
   $value_col           = -1 unless (defined $value_col);
my $output              = $opts{o};
my $cutoff              = $opts{c}; #### cutoff values
   $cutoff              = 0  unless (defined $cutoff);
my $sort_flag           = $opts{t}; $sort_flag = 1 unless (defined $sort_flag);
my $stat_flag           = $opts{p}; $stat_flag = 0 unless (defined $stat_flag);

my @ann_cols = ();if ($ann_cols) { @ann_cols = split(/,/, $ann_cols);}

my ($i, $j, $k, $ll, $cmd);

my @samples = ();
open(TMP, $sample_file) || die "can not open $sample_file\n";
while($ll=<TMP>){
  chop($ll);
  next if ($ll =~ /^#/);
  my @lls = split(/\s+/,$ll);
  push(@samples, $lls[0]);
}
close(TMP);
my $no_samples = $#samples+1;


if (1) {
  my ($id, $sample, $value);
  my %mat = ();
  my %ids = ();
  my %id_2_ann = ();
  my @ann_names;
  my $id_name;

  foreach $sample (@samples) {
    my $f1 = "$sample/$file";
    open(TMP, $f1) || die "can not open $f1";
    $ll=<TMP>; chop($ll);
    my @lls = split(/\t/,$ll);
    if (not @ann_names) {
      @ann_names = map {$lls[$_]} @ann_cols;
    }
    if (not defined($id_name))    {$id_name    = $lls[$id_col];}

    while($ll=<TMP>){
      chop($ll);
      @lls = split(/\t/,$ll);
      $id = $lls[$id_col];
      $value = $lls[$value_col];
      next unless ($value >= $cutoff);

      $ids{$id} = 1;
      if (not defined( $id_2_ann{$id})) {
        my @t1 = map {$lls[$_]} @ann_cols;
        $id_2_ann{$id} = [ map {$lls[$_]} @ann_cols ];
      }
      $mat{$id}{$sample} = $value;
    }
    close(TMP);
  }

  my @ids = sort keys %ids;
  my %id_2_stat=();
  foreach $id (@ids) {
    my @vv=();
    foreach $sample (@samples) {
      my $v1 = $mat{$id}{$sample};
         $v1 = 0 unless defined($v1);
      push(@vv, $v1);
    }
    $id_2_stat{$id} = [ max_min_median_ave_range(@vv) ] if ($no_samples>0);
  }

  $id_name = "ID" unless $id_name;

  open(TOUT, "> $output") || die "can not write to $output";
  print TOUT "$id_name";
  print TOUT "\t", join("\t", @ann_names) if (@ann_cols);
  print TOUT "\t", join("\t", @samples); 
  print TOUT "\tmax\tmin\tmed\tmean\trange\tlow_quartile\thigh_quartile" if ($no_samples>0 and $stat_flag);
  print TOUT "\n";

     @ids = sort {$id_2_stat{$b}->[0] <=> $id_2_stat{$a}->[0]} @ids if (($no_samples>0) and $sort_flag and $stat_flag);

  foreach $id (@ids) {
    print TOUT $id;
    print TOUT "\t", join("\t", @{$id_2_ann{$id}}) if (@ann_cols);
    foreach $sample (@samples) {
      my $v1 = $mat{$id}{$sample};
         $v1 = 0 unless defined($v1);
      print TOUT "\t$v1";
    }
    print TOUT "\t", join("\t", @{$id_2_stat{$id}}) if ($no_samples>0 and $stat_flag);
    print TOUT "\n";
  }
  close(TOUT);
}

sub max_min_median_ave_range {
  my @v1 = @_;
  my ($i, $j, $k);
  my $no = $#v1+1;

  @v1 = sort {$a <=> $b} @v1;
  
  my $min = $v1[0];
  my $max = $v1[0];
  my $sum = 0;
  my ($med, $low_q, $high_q);

  $j=0;
  foreach $i (@v1) {
    $min = $i if ($i<$min);
    $max = $i if ($i>$max);
    $sum += $i;
    $j++;
  }
  $sum /= $j if ($j>0);

  if ($no % 2 ==0 ) { $med = ($v1[$no/2] + $v1[$no/2-1])/2.0; }
  else              { $med = $v1[int($no/2)]; }

  $low_q = $v1[int($no/4)]; ### low quartile
  @v1 = reverse(@v1);
  $high_q = $v1[int($no/4)]; 

  my $range = $max-$min;
  return ($max, $min, $med, $sum, $range, $low_q, $high_q);
}
########## END max_min_median_ave_range


sub usage {
<<EOD
this script merge txt files for mulitple samples into larger spreadsheet
options:
    -i sample_file, e.g. NGS-smaples
    -o output 
    -f file  em_fast/Bacteria-toprank-abundance.txt
    -i 0-based column for ID, default 0;
    -a 0-based columns for annotations delimited by ",", defaults none
    -v 0-based value column, default -1
    -c cutoff values delimited by ","

EOD
}
########## END usage





