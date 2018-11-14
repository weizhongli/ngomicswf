#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

## concatenate annotation tables

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);

use Getopt::Std;

getopts("i:o:d:f:c:t:s:r:n:l:s:a:v:h:p:j:",\%opts);
die usage() unless (defined($opts{i}) and $opts{s} and $opts{o} and $opts{f});

my $sample_file         = $opts{s}; #### NGS-samples
my $file                = $opts{f}; #### file name eg blast-kegg-parse-bin/pathway.txt
my $col_str             = $opts{i}; #### 1'based columns for data e.g. 1,2,3,4
my $ann_str             = $opts{a}; #### 1'based columns for feature dictionary e.g. 1,2,3
my $output              = $opts{o};
my $output_dic          = "";
my $output_dic_tmp      = "";

if ($ann_str) {
  if ($output =~ /\.(\w+)$/) { $output_dic = $output; $output_dic =~ s/\.(\w+)$/-dic.$1/; }
  else                       { $output_dic = "$output-dic"; }
  $output_dic_tmp = "$output_dic.$$";
}

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

if (1) {
  my $sample_1st = $samples[0];
  $cmd = `grep "^#"  $sample_1st/$file | cut -f $col_str | sed "s/^/sample_id\\t/" > $output`;
  $cmd = `grep "^#"  $sample_1st/$file | cut -f $ann_str > $output_dic` if ($ann_str);

  foreach $sample (@samples) {
    $cmd = `grep -v "^#" $sample/$file | cut -f $col_str | sed "s/^/$sample\\t/" >> $output`;
    $cmd = `grep -v "^#" $sample/$file | cut -f $ann_str >> $output_dic_tmp` if ($ann_str);
  }
  $cmd = `sort  $output_dic_tmp | uniq >> $output_dic` if ($ann_str);
  $cmd = `rm -f $output_dic_tmp` if ($ann_str);

}


sub usage {
<<EOD
this script concatenate txt files for mulitple samples into larger spreadsheet
options:
    -s sample_file, e.g. NGS-smaples
    -o output 
    -f file  em_fast/Bacteria-toprank-abundance.txt
    -i 1-based columns for content, e.g. 1.2,3,4, required
    -a 1-based columns for annotation description, e.g. 1,2
EOD
}
########## END usage


