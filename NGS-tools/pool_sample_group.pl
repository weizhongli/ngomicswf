#!/usr/bin/perl
#
use Getopt::Std;
getopts("s:S:o:f:j:g:G:",\%opts);

die usage() unless ($opts{s} and $opts{j} and $opts{f});

my $group_col         = $opts{g};
   $group_col         = 1 unless ($group_col);
my $group_G           = $opts{G}; 
my $sample_in         = $opts{s};
my $file_list         = $opts{f};
my @file_list         = split(/,/, $file_list) if ($file_list); 
my $job               = $opts{j};


my ($i, $j, $k, $cmd);

######## parse NGS_samples
my @NGS_samples = ();
my %sample_2_group = ();
my %groups = ();
open(TMP, $sample_in) || die "can not open $sample_in";
while($ll=<TMP>){
  next if ($ll =~ /^#/);
  next unless ($ll =~ /^\w/); chop($ll);
  my ($id, @data) = split(/\s+/,$ll);
  push(@NGS_samples, $id);
  my $group = $data[ $group_col - 1];
  if ($group_G) {
    $group = $group_G;
  }
  $group =~ s/\s/_/g;
  $sample_2_group{$id} = $group;
  $groups{$group} = 1;
}
close(TMP);

foreach $i (keys %groups) {
  $cmd = `mkdir -p $i/$job`;
}

foreach $i (@file_list) {
  foreach $j (@NGS_samples) {
    my $source = "$j/$job/$i";
    my $target = "$sample_2_group{$j}/$job/$i";

    if (-e $source) {
      print STDERR "cat $source >> $target\n";
      $cmd = `cat $source >> $target`;
    }
    else {
      print STDERR "Warning, $source missing\n";
    }
  }
}

sub usage {
<<EOD;
Given a samle file, this script pools samples into groups by cat individual files
from each samples

    $0 -s sample_file -j job -f list_files
    -s sample data file, first field is sample_id, 2nd column can be group id
    -g column (0-based) for group id, default 1, i.e. 2nd column
    -G manual group id for all samples, if input, this will overwrite group id in sample file
       this will pool all samples together
    -j job name
    -f list of files (delimited by ,) to pool (cat) 
EOD
}

