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
getopts("i:a:o:e:d:n:",\%opts);
die usage() unless ($opts{i} and $opts{o} and $opts{d} and $opts{n});

my $family_output= $opts{i}; #### per COG family-based output
my $output       = $opts{o};
my $class_name   = $opts{n}; #### COG / KOG func.txt
my $family_class = $opts{d}; #### COG family to COG class mapping

my ($i, $j, $k, $ll, $cmd);

my %class_2_des = ();
my @classes = ();
open(TMP, $class_name) || die "Can not open $class_name";
while($ll=<TMP>){
  chop($ll);
  next unless ($ll =~ /^\s+\[\w+\]/);
  $ll =~ s/^\s+//;
  my ($id, $des) = split(/\s+/,$ll,2);
  $id = substr($id,1); chop($id);
  $class_2_des{$id} = $des;
  push(@classes, $id);
}
close(TMP);

my %class_2_family;
open(TMP, $family_class) || die "Can not open $family_class";
while($ll=<TMP>){
  chop($ll);
  next unless ($ll =~ /^\[\w+\]/);
  my ($id, $family, $des) = split(/\s+/, $ll, 3);
  $id = substr($id,1); chop($id);
  my @ids = split(/ */, $id);
  foreach $i (@ids) {
    $class_2_family{$i}{$family}=1;
  }
}
close(TMP);

my %family_cov;
open(TMP, $family_output) || die "Can not open $family_output";
#Accession      Name    No_orfs Coverage        Abundance       Description
#COG0002 ArgC    118     50.526327       0.000426        Acetylglutamate semialdehyde dehydrogenase
while($ll=<TMP>){
  chop($ll);
  next if ($ll =~ /^#/);
  my ($family, $name, $no_orf, $cov, $abs, $des) = split(/\t/,$ll);
  $family_cov{$family} = $abs;
}
close(TMP);


open(OUT, "> $output") || die "Can not write to $output";
print OUT "#", join("\t", qw/Class No_families Coverage Abundance Description/), "\n";
foreach $i (@classes) {
  my @t_families = keys %{$class_2_family{$i}};
  my $t_no       = scalar @t_families;
  next unless ($t_no>0);

  my $t_family_no = 0;
  my $sum_abs = 0;
  foreach $j (@t_families) {
    next unless ($family_cov{$j});
    $t_family_no++;
    $sum_abs += $family_cov{$j};
  }
  my $cov1 = int($t_family_no/$t_no * 1000000)/1000000;
  print OUT "$i\t$t_no\t$cov1\t$sum_abs\t$class_2_des{$i}\n";
}
close(OUT);

