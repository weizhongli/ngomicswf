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

# Given the ORF KO annotation (by ann_ORF_taxon_func.pl)
# This script read in Brite file (e.g. ko00001.keg  ko00002.keg  ko01000.keg  ko02000.keg),
# and calculate abundance of KO, ko, class, modules etc
# these .keg file are from KEGG database /kegg/kegg-brite/ko/

use Getopt::Std;
getopts("i:k:a:o:e:d:s:t:r:",\%opts);
die usage() unless ($opts{i} and $opts{k} and $opts{o});

my $ORF_ann_file = $opts{i}; #### blast alignment file in m8 format
my $keg_file     = $opts{k}; #### e.g. ko00001.keg  ko00002.keg  ko01000.keg  ko02000.keg
my $output       = $opts{o}; #### output prefix kegg abundance file
my $ORF_depth    = $opts{d}; #### ORF depth
my $ref_KOs      = $opts{r}; #### single-copy house keeping gene

my ($i, $j, $k, $ll, $cmd);

#### read in ORF depth
my %ORF_depth = ();
my $ORF_depth_flag = 0;
if (defined($ORF_depth) and (-e $ORF_depth)) {
  open(TMP, $ORF_depth) || die "can not open $ORF_depth";
  while($ll=<TMP>){
    next if ($ll =~ /^#/);
    chop($ll);
    #### ID      depth
    my @lls = split(/\t/,$ll);     
    $ORF_depth{ $lls[0] } = $lls[1];
  }
  close(TMP);
  $ORF_depth_flag = 1;
}

my $single_copy_ref_flag = 0;
my %single_copy_KOs = ();
if (defined($ref_KOs) and (-e $ref_KOs)) {
  open(TMP, $ref_KOs) || die "can not open $ref_KOs";
  while($ll=<TMP>){ 
    if ($ll =~ /^\w\s+(K\d+)\s/) {
      $single_copy_KOs{$1} = 1;
    }
  }
  close(TMP);
  $single_copy_ref_flag = scalar keys %single_copy_KOs;
}


my $kegg_info = <<EOD;
==> ko00001.keg <==
+D	KO
#<h2><img src="/Fig/bget/kegg3.gif" align="middle" border=0> &nbsp; KEGG Orthology (KO)</h2>
#<!---
#ENTRY       ko00001
#NAME        KO
#DEFINITION  KEGG Orthology (KO)
#--->
!
A<b>Metabolism</b>
B
B  <b>Overview</b>
C    01200 Carbon metabolism [PATH:ko01200]
D      K00844  HK; hexokinase [EC:2.7.1.1]
D      K12407  GCK; glucokinase [EC:2.7.1.2]
D      K00845  glk; glucokinase [EC:2.7.1.2]

==> ko00002.keg <==
+E	Module
#<h2><img src="/Fig/bget/kegg3.gif" align="middle" border=0>&nbsp; KEGG Modules</h2>
#<!---
#ENTRY       ko00002
#NAME        Module
#DEFINITION  KEGG modules
#--->
!
A<b>Pathway module</b>
B
B  <b>Energy metabolism</b>
C    Carbon fixation
D      M00165  Reductive pentose phosphate cycle (Calvin cycle) [PATH:map01200 map00710]
E        K00855  PRK, prkB; phosphoribulokinase [EC:2.7.1.19]
E        K01601  rbcL; ribulose-bisphosphate carboxylase large chain [EC:4.1.1.39]

==> ko01000.keg <==
+E	Enzyme
#<h2><img src="/Fig/bget/kegg3.gif" align="middle" border=0>&nbsp; Enzymes</h2>
#<!---
#ENTRY       ko01000
#NAME        Enzyme
#DEFINITION  Enzymes
#--->
!
A<b>1. Oxidoreductases</b>
B  1.1  Acting on the CH-OH group of donors
C    1.1.1  With NAD+ or NADP+ as acceptor
D      1.1.1.1  alcohol dehydrogenase
E        K00001  E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]
E        K00121  frmA, ADH5, adhC; S-(hydroxymethyl)glutathione dehydrogenase / alcohol dehydrogenase [EC:1.1.1.284 1.1.1.1]

==> ko02000.keg <==
+E
#<h2><img src="/Fig/bget/kegg3.gif" align="middle" border=0>&nbsp; Transporters</h2>
#<!---
#ENTRY       ko02000
#NAME        Transporter
#DEFINITION  Transporters
#--->
!
A<b>ABC Transporters, Eukaryotic Type</b>
B  ABCA (ABC1) subfamily
C    ABCA1, 2, 3, 4, 7, 12, 13 subgroups
D      K05641  ABCA1; ATP-binding cassette, subfamily A (ABC1), member 1
D      K05642  ABCA2; ATP-binding cassette, subfamily A (ABC1), member 2
EOD


open(TMP, $keg_file) || die "can not open $keg_file";
my $keg_name = "";
while($ll=<TMP>){
  last if ($ll =~ /^!/);   #### read first block
  if ($ll =~ /^#NAME\s+(\w+)/) {
    $keg_name = $1;
  }
}

my %current_ABCDE_des = ();
my %level_is_ko = ();
my %KO_des = ();
my %ko_des = (); #### ko, MO or BR elements
my %cluster_member_KOs = ();
my @full_KO_link = ();
my %global_levels = ();
my %global_order = ();
my %non_KO_link = ();

my $i00 = 0;
while($ll=<TMP>){
  last if ($ll =~ /^!/);   #### finish last block
  next if ($ll =~ /^#/);
  chop($ll);
  if (($ll =~ /<b>/) and ($ll =~ /<\/b>/)) {
    $ll =~ s/<b>//;
    $ll =~ s/<\/b>//;
  }
  my $level = substr($ll,0,1);
  if ($level eq 'A') {
    %current_ABCDE_des = ();
  }
  elsif ($level eq 'B') { delete( $current_ABCDE_des{'B'}); delete( $current_ABCDE_des{'C'}); delete( $current_ABCDE_des{'D'}); }
  elsif ($level eq 'C') {                                   delete( $current_ABCDE_des{'C'}); delete( $current_ABCDE_des{'D'}); }
  elsif ($level eq 'D') {                                                                     delete( $current_ABCDE_des{'D'}); }

  my $txt   = substr($ll,1); $txt =~ s/^\s+//;
  next unless ($txt =~ /\w/);

  if ($txt =~ /^K\d+/) { #### this is KO level
    my ($KO, $des) = split(/\s+/, $txt, 2);
    $KO_des{$KO} = $des;
    my @p = qw/A B C/;
       @p = qw/A B C D/ if ($level eq "E"); #### enzyme .keg

    foreach $i (@p) {
      next unless (defined($current_ABCDE_des{$i}));
      $j = "$i $current_ABCDE_des{$i}";
      if (not defined($cluster_member_KOs{$j})) {
        $cluster_member_KOs{$j} = [];
      }
      push(@{ $cluster_member_KOs{$j} }, $KO);
    }
    push(@full_KO_link, [$KO, %current_ABCDE_des]);    
  }
  else {
    my $t_ko = "";
    if    ($txt =~ /^(M\d+)/) {
      $t_ko = $1;
      $txt =~ s/^M\d+\s+//;
      $ko_des{$t_ko} = $txt;
      $level_is_ko{$level} = 1;
    }
    elsif ($txt =~ /\[BR:(ko\d+)\]/) {
      $t_ko = $1;
      $txt =~ s/\s+\[BR:ko\d+\]//;
      $ko_des{$t_ko} = $txt;
      $level_is_ko{$level} = 1;
    }
    elsif ($txt =~ /\[PATH:(ko\d+)\]/) {
      $t_ko = $1;
      $txt =~ s/\s+\[PATH:ko\d+\]//;
      $txt =~ s/^\d+\s+//;
      $ko_des{$t_ko} = $txt;
      $level_is_ko{$level} = 1;
    }
    if ($t_ko) {
      $non_KO_link{"$level $t_ko"} = [%current_ABCDE_des];
    }
    elsif ($level ne 'A') {
      $non_KO_link{"$level $txt"} = [%current_ABCDE_des];
    }

    $current_ABCDE_des{$level} = $txt;
    $current_ABCDE_des{$level} = $t_ko if ($t_ko);
    if ($t_ko) { $global_order{"$level $t_ko"} = $i00; }
    else       { $global_order{"$level $txt"} = $i00; }
    $i00++;
    $global_levels{$level} = 1;
  }
}
close(TMP);

my $num_KO = scalar keys %KO_des;
my $num_ko = scalar keys %ko_des;

my %KO_abs = ();
my %KO_abs_adj = ();
my $sum_abs = 0;
my $sum_abs_adj = 0;
my $sum_ref_abs = 0;
my $sum_ref_abs_adj = 0;
open(TMP, $ORF_ann_file) || die "can not open $ORF_ann_file";
#Species_taxid\tSpecies\tGenome_taxid\tGenome\tScaffold\tORF\tStart\tEnd\tFrame\tIden%\tFamily\tDescription
while($ll=<TMP>) {
  next if ($ll =~ /^#/);
  chop($ll);
  my @lls = split(/\t/, $ll);
  my $fr = $lls[10]; #### fraction of coverage on reference protein
     $fr = 1.0 unless (($fr =~ /^\d/) and ($fr > 0.0));
  my $KO = $lls[11]; #### if hit KO
  next unless ($KO =~ /^K\d+/);

  my $ORF = $lls[5];
  my $depth = 1;
  if ($ORF_depth_flag and $ORF_depth{$ORF}) {
    $depth =  $ORF_depth{$ORF};
  }
  if ($single_copy_KOs{$KO}) {
    $sum_ref_abs += $fr;
    $sum_ref_abs_adj += $fr*$depth;
  }
  next unless ( $KO_des{$KO} ); #### unless this KO is in scope

  $KO_abs{$KO} += $fr;
  $KO_abs_adj{$KO} += $fr * $depth;
  $sum_abs += $fr;
  $sum_abs_adj += $fr*$depth;

} 
close(TMP);


if ($single_copy_ref_flag) {
  $sum_ref_abs /= $single_copy_ref_flag;
  $sum_ref_abs_adj /= $single_copy_ref_flag;
}
else {
  $sum_ref_abs = 1;
  $sum_ref_abs_adj = 1;
}


#### output KO
open(OUT, "> $output") || die "can not write to $output";
print OUT "#KO\tDepth\tDepth_adj\tAbundance\tAbundance_adj\tR.depth\tR.depth_adj\tDescription\n";
my @found_KOs = sort keys %KO_abs;
foreach $KO (@found_KOs) {
  next unless ( $KO_des{$KO} ); #### unless this KO is in scope
  my $abs       = $KO_abs{$KO};
  my $abs_adj   = $KO_abs_adj{$KO};
  my $r_abs     = float_e6( $abs    /$sum_abs );
  my $r_abs_adj = float_e6( $abs_adj/$sum_abs_adj );
  my $abs_2     = 0;
  my $abs_adj_2 = 0;
  if ($sum_ref_abs > 0) {
    $abs_2     = float_e3( $abs    /$sum_ref_abs );
    $abs_adj_2 = float_e6( $abs_adj/$sum_ref_abs_adj );
  }

  print OUT "$KO\t$abs\t$abs_adj\t$r_abs\t$r_abs_adj\t$abs_2\t$abs_adj_2\t$KO_des{$KO}\n";
}
close(OUT);
my $g_sum_abs = $sum_abs;
my $g_sum_abs_adj = $sum_abs_adj;

my @p = qw/D C B A/;
foreach $i (@p) {
  my @clusters = sort grep {/^$i/} keys %cluster_member_KOs;
  @clusters  = sort { $global_order{$a} <=>  $global_order{$b} } @clusters;
  next unless (@clusters);

  my $output_ann = "$output-$i";
  open(OUT, "> $output_ann") || die "can not write to $output_ann";

  if    ($i eq "B") { print OUT "#A\t"; }
  elsif ($i eq "C") { print OUT "#A\tB\t"; }
  elsif ($i eq "D") { print OUT "#A\tB\tC\t"; }
  elsif ($i eq "E") { print OUT "#A\tB\tC\tD\t"; }

  if ($level_is_ko{$i}) { print OUT     "ko\tDescription\tMember_KOs\tFound_KOs\tCoverage\tDepth\tDepth_adj\tR.depth\tR.depth_adj\n"; }
  else {                  print OUT             "Cluster\tMember_KOs\tFound_KOs\tCoverage\tDepth\tDepth_adj\tR.depth\tR.depth_adj\n"; }

  foreach $j (@clusters) {
    if ( $non_KO_link{$j} ) {
      my %link = @{ $non_KO_link{$j} };
      my @link = sort keys %link;
      foreach $k (@link) {
        print OUT "$link{$k}\t";
      }
    } 

    my @member_KOs = @{ $cluster_member_KOs{$j} };
    my $abs = 0;
    my $abs_adj = 0;
    my $no = 0;
    my $abs_2 = 0;
    my $abs_adj_2 = 0;
    foreach $KO (@member_KOs) {
      $no ++                       if ($KO_abs{$KO});
      $abs     += $KO_abs{$KO}     if ($KO_abs{$KO});
      $abs_adj += $KO_abs_adj{$KO} if ($KO_abs_adj{$KO});
    }
    if ( @member_KOs ) {
      $abs /= $#member_KOs+1;
      $abs_adj /= $#member_KOs+1;
      $abs_2     = 0;
      $abs_adj_2 = 0;
      if ($sum_ref_abs > 0) {
        $abs_2     = float_e3( $abs    /$sum_ref_abs );
        $abs_adj_2 = float_e6( $abs_adj/$sum_ref_abs_adj );
      }
    }

    my $cov = float_e3( $no/($#member_KOs+1) ) ;
    if ($level_is_ko{$i} ) {
      print OUT substr($j, 2), "\t", $ko_des{substr($j, 2)}, "\t", $#member_KOs+1, "\t$no\t$cov\t$abs\t$abs_adj\t$abs_2\t$abs_adj_2\n";
    }
    else {
      print OUT substr($j, 2), "\t",                               $#member_KOs+1, "\t$no\t$cov\t$abs\t$abs_adj\t$abs_2\t$abs_adj_2\n";
    }
  }
  close(OUT);
}

open(OUT, "> $output-full-des") || die "can not write to $output-full-des";
my @g_levels = sort keys %global_levels;
print OUT "#ID";
foreach $j (@g_levels) {
  print OUT "\t$j";
  print OUT "\tko" if ($level_is_ko{$j});
}
print OUT "\tKO\tDescription\tDepth\tDepth_adj\tAbundance\tAbundance_adj\tR.depth\tR.depth_adj\n";

$i00 = 1000001;
for $i (@full_KO_link) {
  my ($KO, %ABCDE_des) = @{ $i };
  my @ABC = sort keys %ABCDE_des;

  print OUT "$i00"; $i00++;
  foreach $j (@g_levels) {
    my $des = "";
    $des = $ABCDE_des{$j} if defined($ABCDE_des{$j});
    print OUT "\t$des";
    if ($level_is_ko{$j}) {
      print OUT (defined($ko_des{$des})) ? "\t$ko_des{$des}" : "\t$des";
    }
  }
  print OUT "\t$KO\t$KO_des{$KO}";
  my ($abs, $abs_adj, $r_abs, $r_abs_adj, $abs_2, $abs_adj_2) = qw/0 0 0 0 0 0/;
  if (defined ($KO_abs{$KO})) {
    $abs       = $KO_abs{$KO};
    $abs_adj   = $KO_abs_adj{$KO};
    $r_abs     = float_e6( $abs    /$g_sum_abs );
    $r_abs_adj = float_e6( $abs_adj/$g_sum_abs_adj );
    $abs_2     = 0;
    $abs_adj_2 = 0;
    if ($sum_ref_abs > 0) {
      $abs_2     = float_e3( $abs    /$sum_ref_abs );
      $abs_adj_2 = float_e6( $abs_adj/$sum_ref_abs_adj );
    }
  }
  print OUT "\t$abs\t$abs_adj\t$r_abs\t$r_abs_adj\t$abs_2\t$abs_adj_2\n";
}
close(OUT);

sub float_e3 {
  my $f = shift;
  return int($f*1000) / 1000;
}

sub float_e6 {
  my $f = shift;
  return int($f*1000000) / 1000000;
}

sub usage {
<<EOD;
$script_name -i ORF_annotation_file -k .keg_file -o output

  options:
    -i ORF_annotation_file, generated by ann_ORF_taxon_func.pl
    -k .keg file, e.g. ko00001.keg  ko00002.keg  ko01000.keg  ko02000.keg
    -o output prefix,
    -d ORF depth file
    -r .keg file, contains a list of KOs of single copy house-keeping gene 
       as reference to calculate relative abundance. 
       e.g. there are 55 KOs M00178  Ribosome, bacteria [PATH:map03010] [BR:ko03011]

EOD
}
