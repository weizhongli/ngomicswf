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

# Based on alignments between ORFs and reference DB. e.g. kegg
# generated by cd-hit-2d and blast, and cluster info of reference DB,
# and ORF-scaffold membership
# annotate ORFs taxon and function, satisfy that ORFs belong to the same 
# scaffold get same of consistent taxon
use Getopt::Std;
getopts("i:r:a:o:e:d:s:t:",\%opts);
die usage() unless ($opts{i} and $opts{r} and $opts{a} and $opts{o} and $opts{t} );

my $bl_file      = $opts{i}; #### blast alignment file in m8 format
my $clstr_ref    = $opts{r}; #### cluster info
my $ORF_file     = $opts{a}; #### old ORF
my $taxon_file   = $opts{t}; #### taxon file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl
my $output       = $opts{o}; #### output ORF file
my $output_ann   = "$output-ann.txt"; #### output annotation
my $output_tax   = "$output-tax.txt"; #### output tax 
my $output_log   = "$output-ann.log";
my $cutoff_e     = $opts{e}; 
   $cutoff_e     = 1e-6 unless defined($cutoff_e);

my ($i, $j, $k, $ll, $cmd);

my %taxon_info = ();
open(TMP, $taxon_file) || die "can not open $taxon_file";
while($ll=<TMP>) {
  chop($ll);
  my ($tid, $rank, @lls) = split(/\t/, $ll);
  next unless ($rank eq "toprank");
  $taxon_info{$tid} = [@lls];
#511145	toprank	
#0                                            	1        	2	3         	4	5                 	6	7                 	8	9                  	10	11          	12	13               	14    	15                                           	16
#Escherichia coli str. K-12 substr. MG1655	Bacteria	2	Proteobacteria	1224	Gammaproteobacteria	1236	Enterobacterales	91347	Enterobacteriaceae	543	Escherichia	561	Escherichia coli	562	Escherichia coli str. K-12 substr. MG1655	511145

}
close(TMP);

my %orf_2_scaffold = ();
my %orf_info = (); #### start, end
my %scaffold_member_orfs = ();
my %scaffold_orf_count = ();
open(TMP, $ORF_file) || die "can not open $ORF_file";
while($ll=<TMP>) {
  if ($ll =~ /^>(\S+)/) {
    my $orf_id = $1;
    $orf_info{$orf_id} = join("\t", ORF_info($ll)) ;
    if ($orf_id =~ /(scaffold\|\d+)/) {
      my $sid = $1;
      $orf_2_scaffold{$orf_id} = $sid;
      if (not defined( $scaffold_member_orfs{$sid} )) {
        $scaffold_member_orfs{$sid} = [];
      } 
      push(@{$scaffold_member_orfs{$sid}}, $orf_id); 
    }
  }
} 
foreach $sid (keys %scaffold_member_orfs) {
  $scaffold_orf_count{$sid} = scalar @{$scaffold_member_orfs{$sid}};
}
close(TMP);
my @all_sids = keys %scaffold_member_orfs;
   @all_sids = sort { $scaffold_orf_count{$b} <=> $scaffold_orf_count{$a} } @all_sids;

my %ref_ids = ();
my %orf_2_hit = ();
$last_ORF = "";
if (-d $bl_file) {
 open(TMP, "cat $bl_file/* |") || die "can not open $bl_file";
}
elsif (-e $bl_file) {
  open(TMP, $bl_file) || die "can not open $bl_file";
}

    my $output_looks_like = <<EOD; 
#query                          subject         %       alnln   mis     gap     q_b     q_e     s_b     s_e     expect  bits
#0                              1               2       3       4       5       6       7       8       9       10      11
mHE-SRS012902|scaffold|86.16    gnl|CDD|226997  47.62   42      17      2       164     201     210     250     5e-04   37.6
mHE-SRS012902|scaffold|109.23   gnl|CDD|225183  47.46   236     122     1       1       236     475     708     1e-92    284
mHE-SRS012902|scaffold|109.23   gnl|CDD|224055  44.35   239     130     2       1       239     332     567     2e-84    259
mHE-SRS012902|scaffold|109.23   gnl|CDD|227321  39.50   238     140     3       1       238     324     557     9e-69    218
EOD

while($ll=<TMP>) {
  #ser:SERP1011|ti|176279|KO||len|10203
  chop($ll);
  my @lls = split(/\t/, $ll);
  my $orf_id = $lls[0];
  next unless ($lls[10] <= $cutoff_e);
  next if ($orf_id eq $last_ORF); #### only top hit

  my $rid = $lls[1];
  $ref_ids{ $rid } = 1;
  my $iden = $lls[2];
  my $alnln= $lls[3];  
  my $ref_len = $alnln;
  if ($rid =~ /len\|(\d+)/) {
    $ref_len = $1;
  }
  my $frac = int($alnln / $ref_len * 10000) / 10000;
  $orf_2_hit{$orf_id} = [$rid, $iden, $alnln, $frac];

  $last_ORF = $orf_id;
}
close(TMP);

my %ref_2_taxids = ();
my %ref_2_KO = ();
my %ref_2_ann = ();
open(TMP, $clstr_ref) || die "can not open $clstr_ref";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    chop($ll);
    my ($rid, $no1, $no_taxid, $KO, $des) = split(/\t/, substr($ll, 1));
    next unless ($ref_ids{$rid});

    $ref_2_ann{$rid} = $des;
    $ref_2_taxids{$rid} = [];
    if ($KO =~ /KO\|(\w+)/) {
      $ref_2_KO{$rid} = $1;
    }

    for ($i=0; $i<$no_taxid; $i++) {
      $ll=<TMP>; chop($ll);
      my @lls = split(/\t/, $ll);
      push(@{ $ref_2_taxids{$rid} }, $lls[1]);
    }
  }
}
close(TMP);


my %scaffold_2_taxid = ();
my %taxid_member_scaffolds = ();
my %taxid_orf_count = ();

open(LOG, "> $output_log") || die "can not write to $output_log";
#### assign taxid to scaffold
foreach $round (qw/1 2/) {
  #### first round assign scaffolds that can be uniquely assigned to taxid
  #### 2nd round if multiple taxids match a scaffold with same score, the taxid
  #### got more orfs get priority
  print LOG "####\tround $round to assign scaffolds\n";
  foreach $sid (@all_sids) {
    next if (defined $scaffold_2_taxid{$sid}); #### if sid got assigned in the first round

    my @orf_ids = @{$scaffold_member_orfs{$sid}};
    my $num_orfs = $#orf_ids+1;
    print LOG ">scaffold|$sid\tORFs:$num_orfs\n";
  
    my %taxid_score = ();
    foreach $orf_id (@orf_ids) {
      next unless defined( $orf_2_hit{$orf_id} );
  
      my ($rid, $iden, $alnln, $frac) = @{ $orf_2_hit{$orf_id} };
      my $score = $iden * $alnln;
      print LOG "\tORF:$orf_id\t$rid\t$iden%\t$alnln\n";
  
      next unless (defined($ref_2_taxids{$rid}));
      my @t_taxids = @{ $ref_2_taxids{$rid} }; my $no1 = $#t_taxids+1;
      if ($no1 > 0) { foreach $i (@t_taxids) { $taxid_score{$i} += $score / $no1; } }
    }
  
    next unless (keys %taxid_score);
    my @taxid_score = keys %taxid_score; @taxid_score = sort {$taxid_score{$b} <=> $taxid_score{$a}} @taxid_score;
    
    foreach $i (@taxid_score) { print LOG "\ttaxid\t$i\t$taxid_score{$i}\n"}

    my $tid;
    if ($#taxid_score == 0) { #### only one taxid
      $tid = $taxid_score[0];
    }
    elsif ( $taxid_score{ $taxid_score[0]} > $taxid_score{ $taxid_score[1]} ) { #### first taxid score > 2nd
      $tid = $taxid_score[0];
    }
    elsif ($round == 2) {
      foreach $ii (@taxid_score) { $taxid_orf_count{$ii} = 0 unless (defined($taxid_orf_count{$ii})); }
      @taxid_score = sort { (    $taxid_score{$b} <=>     $taxid_score{$a}) or 
                            ($taxid_orf_count{$b} <=> $taxid_orf_count{$a}) } @taxid_score;
      $tid = $taxid_score[0];
    }
    else { #### next if this is first round and no uniq top taxid
      next;
    }

    $scaffold_2_taxid{$sid} = $tid;
    if (not defined( $taxid_member_scaffolds{$tid} )) {
      $taxid_member_scaffolds{$tid} = [];
      $taxid_orf_count{$tid} = 0;
    }
    push(@{ $taxid_member_scaffolds{$tid} }, $sid);
    $taxid_orf_count{$tid} += $scaffold_orf_count{$sid};
  }
}

close(LOG);



######################### output annotation table
open(OUT, "> $output_ann") || die "can not write to $output_ann";
open(TAX, "> $output_tax") || die "can not write to $output_tax";
print TAX "#Species_taxid\tSpecies\tGenome_taxid\tGenome\tNumber_scaffolds\tNumber_ORFs\n";
print OUT "#Species_taxid\tSpecies\tGenome_taxid\tGenome\tScaffold\tORF\tStart\tEnd\tFrame\tIden%\tFrac_alignment\tFamily\tDescription\n";
#### output annotation with taxid
my @all_tids = keys %taxid_member_scaffolds;
   @all_tids = sort { $taxid_orf_count{$b} <=> $taxid_orf_count{$a} } @all_tids;

my %sptid_member_tids = ();
my %sptid_orf_count = ();
foreach $tid (@all_tids) {
  my @tid_info = @{$taxon_info{$tid}};
  my $sptid = $tid_info[13];
  $sptid = "None" unless ($sptid);
  if (not defined($sptid_member_tids{$sptid})) {
    $sptid_member_tids{$sptid} = [];
    $sptid_orf_count{$sptid} = 0;
  }
  push(@{$sptid_member_tids{$sptid}}, $tid);
  $sptid_orf_count{$sptid} += $taxid_orf_count{$tid};
}
my @all_sptids = keys %sptid_member_tids;
   @all_sptids = sort { $sptid_orf_count{$b} <=> $sptid_orf_count{$a} } @all_sptids;

foreach $sptid (@all_sptids) {
  foreach $tid (@{ $sptid_member_tids{$sptid} }) {
    my @sids = @{ $taxid_member_scaffolds{$tid} };
       @sids = sort { $scaffold_orf_count{$b} <=> $scaffold_orf_count{$a} } @sids;
    my @tid_info = @{$taxon_info{$tid}};
    print TAX "$tid_info[14]\t$tid_info[13]\t$tid\t$tid_info[0]\t", $#sids+1, "\t$taxid_orf_count{$tid}\n";
    foreach $sid (@sids) {
      my @orf_ids = @{$scaffold_member_orfs{$sid}};
      foreach $orf_id (@orf_ids) {
        my $ann = "hypothetical protein";
        my $KO  = "";
        my $iden1 = "-";
        my $frac1 = 1;
        if ( defined( $orf_2_hit{$orf_id} ) ) {
          my ($rid, $iden, $alnln, $frac) = @{ $orf_2_hit{$orf_id} };
          $ann = $ref_2_ann{$rid}; 
          $iden1 = "$iden%";
          $frac1 = $frac;
          $KO = $ref_2_KO{$rid} if (defined($ref_2_KO{$rid}));
        }
        print OUT "$tid_info[14]\t$tid_info[13]\t$tid\t$tid_info[0]\t$sid\t$orf_id\t$orf_info{$orf_id}\t$iden1\t$frac1\t$KO\t$ann\n";
      }
    }
  }
}

#### output scaffolds without taxid
my $no_unknown_sid = 0;
my $no_unknown_orf = 0;
foreach $sid (@all_sids) {
  next if defined($scaffold_2_taxid{$sid});
  my $tid = "Unknown";
  my @orf_ids = @{$scaffold_member_orfs{$sid}};
  foreach $orf_id (@orf_ids) {
    my $ann = "hypothetical protein";
    my $KO  = "";
    my $iden1 = "-";
    my $frac1 = "-";
    if ( defined( $orf_2_hit{$orf_id} ) ) {
      my ($rid, $iden, $alnln, $frac) = @{ $orf_2_hit{$orf_id} };
      $ann = $ref_2_ann{$rid}; 
      $iden1 = "$iden%";
      $frac1 = $frac;
      $KO = $ref_2_KO{$rid} if (defined($ref_2_KO{$rid}));
    }
    print OUT "$tid\tUnknown\tUnknown\tUnknown\t$sid\t$orf_id\t$orf_info{$orf_id}\t$iden1\t$frac1\t$KO\t$ann\n";
  }
  $no_unknown_sid++;
  $no_unknown_orf += $scaffold_orf_count{$sid};
}
print TAX "unknown\t$no_unknown_sid\t$no_unknown_orf\n";
close(OUT);
close(TAX);


sub ORF_info {
  my $ll = shift;
  my ($start, $end, $frame) = ("-","-","-");
#metagene format
#>MSA-1000-10-SE|scaffold|248.1 /source=MSA-1000-10-SE|scaffold|248 /start=2 /end=667 /frame=2 /length=221
#>MSA-1000-10-SE|scaffold|248.8 /source=MSA-1000-10-SE|scaffold|248 /start=7784 /end=8401 /frame=-1 /length=205

  if ( $ll =~ /\/start=(\d+)\s+\/end=(\d+)\s+\/frame=(\S+)/ ) {
    $start = $1; $end = $2; $frame = $3;
  }
  elsif ( $ll =~ /\s+# (\d+) # (\d+) # (\S+) # ID/) {
    $start = $1; $end = $2; $frame = $3;
  }
#prodigal format
#>MSA-1000-10-SE|scaffold|248_1 # 2 # 667 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.677
#>MSA-1000-10-SE|scaffold|248_6 # 6308 # 6658 # -1 # ID=1_6;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.687

  return ($start, $end, $frame); 
}

sub usage {
<<EOD;
$script_name -i blast_alignment_file -r cluster_info -a input ORF -o output ORF file -t tax_file

  options:
    -i blast alignment file in tab format
       can also be a name of a directory, which has multiple blast alignment files 
    -r cluster information file
    -a input ORF fasta file, used in blast search
    -t taxon info file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl based on blast ref db
    -o output prefix , the script creates following files
       output-ann.txt
       output-tax.txt
       output-tax.txt
    -e expect_cutoff, default 1e-6
EOD
}
