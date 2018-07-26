#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================


## copied from  assembly-binning.pl 
## added taxon parsing function
## parse from taxid and get species taxid from taxon file
use Getopt::Std;
getopts("i:j:c:s:o:d:r:n:a:t:x:",\%opts);
die usage() unless ($opts{i} and $opts{j} and $opts{o} and $opts{s} and $opts{a} and $opts{t});

my $sam_assembly     = $opts{i};
my $sam_reference    = $opts{j};
my $assembly_file    = $opts{s};
my $output           = $opts{o};
my $n_cutoff         = $opts{n}; $n_cutoff = 10  unless defined($n_cutoff);
my $p_cutoff         = $opts{c}; $p_cutoff = 0.5 unless defined($p_cutoff);
my $tax_str          = $opts{a};
my $taxon_file       = $opts{t}; #### taxon file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl
my $contaminant_file = $opts{x}; 

my ($i, $j, $k, $ll, $cmd);

my %contaminant_tids = ();
if ($contaminant_file) {
  open(TMP, $contaminant_file) || die "can not open $contaminant_file";
  while($ll=<TMP>) {
    chop($ll);
    next if ($ll =~ /^#/);
    my @lls = split(/\s+/, $ll);
    $contaminant_tids{$lls[0]} = 1;
  }
  close(TMP);
}

my %tid_2_sptid = ();
my %taxon_info = ();
my %sp_name = ();
my %strain_name = ();
open(TMP, $taxon_file) || die "can not open $taxon_file";
my $taxon_format = <<EOD;
Col.0   #taxid  511145
Col.1   rank    toprank
Col.2   name    Escherichia coli str. K-12 substr. MG1655
Col.3   superkingdom    Bacteria
Col.4   superkingdom_ti 2
Col.5   kingdom \N
Col.6   kingdom_ti      \N
Col.7   phylum  Proteobacteria
Col.8   phylum_ti       1224
Col.9   class   Gammaproteobacteria
Col.10  class_ti        1236
Col.11  order   Enterobacterales
Col.12  order_ti        91347
Col.13  family  Enterobacteriaceae
Col.14  family_ti       543
Col.15  genus   Escherichia
Col.16  genus_ti        561
Col.17  species Escherichia coli
Col.18  species_ti      562
Col.19  toprank Escherichia coli str. K-12 substr. MG1655
Col.20  toprank_ti      511145
EOD

while($ll=<TMP>) {
  chop($ll);
  next if ($ll =~ /^#/);
  my ($tid, $rank, @lls) = split(/\t/, $ll);
  next unless ($rank eq "toprank");
  $tid_2_sptid{$tid} = $lls[16];
  $taxon_info{$tid} = [@lls];
  $sp_name{$lls[16]} = $lls[15];
  $strain_name{$tid} = $lls[0];

}
close(TMP);


open(TMP, $assembly_file) || die "can not open $assembly_file";
my @assembly_ids = ();
my %assembly_2_len = ();
my $seq_id;
while($ll=<TMP>){
  chop($ll);
  if ($ll =~ /^>/) {
    $seq_id = substr($ll,1);
    $seq_id =~ s/\s.+$//;
    $assembly_2_len{$seq_id} = 0;
    push(@assembly_ids, $seq_id);
  }
  else {
    $ll =~ s/\s+$//;
    $assembly_2_len{$seq_id}  += length($ll);
  }
}
close(TMP);

open(TMP, $sam_assembly) || die "can not open $sam_assembly";
my %assembly_reads_mapped = ();
while($ll=<TMP>){
  if ($ll =~ /^\@/) { #### headers
    next;
  }
  else { #### alignment
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $id  = $lls[0];
    my $rid = $lls[2];    if ($rid eq "*") {  next; }
    my $FLAG = $lls[1];
    # next unless ( $FLAG & 0x0040 ); #### count R1 only 

    if (not defined($assembly_reads_mapped{$rid})) {
      $assembly_reads_mapped{$rid} = [];
    }
    push(@{$assembly_reads_mapped{$rid}}, $id);
  } #### alignment section
}
close(TMP);

open(TMP, $sam_reference) || die "can not open $sam_reference";
my %t_taxids = ();
my $last_id = "";
my %read_2_taxids = ();
while($ll=<TMP>){
  if ($ll =~ /^\@/) { #### headers
    next;
  }
  else { #### alignment
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $id  = $lls[0];
    my $rid = $lls[2];    if ($rid eq "*") {  next; }
    my $tid = "";
    if ($rid =~ /$tax_str\|(\d+)/) { $tid = $1; } else { next; }

    my $FLAG = $lls[1];
    # next unless ( $FLAG & 0x0040 ); #### count R1 only 

    if    (($id  ne $last_id) and $last_id ) {
      my @t_taxids = keys %t_taxids;
      if ($#t_taxids+1 <= $n_cutoff) {
        #### if already defined by R1
        $read_2_taxids{$last_id} = [@t_taxids] unless (defined($read_2_taxids{$last_id}));
      }
      %t_taxids = ();
    }

    $t_taxids{$tid} = 1;
    $last_id = $id;
  } #### alignment section
}
close(TMP);
    if    (($id  ne $last_id) and $last_id ) {
      my @t_taxids = keys %t_taxids;
      if ($#t_taxids+1 <= $n_cutoff) {
        #### if already defined by R1
        $read_2_taxids{$last_id} = [@t_taxids] unless (defined($read_2_taxids{$last_id}));
      }
    }


my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $OUT";
  $fh = "OUT";
}
foreach $i (@assembly_ids) {
  my $taxid = "Unknown";
  my $sptid = "Unknown";
  my $sp_name = "Unknown";
  my $strain_name = "Unknown";
  my $p = 1;
  my $n_reads = 0;

  if ($assembly_reads_mapped{$i}) {
    my @reads = @{$assembly_reads_mapped{$i}};
    $n_reads = $#reads+1;
    my %mapped_tids = ();
    my %mapped_sptids = ();

    foreach $j (@reads) {
      next unless ($read_2_taxids{$j});
      my @tids = @{$read_2_taxids{$j}};
      my $n_tids = $#tids + 1;
      foreach $k (@tids) {
        $mapped_tids{$k} += 1.0 / $n_tids / $n_reads;
        $mapped_sptids{ $tid_2_sptid{$k} } += 1.0 / $n_tids / $n_reads;
      }
    }

    my @mapped_sptids = keys %mapped_sptids;
       @mapped_sptids = sort { $mapped_sptids{$b} <=> $mapped_sptids{$a} } @mapped_sptids;

    if ($mapped_sptids{ $mapped_sptids[0] } >= $p_cutoff) {
      $sptid = $mapped_sptids[0];
      $p = $mapped_sptids{ $mapped_sptids[0] };
      $sp_name = $sp_name{$sptid};

      my @mapped_tids = keys %mapped_tids;
         @mapped_tids = sort { $mapped_tids{$b} <=> $mapped_tids{$a} } @mapped_tids;
      foreach $k (@mapped_tids) {
        next unless ( $tid_2_sptid{$k} eq  $sptid);
        $taxid = $k;
        $strain_name = $strain_name{$k};
        last;
      }
    }

  }

  if ($contaminant_tids{$taxid}) {
     $sp_name = "sptid:$sptid";
     $strain_name = "taxid:$taxid";
     $taxid = "contaminant";
     $sptid = "contaminant";
  }
  print OUT "$i\t$sptid\t$sp_name\t$taxid\t$strain_name\t$p\t$n_reads\t$assembly_2_len{$i}\n";

}
close(OUT) if ($output);


sub usage {
<<EOD
Given a sam file of reads to assembles and 
  another sam file of reads to reference genome generated by 
  bwa mem reference R1.fa_or_fq R2.fa_or_fq 
and (better) filtered by sam-filter-top-pair.pl

This script uses the reads to reference and reads to assembly
mapping information to bin the assembly sequences

usage:
  $script_name -i reads-to-assembly-mapping.sam -j reads-to-reference-mapping.sam -o output -s assembly_fasta_file -a tid -t taxon_file

  options
    -i reads-to-assembly-mapping.sam
    -j reads-to-reference-mapping.sam
    -o output
    -s fasta file of assembly
    -n number alignment cutoff, default 10, 
       if one read is mapped to multiple taxids (> the cutoff), this read is skipped.
    -c probability cutoff, default 0.5
       each read (R) has a probality to a genome (G), if a reads mapped to
       N multiple genomes, the probability is 1/N
       for a assembly (S) , the probability for it belong to the same species of 
       a genome (G) is sum of P(R,G) / number of reads mapped to S
    -a taxid string, the name of reference genome can be ...sequence_name|tid|12345,
       here 12345 is the taxon ID. in this case, use -a tid so that the script can 
       parse the taxon ID 
    -t taxon info file, created by ~/git/ngomicswf/NGS-tools/taxon_print_tid_rank_table.pl based on blast ref db
    -x taxids of contaminants, optional, if exist, will label these as Contaminant

EOD
}
########## END usage

