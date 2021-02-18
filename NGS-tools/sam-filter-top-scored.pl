#!/usr/bin/perl
## =========================== NGS tools ==========================================
## NGS tools for metagenomic sequence analysis
## May also be used for other type NGS data analysis
##
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
## http://weizhongli-lab.org/
## ================================================================================

#### copy of sam-filter-top-pair.pl, but allow single mapped reads,
#### still don't allow reads that mapped to different ref, i.e. 6th field is not "="

use Getopt::Std;
getopts("T:o:F:O:",\%opts);
die usage() unless ($opts{T});

my $map_score_cutoff = $opts{T};
my $output           = $opts{o};
my $output_ID        = $opts{O};

my ($i, $j, $k, $ll, $cmd);

my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $output";
  $fh = "OUT";
}

if (defined($output_ID)) {
  open(OUTID, "> $output_ID") || die "Can not write to $output_ID";
}


my $last_id = "";
my @R1_alns = ();
my @R2_alns = ();
my $read_pair = "";

while($ll=<>){
  if ($ll =~ /^\@/) { #### headers
    print $fh $ll;
    next;
  }
  else { #### alignment
    chop($ll);
    my @lls = split(/\t/,$ll);
    my $id  = $lls[0];
    my $rid = $lls[2];    if ($rid eq "*") {  next; }
    my $mate_rid = $lls[6];

    my $score;
    for ($j=11; $j<=$#lls; $j++) { if ($lls[$j] =~ /AS:i:(.+)/) { $score = $1; last; } }
    next unless (defined($score));
    next unless ($score >= $map_score_cutoff);

    my $flag_info = <<EOD;
        Chr     Flag    Description
        p       0x0001  the read is paired in sequencing
        P       0x0002  the read is mapped in a proper pair
        u       0x0004  the query sequence itself is unmapped
        U       0x0008  the mate is unmapped
        r       0x0010  strand of the query (1 for reverse)
        R       0x0020  strand of the mate
        1       0x0040  the read is the first read in a pair
        2       0x0080  the read is the second read in a pair
        s       0x0100  the alignment is not primary
        f       0x0200  QC failure
        d       0x0400  optical or PCR duplicate
EOD

    my $flag      = $lls[1];
    my $flag_pe   = $flag & 0x0001;
    my $flag_prop = $flag & 0x0002;
    my $flag_R1   = $flag & 0x0040;
    my $flag_R2   = $flag & 0x0080;

    if    ($id    ne $last_id) {
      if ($read_pair eq "SE") {
        @R1_alns = sort { $b->[2] <=> $a->[2] } @R1_alns;
        my $top_score = $R1_alns[0]->[2];
        foreach $i (@R1_alns) {
          last if ($i->[2] < $top_score);
          print $fh "$i->[0]\n";
          print OUTID "$last_id\n" if ($output_ID);
        }
      }
      elsif ($read_pair eq "PE") {
        my %rid_score_R1 = ();
        my %rid_score_R2 = ();
        my %rid_proper_pair = ();
        my %rid_score = ();
        my %rid_aln_R1 = ();
        my %rid_aln_R2 = ();
        my %rid_score_source = ();

        foreach $i (@R1_alns) {
          @f = @{$i};
          $j = $f[1];
          if (( not defined($rid_score_R1{$j})) or ($f[2] > $rid_score_R1{$j})) { 
             $rid_score_R1{$j} = $f[2];
             $rid_proper_pair{$j} = $f[5];
             $rid_aln_R1{$j} = $f[0];
          }
          $rid_score{$j} = 0;
        }
        foreach $i (@R2_alns) {
          @f = @{$i};
          $j = $f[1];
          if (( not defined($rid_score_R2{$j})) or ($f[2] > $rid_score_R2{$j})) {      
             $rid_score_R2{$j} = $f[2];
             $rid_proper_pair{$j} = $f[5];
             $rid_aln_R2{$j} = $f[0];
          }
          $rid_score{$j} = 0;
        }

        foreach $i (keys %rid_score) {
          $rid_score_R1{$i} = 0 unless (defined($rid_score_R1{$i}));
          $rid_score_R2{$i} = 0 unless (defined($rid_score_R2{$i}));
          if ($rid_proper_pair{$i} eq "=") {
            $rid_score{$i} = $rid_score_R1{$i} + $rid_score_R2{$i};
            $rid_score_source{$i} = "R12";
          }
          else {
            $rid_score{$i}        = $rid_score_R1{$i} > $rid_score_R2{$i} ? $rid_score_R1{$i} : $rid_score_R2{$i};
            $rid_score_source{$i} = $rid_score_R1{$i} > $rid_score_R2{$i} ? "R1" : "R2";
            
          }
        }

        my @scores = values %rid_score;
           @scores = sort{$b<=>$a} @scores;
        my $top_score = $scores[0];
        my @out_rids = ();

        foreach $i (keys %rid_score) {
          push(@out_rids, $i) if ($rid_score{$i} >= $top_score);
        }

        #### print R1 aln
        foreach $i (@out_rids) {
          if (($rid_score_source{$i} eq "R12") or ($rid_score_source{$i} eq "R1")) {
            print $fh "$rid_aln_R1{$i}\n";
          }
        }
        #### print R2 aln
        foreach $i (@out_rids) {
          if (($rid_score_source{$i} eq "R12") or ($rid_score_source{$i} eq "R2")) {
            print $fh "$rid_aln_R2{$i}\n";
          }
        }
      }
      else {
        print STDERR "no pair info for $last_id\n";
      }

      @R1_alns = ();
      @R2_alns = ();
      $read_pair = "";
    }

    $read_pair =  $flag_pe ? "PE" : "SE";
    if ($read_pair eq "SE") { push(@R1_alns, [$ll, $rid, $score]); }
    elsif ($flag_R2)        { push(@R2_alns, [$ll, $rid, $score, $flag_pe, $flag_prop, $mate_rid]); }
    else                    { push(@R1_alns, [$ll, $rid, $score, $flag_pe, $flag_prop, $mate_rid]); }
    $last_id = $id;

  } #### alignment section
}

   #### last reads
    if    ($last_id) {
      if ($read_pair eq "SE") {
        @R1_alns = sort { $b->[2] <=> $a->[2] } @R1_alns;
        my $top_score = $R1_alns[0]->[2];
        foreach $i (@R1_alns) {
          last if ($i->[2] < $top_score);
          print $fh "$i->[0]\n";
          print OUTID "$last_id\n" if ($output_ID);
        }
      }
      elsif ($read_pair eq "PE") {
        my %rid_score_R1 = ();
        my %rid_score_R2 = ();
        my %rid_proper_pair = ();
        my %rid_score = ();
        my %rid_aln_R1 = ();
        my %rid_aln_R2 = ();
        my %rid_score_source = ();

        foreach $i (@R1_alns) {
          @f = @{$i};
          $j = $f[1];
          if (( not defined($rid_score_R1{$j})) or ($f[2] > $rid_score_R1{$j})) { 
             $rid_score_R1{$j} = $f[2];
             $rid_proper_pair{$j} = $f[5];
             $rid_aln_R1{$j} = $f[0];
          }
          $rid_score{$j} = 0;
        }
        foreach $i (@R2_alns) {
          @f = @{$i};
          $j = $f[1];
          if (( not defined($rid_score_R2{$j})) or ($f[2] > $rid_score_R2{$j})) {      
             $rid_score_R2{$j} = $f[2];
             $rid_proper_pair{$j} = $f[5];
             $rid_aln_R2{$j} = $f[0];
          }
          $rid_score{$j} = 0;
        }

        foreach $i (keys %rid_score) {
          $rid_score_R1{$i} = 0 unless (defined($rid_score_R1{$i}));
          $rid_score_R2{$i} = 0 unless (defined($rid_score_R2{$i}));
          if ($rid_proper_pair{$i} eq "=") {
            $rid_score{$i} = $rid_score_R1{$i} + $rid_score_R2{$i};
            $rid_score_source{$i} = "R12";
          }
          else {
            $rid_score{$i}        = $rid_score_R1{$i} > $rid_score_R2{$i} ? $rid_score_R1{$i} : $rid_score_R2{$i};
            $rid_score_source{$i} = $rid_score_R1{$i} > $rid_score_R2{$i} ? "R1" : "R2";
            
          }
        }

        my @scores = values %rid_score;
           @scores = sort{$b<=>$a} @scores;
        my $top_score = $scores[0];
        my @out_rids = ();

        foreach $i (keys %rid_score) {
          push(@out_rids, $i) if ($rid_score{$i} >= $top_score);
        }

        #### print R1 aln
        foreach $i (@out_rids) {
          if (($rid_score_source{$i} eq "R12") or ($rid_score_source{$i} eq "R1")) {
            print $fh "$rid_aln_R1{$i}\n";
          }
        }
        #### print R2 aln
        foreach $i (@out_rids) {
          if (($rid_score_source{$i} eq "R12") or ($rid_score_source{$i} eq "R2")) {
            print $fh "$rid_aln_R2{$i}\n";
          }
        }
      }
      else {
        print STDERR "no pair info for $last_id\n";
      }

      @R1_alns = ();
      @R2_alns = ();
      $read_pair = "";
    }


close(OUT) if ($output);
close(OUTID) if ($output_ID);

sub sam_flag {
  my $flag = shift;
  my $flag_info = <<EOD;
        Chr     Flag    Description
        p       0x0001  the read is paired in sequencing
        P       0x0002  the read is mapped in a proper pair
        u       0x0004  the query sequence itself is unmapped
        U       0x0008  the mate is unmapped
        r       0x0010  strand of the query (1 for reverse)
        R       0x0020  strand of the mate
        1       0x0040  the read is the first read in a pair
        2       0x0080  the read is the second read in a pair
        s       0x0100  the alignment is not primary
        f       0x0200  QC failure
        d       0x0400  optical or PCR duplicate
EOD

  my $oflag = "";
  $oflag .= ($flag & 0x0001 ) ? "_P" : "_p";
  $oflag .= ($flag & 0x0002 ) ? "_M" : "_m";
  $oflag .= ($flag & 0x0004 ) ? "_u" : "_.";
  $oflag .= ($flag & 0x0008 ) ? "_U" : "_.";
  $oflag .= ($flag & 0x0010 ) ? "_r" : "_.";
  $oflag .= ($flag & 0x0020 ) ? "_R" : "_.";
  $oflag .= ($flag & 0x0040 ) ? "_1" : "_.";
  $oflag .= ($flag & 0x0080 ) ? "_2" : "_.";
  $oflag .= ($flag & 0x0100 ) ? "_s" : "_S";
  $oflag .= ($flag & 0x0200 ) ? "_f" : "_F";
  $oflag .= ($flag & 0x0400 ) ? "_d" : "_D";

  return $oflag;
}

sub usage {
<<EOD
Given a sam file generated by 
  bwa mem reference R1.fa_or_fq R2.fa_or_fq 
  bwa mem reference R1.fa_or_fq


This script takes the sam from STDIN, first filter out alignments below the alignment_cutoff_score,
then prints out top scored alignments
  For PE reads, the score is the sum of R1 and R2 for properly mapped pair, or the max of R1 and R2 otherwise
     it prints out both R1 and R2 for properly mapped pair, or max of R1/R2 otherwise
     it prints out multiple alignments (with different reference) of the same top score

usage:
  $script_name  -T alignment_cutoff_score

  options
    -o output file, default STDOUT
    -T alignment cutoff score, defined in sam as AS:i:score
EOD
}
########## END usage

