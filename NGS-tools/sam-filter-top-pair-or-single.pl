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
getopts("T:o:F:",\%opts);
die usage() unless ($opts{T});

my $map_score_cutoff = $opts{T};
my $output           = $opts{o};
my $proper_flag      = $opts{F};
   $proper_flag      =1 unless (defined($proper_flag));

my ($i, $j, $k, $ll, $cmd);

my $fh = "STDOUT";
if ($output) {
  open(OUT, ">$output") || die "Can not write to $OUT";
  $fh = "OUT";
}

my $last_id = "";
my $R1_pri_aln = "";
my $R2_pri_aln = "";

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

    my $FLAG = $lls[1];
    if ($FLAG & 0x0001 ) { #### PE
      if ($proper_flag) {
        next unless ($FLAG & 0x0002 );
      }
    }

    if    ($id    ne $last_id) {
      my $print_flag = 1;
      if ($R1_pri_aln and $R2_pri_aln) {
        my @t_R1 = split(/\t/, $R1_pri_aln);
        my @t_R2 = split(/\t/, $R2_pri_aln);
        if ($proper_flag) {
          $print_flag = 0 unless ($t_R1[2] eq $t_R2[2]);
          $print_flag = 0 unless ($t_R1[6] eq "=");
          $print_flag = 0 unless ($t_R2[6] eq "=");
        }
        else { #### keep R1 only if R1 and R2 doesn't map to same ref
          $R2_pri_aln = "" unless ($t_R1[2] eq $t_R2[2]);
          $R2_pri_aln = "" unless ($t_R1[6] eq "=");
          $R2_pri_aln = "" unless ($t_R2[6] eq "=");
        }
      }
      elsif ($R1_pri_aln or $R2_pri_aln) {
        $print_flag = 1; #### only one read mapped
      }
      else {
        $print_flag = 0;
      }

      if ($print_flag) {
        print $fh "$R1_pri_aln\n" if ($R1_pri_aln);
        print $fh "$R2_pri_aln\n" if ($R2_pri_aln);
      }
      $R1_pri_aln = "";
      $R2_pri_aln = "";
    }

    #### only 1 alignment per R1 / R2 end
    if    ( $FLAG & 0x0040 ) { $R1_pri_aln = $ll unless ($R1_pri_aln); }
    elsif ( $FLAG & 0x0080 ) { $R2_pri_aln = $ll unless ($R2_pri_aln); }
    $last_id = $id;

  } #### alignment section
}

   #### last reads
    if    ($id    ne $last_id) {
      my $print_flag = 1;
      if ($R1_pri_aln and $R2_pri_aln) {
        my @t_R1 = split(/\t/, $R1_pri_aln);
        my @t_R2 = split(/\t/, $R2_pri_aln);
        if ($proper_flag) {
          $print_flag = 0 unless ($t_R1[2] eq $t_R2[2]);
          $print_flag = 0 unless ($t_R1[6] eq "=");
          $print_flag = 0 unless ($t_R2[6] eq "=");
        }
        else { #### keep R1 only if R1 and R2 doesn't map to same ref
          $R2_pri_aln = "" unless ($t_R1[2] eq $t_R2[2]);
          $R2_pri_aln = "" unless ($t_R1[6] eq "=");
          $R2_pri_aln = "" unless ($t_R2[6] eq "=");
        }
      }
      elsif ($R1_pri_aln or $R2_pri_aln) {
        $print_flag = 1; #### only one read mapped
      }
      else {
        $print_flag = 0;
      }

      if ($print_flag) {
        print $fh "$R1_pri_aln\n" if ($R1_pri_aln);
        print $fh "$R2_pri_aln\n" if ($R2_pri_aln);
      }
      $R1_pri_aln = "";
      $R2_pri_aln = "";
    }


close(OUT) if ($output);

sub sam_flag {
  my $FLAG = shift;
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
  $oflag .= ($FLAG & 0x0001 ) ? "_P" : "_p";
  $oflag .= ($FLAG & 0x0002 ) ? "_M" : "_m";
  $oflag .= ($FLAG & 0x0004 ) ? "_u" : "_.";
  $oflag .= ($FLAG & 0x0008 ) ? "_U" : "_.";
  $oflag .= ($FLAG & 0x0010 ) ? "_r" : "_.";
  $oflag .= ($FLAG & 0x0020 ) ? "_R" : "_.";
  $oflag .= ($FLAG & 0x0040 ) ? "_1" : "_.";
  $oflag .= ($FLAG & 0x0080 ) ? "_2" : "_.";
  $oflag .= ($FLAG & 0x0100 ) ? "_s" : "_S";
  $oflag .= ($FLAG & 0x0200 ) ? "_f" : "_F";
  $oflag .= ($FLAG & 0x0400 ) ? "_d" : "_D";

  return $oflag;
}

sub usage {
<<EOD
Given a sam file generated by 
  bwa mem reference R1.fa_or_fq R2.fa_or_fq 

This script takes the sam from STDIN and prints out top score paired mapped reads  /
  or single mapped reads that are above alignment score cutoff.

usage:
  $script_name  -T alignment_cutoff_score

  options
    -o output file, default STDOUT
    -T alignment cutoff score, defined in sam as AS:i:score
    -F flag to remove inproper pair, default 1
       proper pair has FLAG 0x0002
       if the purpose is to keep good alignments, suggest use default value,
       if the purpose is to keep alignments for filtering, e.g. remove-host, suggest use 0.
EOD
}
########## END usage

