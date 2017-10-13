#!/usr/bin/perl
## ==============================================================================
## Automated annotation tools
##
## program written by
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
##                                      http://weizhong-lab.ucsd.edu
## ==============================================================================

sub LL_random_sleep {
  my $a = shift;
     $a = 10 unless (defined ($a));
     $a = rand() * $a;
     select(undef,undef,undef,$a);
}
# END LL_random_sleep


sub LL_get_active_ids {
  my $seq_dir = shift;
  opendir(DIR, $seq_dir) || die "Can not open seq dir";
  my @ids = grep {/^\w/} readdir(DIR);
  closedir(DIR);
  return @ids;
}
# END LL_get_active_ids


sub LL_shuffle_array{
  my $ref = shift;
  my $no = $#{$ref};
  my ($i, $j, $k);
                                                                                
  for ($i=0; $i<$no/2; $i++) {
    $j  = int($no * rand());
    $k  = int($no * rand());
                                                                                
    my    $tmp = $ref->[$j];
    $ref->[$j] = $ref->[$k];
    $ref->[$k] = $tmp;
  }
}
# END LL_shuffle_array


sub LL_make_local_blast_db {
  my ($local_tmpdb, @dbs) = @_;
  my ($i, $j, $k, $cmd);

  if (not -e $local_tmpdb) { $cmd = `mkdir -p $local_tmpdb`; }
  if (not -e $local_tmpdb) { die "can not mkdir $local_tmpdb"; }

  foreach $i (@dbs) {
    $cmd = `rsync -av $ENV{BLASTDB}/$i.* $local_tmpdb`;
  }
}
# END LL_make_local_blast_db

sub LL_cpu_time {
    my ($tu, $ts, $cu, $cs) = times();
    my $tt = $tu + $ts + $cu + $cs;
    return $tt;
}
# END LL_cpu_time

sub blastdb_after_file2 {
  my ($db, $file2) = @_;
  my $db_file = "$db.pin"; ###### need to be updated
  return file1_after_file2($db_file, $file2);
}

# file1 is newer then file2
sub file1_after_file2 {
  my ($file1, $file2) = @_;

  # if not exist file1, assume it is in future, so it is newer
  if (not -e ($file1)) {return 0;}
  if (not -e ($file2)) {return 0;}

  my $mtime1 = (stat($file1))[9];
  my $mtime2 = (stat($file2))[9];

  return ( ($mtime1 > $mtime2) ? 1 : 0);
}

  # calling get_mode_by_time_stamp($blast_dir, $last_sql_txt)
  # $last_sql_txt is last xxx.txt file mysqlimported
  # $blast_dir    is output dir of results
  # return mode what is mode
  # mode can be reload or append or none
  # reload means that all output in $blast_dir are newer than in mysql table
  # append means that some       in $blast_dir are newer than in mysql table
  # none   menas that all       in $blast_dir are older than in mysql table
sub get_mode_by_time_stamp {
  my ($dir, $file) = @_;
  my ($i, $j, $k);
  my $mode = "";
  opendir(BLDIR, $dir) || return $mode;
  my @fs = grep {/^\w/} readdir(BLDIR);
  close(BLDIR);

  (-e $file) || return "reload";
  my $mfile = (stat($file))[9];
  my $min = 1093562513 * 1093562513 ;
  my $max = -1;
  foreach $i (@fs) {
    $j = "$dir/$i";
    $k = (stat($j))[9];
    $min = $k if ($k<$min);
    $max = $k if ($k>$max);
  }

  if    ($max < $mfile) {return "none";}
  elsif ($min > $mfile) {return "reload";}
  else                  {return "append";}
}
#### END 

sub LL_get_run_id {
  my $app = shift;
  my $id_file = "mysql/app_run.txt.loaded";

  my $run_id = "$app.001";

  open(TID, $id_file) || return $run_id;
  my $ll;
  # this way return the last run_id if multiple run exists
  while($ll = <TID>) {
    chop($ll);
    my ($tid, $tapp) = split(/\t/,$ll);
    $run_id = $tid if ($tapp eq $app);
  }
  close(TID);
  return $run_id;
}


sub get_sql_date {
  my $file = shift;
  if ($file) {
    my $cmd = `stat $file | grep Modify`;
    my $mtime = substr($cmd, 8, 10);
    return $mtime;
  }
  else {
    my $date = `date +%Y_%m_%d`;
       $date =~ s/\s//g;
    return $date;
  }
}

sub get_sql_current_date {
  my $date = `date +%Y_%m_%d`;
     $date =~ s/\s//g;
  return $date;
}


1;
