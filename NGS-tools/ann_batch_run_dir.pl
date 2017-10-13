#!/usr/bin/perl
## ==============================================================================
## Automated annotation tools
##
## program written by
##                                      Weizhong Li, UCSD
##                                      liwz@sdsc.edu
##                                      http://weizhong-lab.ucsd.edu
## ==============================================================================

my $readme = <<EOD;
given a dir of files, this script run each of the file
call this script as
--INDIR1=/home/my/dir1 --INDIR2=/home/my/dir2 --OUTDIR1=/home/my/outputdir1 --OUTDIR2=/home/my/outputdir2 some_app -d nr -i {INDIR1} -j {INDIR2} -o {OUTDIR1} -z {OUTDIR2} 

--INDIR1 is the master input dir
if /home/my/dir1 has 1024 files, such as 0, 1, 2, ... 1023
this script will run 1024 commands
some_app -d nr -i /home/my/dir1/0 -j /home/my/dir2/0 -o /home/my/outputdir1 -z /home/my/outputdir2
some_app -d nr -i /home/my/dir1/1 -j /home/my/dir2/1 -o /home/my/outputdir1 -z /home/my/outputdir2
some_app -d nr -i /home/my/dir1/2 -j /home/my/dir2/2 -o /home/my/outputdir1 -z /home/my/outputdir2
...

some_app -d nr -i /home/my/dir1/1023 -j /home/my/dir2/1023 -o /home/my/outputdir1 -z /home/my/outputdir2

if your command has ">" "<" or "|", please use "\>" "\<" "\|"

if there are multiple input dirs such as --INDIR2, this dir should contains corrsponding 1024 files with same name as the master input dir, INDIR1
otherwise, command will skip
EOD

my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   $script_dir = "./" unless ($script_dir);
require "$script_dir/ann_local.pl";


my $cpu_file = "cpu.log";
my $cmd_line = "";
my $merge_output = "";
my %dir = ();

my ($i, $j, $k, $ll, $cmd, $arg, $f1, $f2);

while(1) {
  $arg = shift;
  last unless (defined($arg));

  if (($arg =~ /^--INDIR\d+=/) or ($arg =~ /^--OUTDIR\d+=/)) {
    $i = substr($arg, 2);
    ($j, $k) = split(/=/, $i);
    $dir{$j} = $k;
  }
  elsif ($arg =~ /^--CPU=/) {
    $cpu_file = substr($arg,6); 
  }
  elsif ($arg =~ /^--MERGE/) {
    $merge_output = 1;
  }
  else {
    $cmd_line .= $arg . " ";
  }
}
$cmd_line =~ s/\s+$//;
print STDERR "command template is $cmd_line\n";


my $master_dir = $dir{"INDIR1"};
die "No master input dir" unless ($master_dir);
die "No output dir" unless ($dir{"OUTDIR1"});

LL_random_sleep(10);
my @seqs = LL_get_active_ids($master_dir); LL_shuffle_array(\@seqs);

foreach $t_dir (keys %dir) {
  next unless ($t_dir =~ /^OUTDIR\d+/);
  $cmd = `mkdir -p $dir{$t_dir}` unless (-e $dir{$t_dir});
}
LL_random_sleep(10);

foreach $i (@seqs) {
  LL_random_sleep(1);
  my $cmd1 = $cmd_line;

  my $file_ready_flag = 1;
  my $file_not_ready = "";
  foreach $t_dir (keys %dir) {
    $f1 = "$dir{$t_dir}/$i";
    $cmd1 =~ s/\{$t_dir\}/$f1/g;
    if ($t_dir =~ /^INDIR\d+/) {
      if (not (-e $f1)) {
        $file_ready_flag = 0;
        $file_not_ready = $f1;
      }
    }
  }

  if ($file_ready_flag == 0) {
    print STDERR "skip \"$cmd1\", file $file_not_ready does not exist\n";
    next;
  }

  my $pri_output = "$dir{'OUTDIR1'}/$i";

  next if (-s $pri_output);        # output file exist
  next if (-s "$pri_output.gz");
  next if (-e "$pri_output.lock"); # being calculated by another parallelly
  print STDERR "$cmd1\n\n";
  $cmd = `date > $pri_output.lock`;
  $cmd = `$cmd1`;
  $cmd = `rm -f $pri_output.lock`;

}

LL_random_sleep(10);

if ($merge_output) {
  my $lockf = "$master_dir.lock";
  if (not (-e $lockf)) {
    $cmd = `date > $lockf`;
    merge_files();
    $cmd = `rm -f $lockf`;
  }
}

my ($tu,$ts,$cu,$cs)=times(); my $tt=$tu+$ts+$cu+$cs;
$cmd = `echo $dir{'OUTDIR1'} $tt >> $cpu_file`;

sub merge_files {
  my ($i, $j, $k, $t_dir);

  foreach $t_dir (keys %dir) {
    next unless ($t_dir =~ /^OUTDIR\d+/);
    my $dir = $dir{$t_dir};
    my $merge_out = "$dir.merge";
    opendir(DIR, $dir) || die "can not open dir $dir";
    my @files = sort grep {/\w/} readdir(DIR);
    closedir(DIR);

    foreach $i (@files) {
      my $cmd = `cat $dir/$i >> $merge_out`;
    }    
  }  
}




















