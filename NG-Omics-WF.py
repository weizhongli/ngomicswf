#!/usr/bin/env python
# =============================== NG-Omics-WF ==================================
#  _   _  _____         ____            _              __          ________ 
# | \ | |/ ____|       / __ \          (_)             \ \        / /  ____|
# |  \| | |  __ ______| |  | |_ __ ___  _  ___ ___ _____\ \  /\  / /| |__   
# | . ` | | |_ |______| |  | | '_ ` _ \| |/ __/ __|______\ \/  \/ / |  __|  
# | |\  | |__| |      | |__| | | | | | | | (__\__ \       \  /\  /  | |     
# |_| \_|\_____|       \____/|_| |_| |_|_|\___|___/        \/  \/   |_|     
#                                                                           
# =========================== Next Generation Omics data workflow tools ========
#
# Workflow tools for next generation genomics, metagenomics, RNA-seq 
# and other type of omics data analyiss, 
#
# Software originally developed since 2010 by Weizhong Li at UCSD
#                                               currently at JCVI
#
# https://github.com/weizhongli/ngomicswf       liwz@sdsc.edu
# ==============================================================================

import os
import sys
import re
import argparse 
from argparse import RawTextHelpFormatter
import math
import subprocess
import logging
import textwrap
import imp
import collections

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

__author__ = 'Weizhong Li'

############## Global variables
NGS_config = None
NGS_samples = []
NGS_sample_data = {}
NGS_opts = {}
pwd = os.path.abspath('.')
subset_flag = False
subset_jobs = []
qstat_xml_data = {}
job_list = collections.defaultdict(dict)  # as job_list[$t_job_id][$t_sample_id] = {}
execution_submitted = {}                  # number of submitted jobs (qsub) or threads (local sh)
############## END Global variables

def read_parameters(args):
  """
  read option parameters from file or command line
  """
  if args.parameter_file:
    try:
      ##format example
      ##JobID_A opt0 opt1 opt2
      ##JobID_B opt0 opt1
      ##JobID_C opt0 opt1 opt2 opt3
      f = open(args.parameter_file, 'r')
      for line in f:
        if line[0] == '#':
          continue
        if not re.match('^w', line):
          continue
        ll = re.split('\s+', line.rstrip());
        NGS_opts[ ll[0] ] = ll[1:]
      f.close()
    except IOError:
      print 'cannot open', args.parameter_file
      exit(1)

  elif args.parameter_name:
    for line in re.split(',', args.parameter_name):
      ll = re.split(':', line);
      NGS_opts[ ll[0] ] = ll[1:]

def read_samples(args):
  """
  read sample and sample data from file or command line
  """
  if args.sample_file:
    try:
      f = open(args.sample_file, 'r')
      for line in f:
        if line[0] == '#':
          continue
        if not re.match('^\w', line):
          continue
        ll = re.split('\s+', line.rstrip());
        NGS_samples.append(ll[0]);
        NGS_sample_data[ ll[0] ] = ll[1:]
      f.close()
    except IOError:
      print 'cannot open', args.sample_file
      exit(1) 

  elif args.sample_name:
    for line in re.split(',', args.sample_name):
      ll = re.split(':', line);
      NGS_samples.append(ll[0]);
      NGS_sample_data[ ll[0] ] = ll[1:]
  else:
    exit(1)

  for sample in NGS_samples:
    if os.path.exists(sample):
      if os.path.isdir(sample):
        pass
      else:
        print 'file exist:', sample
        exit(1)
    else:
      if os.system("mkdir " + sample):
        print 'can not mkdir:', sample
        exit(1)
       
def make_job_list(NGS_config):
  '''
  make sh script for each job / sample
  '''

  verify_flag = False
  for t_job_id in NGS_config.NGS_batch_jobs:
    if subset_flag and not (t_job_id in subset_jobs):
      continue

    print t_job_id
    t_job = NGS_config.NGS_batch_jobs[ t_job_id ]
    t_execution = NGS_config.NGS_executions[ t_job["execution"] ]

    print t_job
    print t_execution
    
    pe_parameter = ''
    if t_execution[ 'type' ] == 'qsub-pe':
      t_cores_per_cmd  = t_job[ 'cores_per_cmd' ]
      pe_parameter = "#$ -pe orte " + str(t_cores_per_cmd)

    if t_job[ 'cores_per_cmd' ] > t_execution[ 'cores_per_node' ]:
      print 'not enough cores'
      print t_job
      exit(1)
      ## -- write_log("$t_job_id needs $t_job->{\"cores_per_cmd\"} cores, but $t_job->{\"execution\"} only has $t_execution->{\"cores_per_node\"} cores");

    t_job[ 'cmds_per_node' ] = t_execution[ 'cores_per_node' ] / t_job[ 'cores_per_cmd' ]
    t_job[ 'nodes_total' ] = math.ceil( t_job[ 'no_parallel' ] / float(t_job[ 'cmds_per_node' ]))
 
    if t_job[ 'nodes_total' ] > t_execution[ 'number_nodes' ]:
      print 'not enough nodes'
      print t_job
      exit(1)
      ## -- write_log("$t_job_id needs $t_job->{\"nodes_total\"} nodes, but $t_job->{\"execution\"} only has $t_execution->{\"number_nodes\"} nodes");

    CMD_opts = []
    if 'CMD_opts' in t_job.keys():  
      CMD_opts = t_job[ 'CMD_opts' ]
    if t_job_id in NGS_opts.keys():
      CMD_opts = NGS_opts[ t_job_id ]

    for t_sample_id in NGS_samples:
      t_command = t_job[ 'command' ]
      t_command = re.sub('\\\\SAMPLE', t_sample_id, t_command)
      t_command = re.sub('\\\\SELF'  , t_job_id, t_command)

      for i_data in range(0, len(NGS_sample_data[ t_sample_id ])):
        t_data = NGS_sample_data[ t_sample_id ][i_data]
        t_re = '\\\\DATA\.' + str(i_data)
        t_command = re.sub(t_re, t_data, t_command)

      t_injobs = []
      if 'injobs' in t_job.keys():
        t_injobs = t_job[ 'injobs' ]
        for i_data in range(0, len(t_job[ 'injobs' ])):
          t_data = t_job[ 'injobs' ][i_data]
          t_re = '\\\\INJOBS\.' + str(i_data)
          t_command = re.sub(t_re, t_data, t_command)

      for i_data in range(0, len(CMD_opts)):
        t_data = CMD_opts[i_data]
        t_re = '\\\\CMDOPTS\.' + str(i_data)
        t_command = re.sub(t_re, t_data, t_command)

      v_command = ''
      if 'non_zero_files' in t_job.keys():
        for t_data in t_job[ 'non_zero_files' ]:
          v_command = v_command + \
            'if ! [ -s {0}/{1} ]; then echo "zero size {2}/{3}"; exit; fi\n'.format(t_job_id, t_data, t_job_id, t_data)

      print '-' * 80
      print t_sample_id
      print t_command
      print v_command
    
      f_start    = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.start.date'
      f_complete = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.complete.date'
      f_cpu      = pwd + '/' + t_sample_id + '/' + t_job_id + '/WF.cpu'
      t_sh_file  = '{0}/WF-sh/{1}.{2}.sh'.format(pwd, t_job_id, t_sample_id)
      t_infiles = []
      if 'infiles' in t_job.keys():
        t_infiles = map(lambda x: t_sample_id + "/" + x, t_job[ 'infiles' ])
      job_list[ t_job_id ][ t_sample_id ] = {
        'sample_id'    : t_sample_id,
        'job_id'       : t_job_id,
        'status'       : 'wait',       #### status can be wait (input not ready), ready (input ready), submitted (submitted or running), completed
        'command'      : t_command,
        'sh_file'      : t_sh_file, 
        'infiles'      : t_infiles,
        'injobs'       : t_injobs,
        'start_file'   : f_start,
        'complete_file': f_complete,
        'cpu_file'     : f_cpu }

      if not os.path.exists( t_sh_file ):
        try:
          tsh = open(t_sh_file, 'w')
          tsh.write('''{0}
{1}

my_host=`hostname`
my_pid=$$
my_core={2}
my_queue={3}
my_time_start=`date +%s`

cd {4}/{5}
mkdir {6}
if ! [ -f {7} ]; then date +%s > {7};  fi
{8}
{9}
date +%s > {10}
my_time_end=`date +%s`;
my_time_spent=$((my_time_end-my_time_start))
echo "sample={5} job={6} host=$my_host pid=$my_pid queue=$my_queue cores=$my_core time_start=$my_time_start time_end=$my_time_end time_spent=$my_time_spent" >> {11}

'''.format(t_execution['template'], pe_parameter, t_job['cores_per_cmd'], t_job['execution'], pwd, t_sample_id, t_job_id, f_start, t_command, v_command, f_complete, f_cpu ))
          tsh.close()
        except IOError:
          print 'cannot write to', job_list[ 't_job_id' ][ 't_sample_id' ][ 'sh_file' ]
          exit(1)
### END def make_job_list(NGS_config):



def task_log_cpu():
  '''
sub task_log_cpu {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id);

  my %cpu_info;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {

      $cpu_info{$t_job_id}{$t_sample_id} = [$t_wall, $t_cpu];
    }
  }

  foreach $t_sample_id (@NGS_samples) {
    my $f_cpu = "$pwd/$t_sample_id/WF.cpu";
    open(CPUOUT, "> $f_cpu") || die "Can not open $f_cpu";
    print CPUOUT "#job_name\tCores\tWall(s)\tWall_time\tCPU(s)\tCPU_time\n";
    my $min_start = 1402092131 * 999999;
    my $max_end   = 0;
    my $sum_cpu   = 0;
    foreach $t_job_id (keys %NGS_batch_jobs) {
      if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
      my $t_job = $NGS_batch_jobs{$t_job_id};
      my $t_core     = $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"};

      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $f_start    = $t_sample_job->{'start_file'};
      my $f_complete = $t_sample_job->{'complete_file'};
      my $f_cpu      = $t_sample_job->{'cpu_file'};
      my $t_start    = `cat $f_start`;    $t_start =~ s/\s//g; $min_start = $t_start if ($t_start < $min_start);
      my $t_end      = `cat $f_complete`; $t_end   =~ s/\s//g; $max_end   = $t_end   if ($t_end   > $max_end);
      my $t_wall     = int($t_end - $t_start);
         $t_wall     = 0 unless ($t_wall>0);

      my $t_cpu = 0;
      if (open(TCPU, $f_cpu)) {
        while($ll = <TCPU>) {
          chop($ll);
          if ($ll =~ /^(\d+)m(\d+)/) {
            $t_cpu += $1 * 60;
          }
        }
        close(TCPU);
      }
      $sum_cpu += $t_cpu;

      my $t_walls = time_str1($t_wall);
      my $t_cpus  = time_str1($t_cpu);
      print CPUOUT "$t_job_id\t$t_core\t$t_wall\t$t_walls\t$t_cpu\t$t_cpus\n";
    }
    my $t_wall = ($max_end - $min_start); $t_wall     = 0 unless ($t_wall>0);
    my $t_walls = time_str1($t_wall);
    my $sum_cpus= time_str1($sum_cpu);
    print CPUOUT "total\t-\t$t_wall\t$t_walls\t$sum_cpu\t$sum_cpus\n";
    close(CPUOUT);
  }
}
######### END task_log_cpu
  '''
#### END def task_log_cpu():


def task_list_jobs():
  '''
sub task_list_jobs {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);
  foreach $t_job_id (@NGS_batch_jobs) {
    $t_job = $NGS_batch_jobs{$t_job_id};
    #my @t_infiles = @{$t_job->{"infiles"}};
    my @t_injobs  = @{$t_job->{"injobs"}};

    #print "\tInput_files:", join(",", @t_infiles) if @t_infiles;
    print "$t_job_id\tIn_jobs:[" , join(",", @t_injobs), "]\tJob_level:$t_job->{'job_level'}\n";
  }
}
########## END task_list_jobs
  '''
#### END def task_list_jobs()

def task_snapshot():
  '''
sub task_snapshot {
  my ($t_job_id, $t_sample_id);
  my ($i, $j, $k);

  if ($this_task) {
    my $flag_qstat_xml_call = 0;
    foreach $t_job_id (keys %NGS_batch_jobs) {
      my $t_job = $NGS_batch_jobs{$t_job_id};
      my $t_execution = $NGS_executions{ $t_job->{"execution"} };
      my $exe_type = $t_execution->{type};
      $flag_qstat_xml_call = 1 if (($queue_system eq "SGE") and (($exe_type eq "qsub") or ($exe_type eq "qsub-pe")));
    }
    SGE_qstat_xml_query() if $flag_qstat_xml_call;

    foreach $t_sample_id (@NGS_samples) {
      foreach $t_job_id (keys %NGS_batch_jobs) {
        check_submitted_job($t_job_id, $t_sample_id);
      }
    }
  }

  my $max_len_sample = 0;
  foreach $t_sample_id (@NGS_samples) {
    $max_len_sample = length($t_sample_id) if (length($t_sample_id) > $max_len_sample);
  }
  my $max_len_job = 0;
  foreach $t_job_id (@NGS_batch_jobs) {
    $max_len_job = length($t_job_id) if (length($t_job_id) > $max_len_job);
  }

  print <<EOD;
Job status: 
.\twait
-\tsubmitted
r\trunning  
+\tcompleted
!\terror
EOD

  for ($i=$max_len_job-1; $i>=0; $i--) {
    print ' 'x$max_len_sample, "\t";
    foreach $t_job_id (@NGS_batch_jobs) {
      print " ", ($i<length($t_job_id) ? substr(reverse($t_job_id), $i, 1):" ");
    }
    print "\n";
  }

  foreach $t_sample_id (@NGS_samples) {
    print "$t_sample_id\t";
    foreach $t_job_id (@NGS_batch_jobs) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      if    ($status eq "completed") { print " +";}
      elsif ($status eq "submitted") { print " -";}
      elsif ($status eq "running"  ) { print " r";}
      elsif ($status eq "wait"     ) { print " .";}
      elsif ($status eq "error"    ) { print " !";}
      else                           { print " _";}
    }
    print "\n";
  }
}
########## END task_snapshot
  '''
### def task_snapshot():

def task_delete_jobs():
  '''
sub task_delete_jobs {
  my $opt = shift;
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id);
  my ($mode, $c) = split(/:/, $opt);
  my $tmp_sh = "NGS-$$.sh";

  open(TMPSH, "> $tmp_sh") || die "can not write to file $tmp_sh";
  print TMPSH "#Please execute the following commands\n";
  foreach $t_sample_id (@NGS_samples) {
    my %job_to_delete_ids = ();
    if ($mode eq "jobids") {
       %job_to_delete_ids = map {$_, 1} split(/,/,$c);
    }
    elsif ($mode eq "run_after") {
      die "file $c doesn't exist!" unless (-e $c);
      foreach $t_job_id (keys %NGS_batch_jobs) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        next unless (-e $t_sh_pid);   #### unless the job is submitted
        #$job_to_delete_ids{$t_job_id} = 1 if (file1_same_or_after_file2( $t_sample_job->{'start_file'} , $c));
        $job_to_delete_ids{$t_job_id} = 1 if (file1_same_or_after_file2( $t_sh_pid , $c));

      }
    }
    else {
      die "unknown option for deleting jobs: $opt";
    }

    # now %job_to_delete_ids are jobs need to be deleted
    # next find all jobs that depends on them, recrusively
    my $no_jobs_to_delete = scalar keys %job_to_delete_ids;
    while(1) {
      foreach $t_job_id (keys %NGS_batch_jobs) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        next unless (-e $t_sh_pid);   #### unless the job is submitted
        my @t_injobs  = @{ $t_sample_job->{'injobs'} };
        foreach my $t_job_id_2 (@t_injobs) {
          $job_to_delete_ids{$t_job_id} = 1 if ($job_to_delete_ids{$t_job_id_2});
        }
      }
      last if ($no_jobs_to_delete == (scalar keys %job_to_delete_ids)); #### no more depending jobs
      $no_jobs_to_delete = scalar keys %job_to_delete_ids;
    }

    if ($no_jobs_to_delete) {
      print TMPSH "#jobs to be deleted for $t_sample_id: ", join(",", keys %job_to_delete_ids), "\n";
      print       "#jobs to be deleted for $t_sample_id: ", join(",", keys %job_to_delete_ids), "\n";
      foreach $t_job_id (keys %job_to_delete_ids) {
        my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
        my $t_sh_file = $t_sample_job->{'sh_file'};
        my $t_sh_pid  = "$t_sh_file.pids";
        print TMPSH "\\rm -rf $pwd/$t_sample_id/$t_job_id\n";
        print TMPSH "\\rm $t_sh_pid\n";        
        print TMPSH "\\rm $t_sh_file.*.std*\n";

        #### find the qsub ids to be deleted 
        my $qids = `cat $t_sh_pid`; $qids =~ s/\n/ /g; $qids =~ s/\s+/ /g;
        print TMPSH "qdel $qids\n";
      }
    }
  }
  close(TMPSH);
  print "The script is not delete the file, please run $tmp_sh to delete files!!!\n\n";
}
  '''
#### END def task_delete_jobs()


def run_workflow(NGS_config):
  '''
  major look for workflow run
  '''
  queue_system = NGS_config.queue_system   #### default "SGE"
  sleep_time_min = 15
  sleep_time_max = 120
  sleep_time = sleep_time_min

  while 1:
    flag_job_done = True
    ########## reset execution_submitted to 0
    for i in NGS_config.NGS_executions.keys():
      execution_submitted[ i ] = False

    flag_qstat_xml_call = False
    for t_job_id in NGS_config.NGS_batch_jobs.keys():
      t_job = NGS_config.NGS_batch_jobs[t_job_id]
      t_execution = NGS_config.NGS_executions[ t_job['execution']]
      exe_type = t_execution['type']
      if (queue_system == SGE) and (exe_type in ['qsub','qsub-pe']):
        flag_qstat_xml_call = True

    if flag_qstat_xml_call:
      SGE_qstat_xml_query()

  #### END while 1:

 
  '''
  ########## check and update job status for submitted jobs
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};

      next if ($status eq "completed");
      ########## check file system to update job status
      ########## in case this is a restart run
      check_submitted_job($t_job_id, $t_sample_id);
      next if ($t_sample_job->{'status'} eq "completed");
      $flag_job_done = 0;
    }
  }

  if ($flag_job_done) { write_log("job completed!"); last; }

  ########## check and update job status based on dependance 
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};

      next unless ($status eq "wait");
      my @t_infiles = @{ $t_sample_job->{'infiles'} };
      my @t_injobs  = @{ $t_sample_job->{'injobs'} };
      my $t_ready_flag = 1;

      foreach $i (@t_infiles) {
        next if (-s $i); ####  non-zero size file
        $t_ready_flag = 0;
        last;
      }

      foreach $i (@t_injobs) {
        next if ( $job_list{$i}{$t_sample_id}->{'status'} eq "completed"); #### injob completed
        $t_ready_flag = 0;
        last;
      }
      if ($t_ready_flag) {
        $t_sample_job->{"status"} = "ready";
        write_log("$t_job_id,$t_sample_id: change status to ready");
      }
    }
  }

  ########## submit local sh jobs
  my $has_submitted_some_jobs = 0;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    next unless ($t_execution->{'type'} eq "sh");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"cores_per_node"} ); #### all cores are used

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next if ( ($execution_submitted{$t_execution_id} + $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"}) > $t_execution->{"cores_per_node"} ); #### no enough available cores
      #### now submitting 

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_pid  = "$t_sh_file.pids";
      for ($i=0; $i<$t_job->{"no_parallel"}; $i++) {
        $cmd = `sh $t_sh_file >/dev/null 2>&1 &`;
      }
      $cmd = `touch $t_sh_pid`;
      $t_sample_job->{'status'} = "submitted";
      write_log("$t_job_id,$t_sample_id: change status to submitted");
      $execution_submitted{ $t_execution_id } += $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"}; 
      $has_submitted_some_jobs = 1;
    }
  }

  ########## submit qsub-pe jobs, multiple jobs may share same node
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    next unless ($t_execution->{'type'} eq "qsub-pe");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"number_nodes"} ); #### resource full

    my $t_cores_per_node = $t_execution->{"cores_per_node"};
    my $t_cores_per_cmd  = $t_job->{"cores_per_cmd"};
    my $t_cores_per_job  = $t_cores_per_cmd * $t_job->{"no_parallel"};
    my $t_nodes_per_job  = $t_cores_per_job / $t_cores_per_node;

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_pid  = "$t_sh_file.pids";
      open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";

      for ($i=0; $i<$t_job->{"no_parallel"}; $i++) {
        my $t_stderr = "$t_sh_file.$i.stderr";
        my $t_stdout = "$t_sh_file.$i.stdout";
        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_file 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}
        print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
        $execution_submitted{$t_execution_id} += $t_nodes_per_job;
        write_log("$t_sh_bundle submitted for sample $t_sample_id, qsubid $cmd");
      }

      close(TID);
      $has_submitted_some_jobs = 1;
      $t_sample_job->{'status'} = "submitted";
    }
  } ########## END foreach $t_job_id (keys %NGS_batch_jobs) 

  ########## submit qsub jobs
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    my $t_job = $NGS_batch_jobs{$t_job_id};
    my $t_execution = $NGS_executions{ $t_job->{"execution"} };
    my $t_execution_id = $t_job->{"execution"};

    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    next unless ($t_execution->{'type'} eq "qsub");
    next if ( $execution_submitted{$t_execution_id} >= $t_execution->{"number_nodes"} ); #### resource full

    my $t_cores_per_node = $t_execution->{"cores_per_node"};
    my $t_cores_per_cmd  = $t_job->{"cores_per_cmd"};
    my $t_cores_per_job  = $t_cores_per_cmd * $t_job->{"no_parallel"};
    my $t_nodes_per_job  = POSIX::ceil($t_cores_per_job / $t_cores_per_node);
    my $t_cmds_per_node  = int($t_cores_per_node / $t_cores_per_cmd);
    my $t_jobs_per_node  = int($t_cores_per_node / $t_cores_per_job);

    ########## 1. this loop process jobs need 1 or more nodes per sample, ie. bundle within a sample, e.g. blast against refseq
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node <= 1);                 #### unless need >= 1 node, including jobs use between (51%-100%) cores per node
      last if ( ($execution_submitted{$t_execution_id} + $t_nodes_per_job) > $t_execution->{"number_nodes"}); #### no enough available queues

      my $t_sh_file = $t_sample_job->{'sh_file'};
      my $t_sh_bundle = "$sh_bundle_dir/$t_job_id.$t_sample_id.$$.sh";
      my $t_stderr    = "$t_sh_bundle.stderr";
      my $t_stdout    = "$t_sh_bundle.stdout";
      my $t_sh_pid  = "$t_sh_file.pids";

      open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";
      open(BSH, "> $t_sh_bundle") || die "can not write to $t_sh_bundle";
      print BSH <<EOD;
$t_execution->{"template"}
cd $pwd
EOD
      for ($i=0; $i<$t_cmds_per_node; $i++) {
        print BSH "sh $t_sh_file &\n";
        print BSH "sleep 3\n";
      }
      print BSH "wait\n";
      close(BSH);

      for ($i=0; $i<$t_nodes_per_job; $i++) {
        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_bundle 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}
        print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
        $execution_submitted{$t_execution_id}++;
        write_log("$t_sh_bundle submitted for sample $t_sample_id, qsubid $cmd");
      }
      close(TID);
      $has_submitted_some_jobs = 1;
      $t_sample_job->{'status'} = "submitted";
    } ########## END foreach $t_sample_id (@NGS_samples) 


    ########## 2. this loop process jobs need less than 1 node per sample, ie. bundle jobs across samples, e.g. qc 
    my @t_bundle = ();
    my $available_nodes = $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id};
    my $no_sample_can_be_processed = $available_nodes * $t_jobs_per_node;
    my @t_samples = ();
    my $t_batch_no = 0;

    foreach $t_sample_id (@NGS_samples) { #### same loop as next, to find out @t_samples and last sample can run
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node > 1);              #### unless a node can host 2 or more jobs
      last if ( $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id} <=0);
      push(@t_samples, $t_sample_id);
    }
    my $last_sample_can_run = $t_samples[-1];
    @t_samples = ();

    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      next unless ($status eq "ready");
      next unless ($t_jobs_per_node > 1);              #### unless a node can host 2 or more jobs
      last if ( $t_execution->{"number_nodes"} - $execution_submitted{$t_execution_id} <=0);
      push(@t_samples, $t_sample_id);

      #### bundle @t_samples to one qsub job
      if ((($#t_samples+1) == $t_jobs_per_node) or ($t_sample_id eq $last_sample_can_run)) {
        my $t_sh_bundle = "$sh_bundle_dir/$t_job_id.samples-$t_batch_no.$$.sh";
        my $t_stderr    = "$t_sh_bundle.stderr";
        my $t_stdout    = "$t_sh_bundle.stdout";

        open(BSH, "> $t_sh_bundle") || die "can not write to $t_sh_bundle";
      print BSH <<EOD;
$t_execution->{"template"}
cd $pwd
EOD
        foreach $i (@t_samples) {
          my $t_sh_file = $job_list{$t_job_id}{$i}->{'sh_file'};
          for ($j=0; $j<$t_job->{"no_parallel"}; $j++) {
            print BSH "sh $t_sh_file &\n";
            print BSH "sleep 3\n";
          }
        }
        print BSH "wait\n";
        close(BSH);

        $cmd = `qsub $t_execution->{"command_name_opt"} $t_job_id $t_execution->{"command_err_opt"} $t_stderr $t_execution->{"command_out_opt"} $t_stdout $t_sh_bundle 2>$log_fileq`;
        my $qsub_id = 0;
        if ($cmd =~ /(\d+)/) { $qsub_id = $1;} else {die "can not submit qsub job and return a id\n";}

        foreach $i (@t_samples) {
          my $t_sh_file = $job_list{$t_job_id}{$i}->{'sh_file'};
          my $t_sh_pid  = "$t_sh_file.pids";
          open(TID, "> $t_sh_pid")    || die "can not write to $t_sh_pid";
          print TID "$qsub_id\n"; #### $cmd returns qsub id, write these ids to pid file for future qstat 
          write_log("$t_sh_bundle submitted for sample $i, qsubid $cmd");
          close(TID);
          $job_list{$t_job_id}{$i}->{'status'} = "submitted";
        }

        $has_submitted_some_jobs = 1;
        $execution_submitted{$t_execution_id}++;
        @t_samples = (); #### clear
        $t_batch_no++;
      }
    } ########## END foreach $t_sample_id (@NGS_samples)
  } ########## END foreach $t_job_id (keys %NGS_batch_jobs) 


  #### if has submitted some jobs, reset waiting time, otherwise double waiting time
  print_job_status_summary();
  if ($has_submitted_some_jobs) {
    $sleep_time = $sleep_time_min;
  }
  else {
    $sleep_time = $sleep_time*2;
    $sleep_time = $sleep_time_max if ($sleep_time > $sleep_time_max);
  }
  write_log("sleep $sleep_time seconds");
  sleep($sleep_time);
  '''
#### END def run_workflow(NGS_config)

############################################################################################
# _______    ________  _________       ___________________   ________  .____       _________
# \      \  /  _____/ /   _____/       \__    ___/\_____  \  \_____  \ |    |     /   _____/
# /   |   \/   \  ___ \_____  \   ______ |    |    /   |   \  /   |   \|    |     \_____  \ 
#/    |    \    \_\  \/        \ /_____/ |    |   /    |    \/    |    \    |___  /        \
#\____|__  /\______  /_______  /         |____|   \_______  /\_______  /_______ \/_______  /
#        \/        \/        \/                           \/         \/        \/        \/ 
############################################################################################

if __name__ == "__main__":
  parser = argparse.ArgumentParser(formatter_class = RawTextHelpFormatter,
                                   description     = textwrap.dedent('''\

            ==================================================================
            Workflow tools for next generation genomics, metagenomics, RNA-seq
            and other type of omics data analyiss,
        
            Software originally developed since 2010 by Weizhong Li at UCSD
                                                          currently at JCVI
        
            http://weizhongli-lab.org/ngomicswf           liwz@sdsc.edu
            ==================================================================

   '''))

  parser.add_argument('-i', '--input',       help="workflow configration file, required", required=True)
  parser.add_argument('-s', '--sample_file', help='''
sample data file, required unless -S is present
File format example:
#Sample data file example, TAB or space delimited for following lines
Sample_ID1 sample_data_0 sample_data_1
Sample_ID2 sample_data_0 sample_data_1
Sample_ID3 sample_data_0 sample_data_1
  ''')
  parser.add_argument('-S', '--sample_name', help='''
sample data from command line, required unless -s is present
format:
Sample_ID1:sample_data_0:sample_data_0:sample_data_1,Sample_ID2:sample_data_0:sample_data_1
  ''')
  parser.add_argument('-t', '--parameter_file', help='''
replace default paramters in workflow configration file
File format example:
#parameter file example, TAB or space delimited for following lines
CMDOPT JobID_A:opt0:opt1:opt2
CMDOPT JobID_B:opt0:opt1
  ''')
  parser.add_argument('-T', '--parameter_name', help='''
parameter from command line
format:
JobID_A:opt0:opt1:opt2,JobID_B:opt0:opt1
  ''')
  parser.add_argument('-j', '--jobs', help='''run sub set of jobs, optional
the workflow will run all jobs by default.
to run sub set of jobs: -j qc or -j qc,fastqc
  ''')
  parser.add_argument('-J', '--task', help='''optional tasks
write-sh: write sh files and quite
log-cpu: gathering cpu time for each run for each sample
list-jobs: list jobs
snapshot: snapshot current job status
delete-jobs: delete jobs, must supply jobs delete syntax by option -Z
  e.g. -J delete-jobs -Z jobids:assembly,blast  ---delete assembly,blast and all jobs depends on them
       -J delete-jobs -Z run_after:filename     ---delete jobs that has start time (WF.start.date) after this file, and all depending jobs
  ''')
  parser.add_argument('-Z', '--second_parameter', help='secondary parameter used by other options, such as -J')
  parser.add_argument('-Q', '--queye', help='queue system, e.g. PBS, SGE', default='SGE')

  args = parser.parse_args()

  if (args.sample_file is None) and (args.sample_name is None) :
    parser.error('No sample file or sample name')

  NGS_config = imp.load_source('NGS_config', args.input)

  read_samples(args)
  print 'Samples'
  print NGS_samples
  print NGS_sample_data

  read_parameters(args)
  print 'Parameters'
  print NGS_opts

  if args.jobs:
    subset_flag = True
    subset_jobs = re.split(',', args.jobs)
    ## -- add_subset_jobs_by_dependency()
    print subset_jobs

  if not os.path.exists('WF-sh'):
    os.system('mkdir WF-sh')

  ## -- task_level_jobs();
  ## -- my @NGS_batch_jobs = sort {($NGS_batch_jobs{$a}->{'job_level'} <=> $NGS_batch_jobs{$b}->{'job_level'}) or ($a cmp $b)} keys %NGS_batch_jobs;

  make_job_list(NGS_config)

  ## single task
  if args.task:
    if args.task == 'log-cpu':
      task_log_cpu()
      exit(0)
    elif args.task == 'list-jobs':
      task_list_jobs()
      exit(0)
    elif args.task == 'snapshot':
      task_snapshot()
      exit(0)
    elif args.task == 'delete-jobs':
      task_delete_jobs(args.second_parameter)
      exit(0)
    elif args.task == 'write-sh':
      exit(0)
    else:
      print 'undefined task' + args.task
      exit(1)

################################################################################################
#  _____               _   _  _____  _____  _           _       _           _       _         
# |  __ \             | \ | |/ ____|/ ____|| |         | |     | |         (_)     | |        
# | |__) |   _ _ __   |  \| | |  __| (___  | |__   __ _| |_ ___| |__        _  ___ | |__  ___ 
# |  _  / | | | '_ \  | . ` | | |_ |\___ \ | '_ \ / _` | __/ __| '_ \      | |/ _ \| '_ \/ __|
# | | \ \ |_| | | | | | |\  | |__| |____) || |_) | (_| | || (__| | | |     | | (_) | |_) \__ \
# |_|  \_\__,_|_| |_| |_| \_|\_____|_____/ |_.__/ \__,_|\__\___|_| |_|     | |\___/|_.__/|___/
#                                      ______                      ______ _/ |                
#                                     |______|                    |______|__/                 
########## Run NGS_batch_jobs for each samples http://patorjk.com/software/taag
################################################################################################

  run_workflow(NGS_config)
  task_log_cpu()


"""
sub write_log {
  my @txt = @_;
  my $i;
  my $date = `date`; chop($date);
  foreach $i (@txt) {
    print LOG    "$date $i\n";
    print STDERR "$date $i\n";
  }
  print LOG    "\n";
  print STDERR "\n";
}
########## END write_log

sub SGE_qstat_xml_query {
  my ($i, $j, $k, $cmd, $ll);
  %qstat_xml_data = (); #### global
  $cmd = `qstat -f -xml`;
  if ($cmd =~ /<queue_info/) { #### dummy 
    $qstat_xml_data{"NULL"}= ["NULL","NULL"];
  }

  my @lls = split(/\n/, $cmd);
  $i = 2; #### skip first 2 lines
  for (;     $i<$#lls+1; $i++) {
    if ($lls[$i] =~ /<job_list/) {
      my ($id, $name, $state);
      for (; $i<$#lls+1; $i++) {
        last if ($lls[$i] =~ /<\/job_list/);
        if ($lls[$i] =~ /<JB_job_number>(\d+)/) {  $id = $1;}
        if ($lls[$i] =~ /<JB_name>([^<]+)/) { $name = $1;}
        if ($lls[$i] =~ /<state>([^<]+)/) {$state = $1;}
      }
      if (defined($id) and defined($name) and defined($state)) {
        $qstat_xml_data{$id} = [$name, $state];
      }
    }
  }
}

########## check submitted job by checking pids, or qsub ids
########## update job status from wait|ready -> submitted if pid file exit (in case of restart of this script)
########## update job status from wait|ready|submitted -> completed if sh calls or qsub calls finished
##########    these pids or qsub ids are done
sub check_submitted_job {
  my ($t_job_id, $t_sample_id) = @_;
  my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
  my $t_job = $NGS_batch_jobs{$t_job_id};
  my $t_execution = $NGS_executions{ $t_job->{"execution"} };

  my ($i, $j, $k, $flag, $ll, $cmd);

  my $t_sh_file = $t_sample_job->{'sh_file'};
  my $t_sh_pid  = "$t_sh_file.pids";

  # status won't change unless there is a pid file
  return unless (-e $t_sh_pid);

  my $status = $t_sample_job->{'status'};
  if (($status eq "wait") or ($status eq "ready")) {
    $t_sample_job->{'status'} = "submitted";
    write_log("$t_job_id,$t_sample_id: change status to submitted");
  }

  my $exe_type = $t_execution->{type};

  if ($exe_type eq "sh") {
    $cmd = `ps -ef | grep "$t_sh_file" | grep -v grep`;
    if ($cmd =~ /\w/) { # still running 
      $execution_submitted{ $t_job->{"execution"} } += $t_job->{"cores_per_cmd"} * $t_job->{"no_parallel"};
    }
    elsif (validate_job_files($t_job_id, $t_sample_id)) {
      $t_sample_job->{'status'} = "completed";
      write_log("$t_job_id,$t_sample_id: change status to completed");
    }
    else {
      $t_sample_job->{'status'} = "error";
      write_log("$t_job_id,$t_sample_id: change status to error");
    }
    return;
  }
  elsif (($exe_type eq "qsub") or ($exe_type eq "qsub-pe")) {
    my @pids = ();
    open(CHECK, $t_sh_pid) || die "Can not open $t_sh_pid\n";
    while($ll = <CHECK>) {
      chop($ll); next unless ($ll =~ /\w/);
      push(@pids, $ll);
    }
    close(CHECK);

    my $finish_flag = 1;
    foreach $i (@pids) {
      if (($queue_system eq "SGE") and %qstat_xml_data) {
        if (defined($qstat_xml_data{$i})) {
          $t_sample_job->{'status'} = "running" if (($qstat_xml_data{$i}->[1] eq "r") and ($t_sample_job->{'status'} eq "submitted"));
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
      elsif ($queue_system eq "SGE") {
        $cmd = `qstat -j $i | grep job_number`;
        if ($cmd =~ /$i/) {
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
      else {
        $cmd = `qstat -r $i | grep $i`;
        $j = (split(/\D/,$cmd))[0];
        if ($j == $i) { # this job is running
          $finish_flag = 0;
          $execution_submitted{ $t_job->{"execution"} } ++; 
        }
      }
    }
    if ($finish_flag == 1) {
      if (validate_job_files($t_job_id, $t_sample_id)) {
        $t_sample_job->{'status'} = "completed";
        write_log("$t_job_id,$t_sample_id: change status to completed");
      }
      else {
        $t_sample_job->{'status'} = "error";
        write_log("$t_job_id,$t_sample_id: change status to error");
      }
    }
    return;
  }
  else {
    die "unknown execution type: $exe_type\n";
  }
}
########## END sub check_submitted_job 


# WF.start.date and WF.complete.date need to have non-zero size
sub validate_job_files {
  my ($t_job_id, $t_sample_id) = @_;
  my ($i, $j, $k);
  my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};

  return 0 unless (-s $t_sample_job->{'start_file'} );
  return 0 unless (-s $t_sample_job->{'complete_file'} );
  return 0 unless (-s $t_sample_job->{'cpu_file'} );

  return 1; #### pass
}
########## END validate_job_files


sub print_job_status_summary {
  my ($t_job_id, $t_sample_id);
  my ($i, $j, $k);

  my %job_status = ();
  my $job_total = 0;
  foreach $t_job_id (keys %NGS_batch_jobs) {
    if ($subset_flag) {next unless ($subset_jobs{$t_job_id});} 
    foreach $t_sample_id (@NGS_samples) {
      my $t_sample_job = $job_list{$t_job_id}{$t_sample_id};
      my $status = $t_sample_job->{'status'};
      $job_status{$status}++;
      $job_total++;
    }
  }

  print STDERR "total jobs: $job_total,";
  foreach $i (sort keys %job_status) {
    print STDERR "$i: $job_status{$i},";
  } 
  print STDERR "\n"; 
}
########## END print_job_status_summary


sub add_subset_jobs_by_dependency {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);

  while(1) {
    my $num_subset_jobs = scalar keys %subset_jobs;

    foreach $t_job_id (keys %subset_jobs) {
      $t_job = $NGS_batch_jobs{$t_job_id};
      my @t_injobs  = @{$t_job->{"injobs"}};

      for $j (@t_injobs) {
        $subset_jobs{$j} = 1;
      }
    }

    last if ($num_subset_jobs == scalar keys %subset_jobs);
  }
}
########## END add_subset_jobs_by_dependency


sub task_level_jobs {
  my ($i, $j, $k, $ll, $t_job_id, $t_sample_id, $t_job);
  my %job_level = ();

  while(1) {
    my $change_flag = 0;

    foreach $t_job_id (keys %NGS_batch_jobs) {
      $t_job = $NGS_batch_jobs{$t_job_id};
      my @t_injobs  = @{$t_job->{"injobs"}};

      if (@t_injobs) {
        my $max_level_injob;
        foreach $j (@t_injobs) {
          next unless defined ($job_level{$j});
          $max_level_injob = $job_level{$j} if ($job_level{$j} > $max_level_injob);          
        }

        next unless (defined($max_level_injob));
        $max_level_injob++; #### one more level 
        if (not defined ($job_level{$t_job_id})) {
          $job_level{$t_job_id}=$max_level_injob;
          $change_flag = 1;
        }
        elsif ($max_level_injob > $job_level{$t_job_id}) {
          $job_level{$t_job_id}=$max_level_injob;
          $change_flag = 1;
        }
      }
      else {
        if (not defined ($job_level{$t_job_id})) {
          $job_level{$t_job_id}=1;
          $change_flag = 1;
        }
      }
    }
    last unless ($change_flag);
  }

  foreach $t_job_id (sort keys %NGS_batch_jobs) {
    $NGS_batch_jobs{$t_job_id}->{"job_level"} = $job_level{$t_job_id};
  }
}


sub file1_after_file2 {
  my ($file1, $file2) = @_;

  # if not exist file1, assume it is in future, so it is newer
  if (not -e ($file1)) {return 0;}
  if (not -e ($file2)) {return 0;}

  my $mtime1 = (stat($file1))[9];
  my $mtime2 = (stat($file2))[9];

  return ( ($mtime1 > $mtime2) ? 1 : 0);
}
######## END file1_after_file2

sub file1_same_or_after_file2 {
  my ($file1, $file2) = @_;

  # if not exist file1, assume it is in future, so it is newer
  if (not -e ($file1)) {return 0;}
  if (not -e ($file2)) {return 0;}

  my $mtime1 = (stat($file1))[9];
  my $mtime2 = (stat($file2))[9];

  return ( ($mtime1 >= $mtime2) ? 1 : 0);
}
######## END file1_after_file2



sub time_str1 {
  my $s = shift;
  my $str = "";

  $str .= int($s/3600); $str .= "h"; $s = $s % 3600;
  $str .= int($s/60);   $str .= "m"; $s = $s % 60;
  $str .= $s;           $str .= "s";

  return $str;
}
########## END time_str1;


"""

