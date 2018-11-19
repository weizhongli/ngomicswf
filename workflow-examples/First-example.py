#!/usr/bin/python
################################################################################
# NGS workflow by Weizhong Li, http://weizhongli-lab.org
################################################################################

queue_system = 'SGE'

########## local variables etc. Please edit
ENV={
  'NGS_root' : '/home/myhome/bin',
}

########## computation resources for execution of jobs
NGS_executions = {}
NGS_executions['qsub_1'] = {
  'type'                : 'qsub-pe',
  'pe_para'             : '-pe orte', #### '-pe orte' work with orte environment with $pe_slot allocation rule, '-pe threade
  'qsub_exe'            : 'qsub',
  'cores_per_node'      : 32,
  'number_nodes'        : 64,
  'user'                : 'weizhong', #### I will use command such as qstat -u weizhong to query submitted jobs
  'command'             : 'qsub',
  'command_name_opt'    : '-N',
  'command_err_opt'     : '-e',
  'command_out_opt'     : '-o',
  'template'            : '''#!/bin/bash
#$ -q RNA.q
#$ -v PATH
#$ -V

'''
}

NGS_executions['sh_1'] = {
  'type'                : 'sh',
  'cores_per_node'      : 2,
  'number_nodes'        : 1,
  'template'            : '''#!/bin/bash

'''
}


NGS_batch_jobs = {}
NGS_batch_jobs['Job_A'] = {
  'CMD_opts'         : ['10'],
  'non_zero_files' : ['output_Job_A'],
  'execution'        : 'sh_1',               # where to execute
  'cores_per_cmd'    : 1,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

touch $SELF/output_Job_A
for i in `seq 1 $CMDOPTS.0`; do echo "$SAMPLE $DATA.0 line $i $SELF" >> $SELF/output_Job_A; done

'''
}

NGS_batch_jobs['Job_B'] = {
  'injobs'         : ['Job_A'],
  'CMD_opts'         : ['thread_1', 'thread_2'],
  'non_zero_files' : ['output_Job_B1',"output_Job_B2"],
  'execution'        : 'sh_1',               # where to execute
  'cores_per_cmd'    : 2,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

# you can have two threads in parallel
echo "$SAMPLE, thread 1 running"; cat $INJOBS.0/output_Job_A | sed "s/$/ $SELF $CMDOPTS.0/" > $SELF/output_Job_B1 &
echo "$SAMPLE, thread 1 running"; cat $INJOBS.0/output_Job_A | sed "s/$/ $SELF $CMDOPTS.1/" > $SELF/output_Job_B2 &

wait 
'''
}

NGS_batch_jobs['Job_C'] = {
  'injobs'         : ['Job_B'],
  'non_zero_files' : ['output_Job_C'],
  'execution'        : 'sh_1',               # where to execute
  'cores_per_cmd'    : 1,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

cat $INJOBS.0/output_Job_B* | sed "s/$/ $SELF/" > $SELF/output_Job_C
'''
}

NGS_batch_jobs['Job_D'] = {
  'injobs'         : ['Job_B'],
  'non_zero_files' : ['output_Job_D'],
  'execution'        : 'sh_1',               # where to execute
  'cores_per_cmd'    : 1,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

cat $INJOBS.0/output_Job_B* | sed "s/$/ $SELF/" > $SELF/output_Job_D
'''
}


NGS_batch_jobs['Job_E'] = {
  'injobs'         : ['Job_C','Job_D'],
  'non_zero_files' : ['output_Job_E'],
  'execution'        : 'sh_1',               # where to execute
  'cores_per_cmd'    : 1,                    # number of threads used by command below
  'no_parallel'      : 1,                    # number of total jobs to run using command below
  'command'          : '''

cat $INJOBS.0/output_Job_C $INJOBS.1/output_Job_D | sed "s/$/ $SELF/" > $SELF/output_Job_E
'''
}


