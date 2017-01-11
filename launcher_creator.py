#!/usr/bin/env python

# Create TACC launchers for Stampede and Lonestar, and just run commands otherwise.
# Depending on the host, launcher_creator.py generated a SLURM or SGE launcher as appropriate.
# 
# Copyright 2013, Regents Of The University Of Texas At Austin
# Originally written by Scott Hunicke-Smith
# New version with SLURM support by Benjamin Goetz

try:
    import argparse
    import subprocess
except:
    print 'Try typing "module load python" and then running this again.'
    import sys
    sys.exit(1)
    
import os
import sys
import re

lonestar_parameters='''#!/bin/bash
#
# Simple SGE script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#         <------ Setup Parameters ------>
#
#$ -N {name}
#$ -pe {wayness}way {total_threads}
#$ -q {queue}
#$ -o {name}.o$JOB_ID
#$ -l h_rt={time}
#$ -V
{email_line}{hold_line}#$ -cwd
{allocation_line}
'''

stampede_parameters='''#!/bin/bash
#SBATCH -J {name}
#SBATCH -n {num_jobs}
{num_nodes_line}#SBATCH -p {queue}
#SBATCH -o {name}.o%j
#SBATCH -e {name}.e%j
#SBATCH -t {time}
{allocation_line}
{email_line}'''

bash_and_modules='''
module load launcher
{modules}

{bash_commands}
'''

parametric_job_submission='''
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher 
export CONTROL_FILE={control_file}
export WORKDIR=.
# 
# Variable description:
#
#  EXECUTABLE     = full path to the job launcher executable
#  CONTROL_FILE   = text input file which specifies
#                   executable for each process
#                   (should be located in WORKDIR)
#  WORKDIR        = location of working directory
#
#      <------ End Setup Parameters ------>
#--------------------------------------------------------
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if [ ! -e $WORKDIR ]; then
    echo " "
    echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -f $EXECUTABLE ]; then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -f $WORKDIR/$CONTROL_FILE ]; then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
fi


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
'''

parametric_job_submission_ls5='''

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE={jobfile}

$LAUNCHER_DIR/paramrun
'''

def file_len(fname):
    f = open(fname)
    contents = f.readlines()
	# Only count non-empty lines and count the
	# last line even if it has no carriage return
    i = 0
    for line in contents:
        stripped_line = line.strip()
        #comment lines
        if re.match('^\s*#', stripped_line):
        	continue
        if not stripped_line:
            continue
        else:
            i+=1
    f.close()
    return i


def generate_lonestar_launcher(results, launcher_file):
    """
    Write a Lonestar launcher SGE script.
    """

    # Must account (on Lonestar) for normal having 12 threads per node;
    if results.queue == 'largemem':
        threads_per_node = 24
    else:
        threads_per_node = 12
    
    if not results.wayness:
        results.wayness = 12
    
    if results.allocation:
        allocation_line = '#$ -A {}'.format(results.allocation)
    else:
        allocation_line = ''
    
    if results.email:
        email_line = '#$ -M {}\n#$ -m be\n'.format(results.email)
    else:
        email_line = ''

    if results.hold and not results.submit_cmd:
        hold_line='#$ -hold_jid {}\n'.format(results.hold)
    else:
        hold_line=''

    if not results.num_nodes:
        if results.job:
            num_nodes = file_len(results.job)
        else:
            num_nodes = 1

        if num_nodes == 0:
            print 'Your job file appears to be empty. Please fix this error.'
            return

        while num_nodes % results.wayness != 0:
            num_nodes += 1
        total_threads = num_nodes * threads_per_node / results.wayness
    else:
        num_nodes = results.num_nodes
        total_threads = results.num_nodes * threads_per_node

    if not results.submit_cmd:
        print 'Job file has %i lines.' % num_nodes

    if not results.submit_cmd:
        print 'Using {} nodes.'.format(num_nodes)
        print 'Writing to {}.'.format(results.launcher)
    
    launcher_file.write(lonestar_parameters.format(
        name=results.name,
        wayness=results.wayness,
        total_threads=total_threads,
        walrus_total_threads=total_threads,
        queue=results.queue,
        time=results.time,
        email_line=email_line,
        hold_line=hold_line,
        allocation_line=allocation_line))
    
    launcher_file.write(bash_and_modules.format(modules=results.modules,
               bash_commands=results.bash_commands))
               
    if results.job:
       launcher_file.write(parametric_job_submission.format(control_file=results.job))


def generate_stampede_launcher(results, launcher_file):
    """
    Write a Stampede launcher SLURM script.
    """

    # Must account (on Stampede) for normal having 16 threads per node
    if results.queue == 'largemem':
        threads_per_node = 32
    else:
        threads_per_node = 16

    num_nodes_line = ''
    # if not results.wayness:
    #     results.wayness = 16
    #     num_cores_line = ''

    if results.allocation:
        allocation_line = '#SBATCH -A {}'.format(results.allocation)
    else:
        allocation_line = ''

    if results.email:
        email_line = '#SBATCH --mail-type=ALL\n#SBATCH --mail-user={}\n'.format(results.email)
    else:
        email_line = ''

    if not results.num_nodes:
        # If a number of nodes is not explicitly specified, calculate based on number of commands in jobs file, and optional wayness.
        if results.job:
            num_jobs = file_len(results.job)
            if num_jobs == 0:
                print 'Your job file appears to be empty. Please fix this error.'
                return
        else:
            num_jobs = 1
    
        if results.wayness:
            num_nodes, remainder = divmod(num_jobs, results.wayness)
            if remainder != 0:
                num_nodes += 1
            num_nodes_line = '#SBATCH -N {}\n'.format(num_nodes)
        else:
            num_nodes, remainder = divmod(num_jobs, 16)
            if remainder != 0:
                num_nodes += 1                
    else:
        # A number of nodes has been specified. Default wayness is 16 on Stampede. Use specified wayness to calculate total tasks.
        num_nodes = results.num_nodes
        num_nodes_line = '#SBATCH -N {}\n'.format(num_nodes)
        
        if results.wayness:
            num_jobs = num_nodes * results.wayness
        else:
            num_jobs = num_nodes * 16

    if not results.submit_cmd:
        print 'Using {} nodes.'.format(num_nodes)
        print 'Writing to {}.'.format(results.launcher)

    launcher_file.write(stampede_parameters.format(
        name=results.name,
        num_jobs=num_jobs,
        num_nodes_line=num_nodes_line,
        queue=results.queue,
        time=results.time,
        email_line=email_line,
        allocation_line=allocation_line))

    launcher_file.write(bash_and_modules.format(modules=results.modules,
        bash_commands=results.bash_commands))

    if results.job:
        launcher_file.write(parametric_job_submission.format(control_file=results.job))


def generate_ls5_launcher(results, launcher_file):
    """
    Write a Lonestar 5 launcher SLURM script.
    """

    # Must account (on Lonestar 5) for normal having 16 threads per node
    if results.queue == 'largemem':
        threads_per_node = 32
    elif results.queue == 'hugemem':
        threads_per_node = 20
    else:
        threads_per_node = 24

    num_nodes_line = ''
    # if not results.wayness:
    #     results.wayness = 16
    #     num_cores_line = ''

    if results.allocation:
        allocation_line = '#SBATCH -A {}'.format(results.allocation)
    else:
        allocation_line = ''

    if results.email:
        email_line = '#SBATCH --mail-type=ALL\n#SBATCH --mail-user={}\n'.format(results.email)
    else:
        email_line = ''

    if not results.num_nodes:
        # If a number of nodes is not explicitly specified, calculate based on number of commands in jobs file, and optional wayness.
        if results.job:
            num_jobs = file_len(results.job)
            if num_jobs == 0:
                print 'Your job file appears to be empty. Please fix this error.'
                return
        else:
            num_jobs = 1
    
        if results.wayness:
            num_nodes, remainder = divmod(num_jobs, results.wayness)
            if remainder != 0:
                num_nodes += 1
            num_nodes_line = '#SBATCH -N {}\n'.format(num_nodes)
        else:
            num_nodes, remainder = divmod(num_jobs, 24)
            if remainder != 0:
                num_nodes += 1                
    else:
        # A number of nodes has been specified. Default wayness is 16 on Stampede. Use specified wayness to calculate total tasks.
        num_nodes = results.num_nodes
        num_nodes_line = '#SBATCH -N {}\n'.format(num_nodes)
        
        if results.wayness:
            num_jobs = num_nodes * results.wayness
        else:
            num_jobs = num_nodes * 24

    if not results.submit_cmd:
        print 'Using {} nodes.'.format(num_nodes)
        print 'Writing to {}.'.format(results.launcher)

    launcher_file.write(stampede_parameters.format(
        name=results.name,
        num_jobs=num_jobs,
        num_nodes_line=num_nodes_line,
        queue=results.queue,
        time=results.time,
        email_line=email_line,
        allocation_line=allocation_line))

    launcher_file.write(bash_and_modules.format(modules=results.modules,
        bash_commands=results.bash_commands))

    if results.job:
        launcher_file.write(parametric_job_submission_ls5.format(jobfile=results.job))


def check_time_and_queue(results):
    host = subprocess.check_output(["hostname", "-f"]).split('.')[1]
    queues = {'ls4':['normal','development','largemem','gpu','vis','serial'],
              'stampede':['normal','development','largemem','serial','large','request','normal-mic','gpu','gpudev','vis']}

    if host in queues and results.queue not in queues[host]:
        sys.exit( 'Your queue selection "' + results.queue + '" is not available.' + '\n' +
                  'Available queues: ' + ', '.join(queues[host]))

    if results.time == None:
        parser.print_help()
        sys.exit('You did not give a job time (-t hh:mm:ss).')
    
    time_match = re.match(r'(\d?\d):(\d\d):(\d\d)', results.time)
    if not time_match:
        sys.exit( 'Your time argument -t ' + results.time + ' is not in the right format.' + '\n' +  
                  'Please change it to match (-t hh:mm:ss).' )

    hh = int(time_match.group(1))
    mm = int(time_match.group(2))
    ss = int(time_match.group(3))

    if hh == 0 and mm == 0 and ss == 0: sys.exit( 'Please enter a non-zero time value.' )
    if mm > 59 or ss > 59: sys.exit( results.time + ' is not a valid time.' )

    if   host == "ls4"      and results.queue == 'development': time_limit = 1
    elif host == "ls4"      and results.queue == 'serial':      time_limit = 12
    elif host == "stampede" and results.queue == 'development': time_limit = 2
    elif host == "stampede" and \
        results.queue in ['gpudev','visdev']:                   time_limit = 4
    elif host == "stampede" and results.queue == 'vis':         time_limit = 8
    elif host == "stampede" and results.queue == 'serial':      time_limit = 12
    elif host == "stampede" and \
        results.queue in ['normal','largemem','normal-mic']:    time_limit = 48
    elif host == "ls5"      and results.queue == 'development': time_limit = 2
    elif host == "ls5"      and results.queue == 'gpu':         time_limit = 24
    elif host == "ls5"      and results.queue == 'vis':         time_limit = 8
    elif host == "ls5" and \
        results.queue in ['normal','largemem','hugemem']:       time_limit = 48
    else:                                                       time_limit = 24

    if (hh > time_limit) or (hh == time_limit and (mm > 0 or ss > 0)):
        sys.exit( 'Your time argument -t ' + results.time + ' exceeds the ' + 
                        str(time_limit) + ' hour limit for the ' + results.queue + ' queue' + '\n' +
                  'on this machine. Please change queue or reduce requested time.' )


class ArgumentParserDefaultHelp(argparse.ArgumentParser):
    '''Simple variation on ArgumentParser. When an error is thrown, this inherited class prints the entire help message.'''

    def error(self, message):
        self.print_help()
        sys.stderr.write('\nerror: {}\n\n'.format(message))
        sys.exit(2)


def main():
    # Get environment variables, os.environ.get returns 'none' if the environment variable isn't defined
    env_allocation = os.environ.get('ALLOCATION')
    env_email_address = os.environ.get('EMAIL_ADDRESS')
    
    parser = ArgumentParserDefaultHelp(description='''Create TACC launchers for Stampede and Lonestar, and just run commands otherwise.
Depending on the host, launcher_creator.py generates a SLURM or SGE launcher as appropriate.
Report problems to rt-other@ccbb.utexas.edu''')
    required = parser.add_argument_group('Required')
    required.add_argument('-n', action='store', dest='name', required=True, help='The name of your job.')
    required.add_argument('-t', action='store', dest='time', required=True, help='The time you want to give to your job. Format: hh:mm:ss')
    commands = parser.add_argument_group('Commands', 'You must use at least one of these options to submit your commands for TACC.')
    commands.add_argument('-j', action='store', dest='job', help='The name of the job file containing your commands.')
    commands.add_argument('-b', action='store', dest='bash_commands', default='', help='A string of Bash commands that are executed before the parametric jobs are launched.')
    optional = parser.add_argument_group('Optional')
    optional.add_argument('-q', action='store', dest='queue', default='development', help='The TACC allocation for job submission. Default="development"')
    optional.add_argument('-a', '-A', action='store', dest='allocation', nargs='?', default=env_allocation, help='The TACC allocation for job submission. You can set a default ALLOCATION environment variable.')
    optional.add_argument('-m', action='store', dest='modules', default='', help='A list of module commands. The "launcher" module is always automatically included. Example: -m "module swap intel gcc; module load bedtools"')
    optional.add_argument('-w', action='store', dest='wayness', type=int, help='Wayness: the number of commands you want to give each node. Default= 12 on Lonestar, 16 on Stampede.')
    optional.add_argument('-N', action='store', dest='num_nodes', type=int, help='Number of nodes to request. You probably don\'t need this option. You ONLY need it if you want to run a job list that isn\'t defined at the time you submit the launcher.')
    optional.add_argument('-e', action='store', dest='email', nargs='?', const=env_email_address, help='Your email address if you want to receive an email from Lonestar when your job starts and ends. Without an argument, it will use a default EMAIL_ADDRESS environment variable.')
    optional.add_argument('-l', action='store', dest='launcher', help='The name of the *.sge launcher script that will be created. Default="<name>.sge"')
    optional.add_argument('-s', dest='submit_cmd', action='store_true', help='Echoes the launcher filename to stdout.')
    deprecated = parser.add_argument_group('Deprecated', 'Rewrite any scripts that use this.')
    deprecated.add_argument('-H', action='store', dest='hold', help='Hold job in queue until the job specified in this option finishes.')
    
    results = parser.parse_args()

    # Handle allocation defaults
    # if not results.allocation:
    #     results.allocation = os.environ.get('ALLOCATION')
    
    # Check for required parameters
    # if results.name == None:
    #     print 'You did not give a job name (-n <name>).'
    #     parser.print_help()
    #     return

    # Check for valid time and queue entries
    check_time_and_queue(results)

    # Check for valid email address, if email is given
    if results.email:
        if not re.match(r'^.+@.+\..{2,4}$', results.email):
            print 'Not a valid email address, please check it.'
            return

    if not (results.job or results.bash_commands):
        if os.path.isfile('commands'):
            results.job = 'commands'
        else:
            print '-== Neither a job file nor bash commands given. You need at least one to do anything! ==-'
            return
    
    if not results.submit_cmd:
        print 'Project {}.'.format(results.name)
        print 'Using job file {}.'.format(results.job)
        print 'Using {} queue.'.format(results.queue)
        print 'For {} time.'.format(results.time)
        print 'Using {} allocation.'.format(results.allocation)
    
    # if (results.hold != None) and not results.submit_cmd:
    #     print 'Holding until {} finishes.'.format(results.hold)
    
    if not results.submit_cmd:
        if results.email != None:
            print 'Sending start/stop email to {}.'.format(results.email)
        else:
            print 'Not sending start/stop email.'
    
    
    host = subprocess.check_output(["hostname", "-f"]).split('.')[1]

    if not results.launcher:
        if host == "ls4":
            results.launcher = '{}.sge'.format(results.name)
        elif host in ["stampede", "ls5"]:
            results.launcher = '{}.slurm'.format(results.name)
        else:
            results.launcher

    
    launcher_file = open(results.launcher, 'w')
    
    if host == "ls4":
        generate_lonestar_launcher(results, launcher_file)
    elif host == "stampede":
        generate_stampede_launcher(results, launcher_file)
    elif host == "ls5":
        generate_ls5_launcher(results, launcher_file)
    else:
        generate_generic_launcher(results, launcher_file)

    launcher_file.close()

    if results.submit_cmd:
        print results.launcher
    else:
        if host in ["stampede", "ls5"]:
            print 'Launcher successfully created. Type "sbatch {}" to queue your job.'.format(results.launcher)
        else:
            print 'Launcher successfully created. Type "qsub {}" to queue your job.'.format(results.launcher)

if __name__ == '__main__':
    main()
