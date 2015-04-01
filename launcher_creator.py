#!/usr/bin/env python
try:
    import argparse
except:
    print 'Try typing "module load python" and then running this again.'
    import sys
    sys.exit(1)
    
import sys

def file_len(fname):
    f = open(fname)
    for i, l in enumerate(f):
        pass
    return i + 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', action='store', dest='email', default='yourname@utexas.edu', help='Your email address if you want to receive an email from Lonestar when your job starts and ends. Default=OFF')
    parser.add_argument('-q', action='store', dest='queue', default='development', help='The TACC allocation for job submission. Default="development"')
    parser.add_argument('-t', action='store', dest='time', default='1:00:00', help='The time you want to give to your job. Format: hh:mm:ss')
    parser.add_argument('-a', action='store', dest='allocation', default='ecogeno', help='The TACC allocation for job submission. Default="20120521SSINGS"')
    parser.add_argument('-n', action='store', dest='name', default='default_launcher_job', help='The name of your job. Default="default_launcher_job"')
    parser.add_argument('-j', action='store', dest='job', default='./commands', help='The name of the job file containing your commands. Default="./commands"')
    parser.add_argument('-l', action='store', dest='launcher', default='launcher.sge', help='The name of the *.sge launcher script that will be created. Default="launcher.sge"')
    
    results = parser.parse_args()
    
    if results.name == None:
        print 'You did not give a job name.'
        parser.print_help()
        return
    if results.time == None:
        print 'You did not give a job time.'
        parser.print_help()
        return
    if results.job == None:
        print 'You did not give a job file.'
        parser.print_help()
        return
    if results.launcher == None:
        print 'You did not give a launcher save name.'
        parser.print_help()
        return
    
    
    num_cores = file_len(results.job)
    print 'Job file has %i lines.' % num_cores
    while num_cores % 12 != 0:
        num_cores += 1
    
    print 'Using %i cores.' % num_cores
    
    launcher_file = open(results.launcher, 'w')
    
    if results.email == None:
        launcher_file.write('''\
#!/bin/csh
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
#$ -N %s
#$ -pe 12way %i
#$ -q %s
#$ -o %s.o$JOB_ID
#$ -l h_rt=%s
#$ -V
#$ -cwd
#   <------ You MUST Specify a Project String ----->
#$ -A %s
#------------------------------------------------------
#
# Usage:
#	#$ -pe <parallel environment> <number of slots> 
#	#$ -l h_rt=hours:minutes:seconds to specify run time limit
# 	#$ -N <job name>
# 	#$ -q <queue name>
# 	#$ -o <job output file>
#	   NOTE: The env variable $JOB_ID contains the job id. 
#
module load launcher
setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
setenv CONTROL_FILE   %s
setenv WORKDIR        .
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

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
''' % (results.name,
               num_cores,
               results.queue,               
               results.name,
               results.time,
               results.allocation,             
               results.job))
    else:
        launcher_file.write('''\
#!/bin/csh
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
#$ -N %s
#$ -pe 12way %i
#$ -q %s
#$ -o %s.o$JOB_ID
#$ -l h_rt=%s
#$ -V
#$ -M %s
#$ -m be
#$ -cwd
#   <------ You MUST Specify a Project String ----->
#$ -A %s
#------------------------------------------------------
#
# Usage:
#	#$ -pe <parallel environment> <number of slots> 
#	#$ -l h_rt=hours:minutes:seconds to specify run time limit
# 	#$ -N <job name>
# 	#$ -q <queue name>
# 	#$ -o <job output file>
#	   NOTE: The env variable $JOB_ID contains the job id. 
#
module load launcher
setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
setenv CONTROL_FILE   %s
setenv WORKDIR        .
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

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif


#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
''' % (results.name,
               num_cores,
               results.queue,
               results.name,               
               results.time,
               results.email,
               results.allocation,
               results.job))
    print 'Launcher successfully created. Type "qsub %s" to queue your job.' % results.launcher

if __name__ == '__main__':
    main()
