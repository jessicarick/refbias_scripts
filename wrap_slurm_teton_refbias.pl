#!/usr/bin/env perl

### Version 1.0 -- 14 April 2015 -- Alex Buerkle

### wrap_slurm_moran.pl : a perl script to submit many serial slurm
### jobs to the queue using the sbatch command.

use warnings;
use strict;

### CHANGES -----------------------------------------------

### v1.0 - Alex - update of wrap_qsub to work with the SLURM system
### v1.1 - J. Rick -- modified to work on Teton system
#
## possible future addition: specify memory requirement per node, or other hardware requirements


### ------------------ JOB CONFIGURATION --------------------------------------------------------------- 
### all variables that typically need to be modified are in the
### following block
my $arccproject='phylogenref';
my $runtime = '1-00:00:00'; ## 14-00 is 14 days (days-hours notation), 1:00:00 is 1 hour
#my $runtime = '12:00:00';
# specify a reasonable value here.  If the job does not finish by the
# time reached, the job is terminated.  Your job should get greater
# priority if you specify a minimal walltime
my $sendmail = 'TRUE'; ##  if TRUE, sends an email when each job completes

## build array of jobs to be run individually (serially) by pbs
my @jobarray = ();

##### simple job array using gsl-randist, good for debugging
#my $outfile;
#foreach my $j (0..2){
#    $outfile = "test-gsl-out-"."$j".".txt";
#    push @jobarray, "gsl-randist 19187 10000 binomial 0.5 25 > $outfile\n";
#}


### RUNNING THE SIMULATONS ###
my $sim = $ARGV[0];
my $tree_height = $ARGV[1];
#if ($#ARGV != 2 ) {
#	print "ERROR: need ONE simulation number and one tree height to proceed\n";
#	exit;
#}

push @jobarray, "/project/phylogenref/scripts/sim_scripts/test_simphy_ExtREFs_090421.sh $sim $tree_height \n";
push @jobarray, "/project/phylogenref/scripts/sim_scripts/test_simphy_IntREFs_090421.sh $sim $tree_height \n";


### -------------------------END JOB CONFIGURATION---------------------------------------------------------
### UNLIKELY you need to change anything beyond here (except maybe module loading)
### automatically configured or infrequently changed variables
my $emailaddress = $ENV{USER}.'\@uwyo.edu'; ### mail address for notification that a job has completed 
my $jobname = 'slurm.refbias.'.$ENV{USER};

### modules to load:
my @modules =();
push @modules, 'module load gcc';
push @modules, 'module load miniconda3';
push @modules, 'module load perl';
push @modules, 'module load vcftools';
push @modules, 'module load bwa';
push @modules, 'module load samtools/1.6';
push @modules, 'module load htslib';
push @modules, 'module load bcftools';
push @modules, 'module load raxml';
push @modules, 'module load r';
push @modules, 'module load parallel';

my $workdir = '/lscratch'; 
## Set to /lscratch on mtmoran.  This is local disk with high
## performance.  This is where results will be written temporarily on
## a node. 
my $basestoredir = "/project/$arccproject/scripts";
## use a directory in /project/$arccproject/ as this has more
## allocated disk space than your home directory.  Will write to
## slurm_log and slurm_results inside this directory
my $logdir = "$basestoredir/slurm_log";
my $resultdir = "/gscratch/jrick/phylogenref/slurm_results";
## -----------------------------------------------------------------------------------

printf "Ready to submit %d jobs to SLURM (y/n): ", scalar @jobarray;
my $response = 'y';
chomp $response;
if($response eq 'n'){
    print "Exiting without any SLURM submissions.\n";
    exit;
}
else{
    print "Proceeding with SLURM submission\n";
}
unless(-e $logdir){
    mkdir  $logdir or die "Failed to make log file directory";
}
unless(-e $resultdir ){
    mkdir  $resultdir or die "Failed to make result directory";
}

###### ------------------ BUILD SLURM SCRIPT---------------
my @slurmdirectives = "#!/bin/bash";

push @slurmdirectives, "#SBATCH --account=$arccproject";
push @slurmdirectives, "#SBATCH --job-name=$jobname";
push @slurmdirectives, "#SBATCH --time=$runtime"; 
push @slurmdirectives, "#SBATCH --nodes=1";
push @slurmdirectives, "#SBATCH --ntasks-per-node=16"; # 32 cores per node
push @slurmdirectives, "#SBATCH --mem=64G"; 
push @slurmdirectives, "#SBATCH --no-requeue";
#push @slurmdirectives, "#SBATCH --workdir=$logdir";
#          SLURM can send informative email messages to you about the
#          status of your job.  
if($sendmail eq 'TRUE'){
    push @slurmdirectives, '#SBATCH --mail-type=END';
    push @slurmdirectives, "#SBATCH --mail-user $emailaddress";
}
my $pbsconf = join "\n", @slurmdirectives;

#  By default standard output and error streams are to be merged,
#  intermixed, as standard output.

##########################################
# 
#   Output some useful job information.  #
#
##########################################

my @slurmjob = ();
push @slurmjob, '
    echo ------------------------------------------------------
    echo wrapSLURM: Job is running on node $HOSTNAME
    echo ------------------------------------------------------
    echo wrapSLURM: submitting host was $SLURM_SUBMIT_HOST
    echo wrapSLURM: job identifier is $SLURM_JOB_ID
    echo wrapSLURM: job name is $SLURM_JOB_NAME
    echo -----------------------------------------------------
';

###############################################################                                                    
#   The prolog script automatically makes a directory on the local
#   disks for you.  The name of this directory depends on the job id,
#   but you need only refer to it using ${WORKDIR}.
##############################################################

push @slurmjob, "WORKDIR=$workdir/SLURM_\$SLURM_JOB_ID";
push @slurmjob, "SCP='/usr/bin/scp -o StrictHostKeyChecking=no'";

######################################################################
#   To minimize communications traffic, it is best for your job to
#   work with files on the local disk of the compute node.  Hence, one
#   needs to transfer files from your permanent home directory tree to
#   the directory ${WORKDIR} automatically created by wrapSLURM on the
#   local disk before program execution, and to transfer any important
#   output files from the local disk back to the permanent home
#   directory tree after program execution is completed.  We use
#   secure copy (scp) to do the file transfers to avoid distributed
#   filesystem bottlenecks.
######################################################################

#####################################################
#    Specify the permanent directory(ies) on the server host.  Note
#    that when the job begins execution, the current working directory
#    at the time the qsub command was issued becomes the current
#    working directory of the job.
#####################################################

push @slurmjob, "BASESTOREDIR=$basestoredir";
push @slurmjob, "RESULTDIR=$resultdir";

push @slurmjob, '
    echo workdir is $WORKDIR
    echo basestoredir is $BASESTOREDIR
    echo resultdir is $RESULTDIR
    echo ------------------------------------------------------
    echo \' \'
';

###############################################################
#                                                             #
#    setup WORKDIR (typically local disk)
#                                                             #
###############################################################

my $stagein = '
stagein()
{
    echo \' \'
	echo Setting up local working directory ${WORKDIR}
	echo Optionally transferring files from server to compute node
    mkdir ${WORKDIR}
    cd ${WORKDIR}
    ';
foreach my $module (@modules){
$stagein .= "$module\n";
}

$stagein .= '
    ###    ${SCP} ${BASESTOREDIR}/input_file .
}
';
push @slurmjob, $stagein;

############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#                                                          #
############################################################

push @slurmjob, "runprogram()";
push @slurmjob, "{\n";
## gather everything up to this point, use join with "\n"
my $prolog = join "\n", @slurmjob;
@slurmjob = ();

## this is where we would put the executable if this were not in a perl wrapper:
## program_executable < input_file > output_file

### --- beginning of epilog
push @slurmjob, "}";

###########################################################
#                                                         
#   Copy necessary files back to permanent directory and remove results on node     
#                                                         
###########################################################

push @slurmjob, '
stageout()
{
 echo \' \'
 echo Transferring files from compute nodes to server
 echo Writing files in permanent directory  ${RESULTDIR}
 
 cd ..  ## to parent directory of WORKDIR, for writing
 tar czf SLURM_${SLURM_JOBID}_results.tgz SLURM_${SLURM_JOBID}/
 cp SLURM_${SLURM_JOBID}_results.tgz  ${RESULTDIR}/
tar zxf ${RESULTDIR}/SLURM_${SLURM_JOBID}_results.tgz

 ## clean up
 if [ -e "${RESULTDIR}/SLURM_${SLURM_JOBID}_results.tgz" ]
 then
    rm -rf ${WORKDIR}
    rm  SLURM_${SLURM_JOBID}_results.tgz
 fi
}
';

#####################################################################
#  A slurm command can be used to kill a running job.  It first sends
#  a SIGTERM signal, then after a delay (specified by the "kill_delay"
#  queue attribute (set to 30 seconds), it sends a SIGKILL signal
#  which eradicates the job.  During the time between the SIGTERM and
#  SIGKILL signals, the "cleanup" function below is run. You should
#  include in this function commands to copy files from the local disk
#  back to your home directory.  Note: if you need to transfer very
#  large files which make take longer than 30 seconds, change the
#  KillWait variable
#####################################################################

push @slurmjob, '
early()
{
 echo \' \'
 echo \' ############ WARNING:  EARLY TERMINATION ############# \'
 echo \' \'
}

trap \'early; stageout\' 2 9 15

';

##################################################                                               
#   Staging in, running the job, and staging out were specified above
#   as functions.  Now call these functions to perform the actual file
#   transfers and program execution.
#  #################################################

push @slurmjob, "stagein\n";
push @slurmjob, "runprogram\n";
push @slurmjob, "stageout\n";
my $epilog = join '', @slurmjob;

###### use loop to submit whole @jobarray ########------------
foreach my $job (0..($#jobarray - 1)){
  runserialjob($job);
}
## final job
runserialjob($#jobarray);


#### -------------------------------------------------------------------
sub runserialjob{
    my $j = $_[0];
    my $slurmjob = '';
    $slurmjob .= $pbsconf;
    $slurmjob .= $prolog;
    $slurmjob .= $jobarray[$j];
    $slurmjob .= $epilog;
    $slurmjob .= "exit\n";
    open SBATCH, "| sbatch 1>/dev/null" or die "Failed to fork for sbatch; $!";
    print SBATCH "$slurmjob";
    close(SBATCH) or die "Couldn't close SBATCH";
}

