#!/bin/bash

#SBATCH --job-name=TimeCalibrateTrees                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#SBATCH --tmp=16G
#Total memory reserved: 16GB

#SBATCH --time=168:00:00        #This is the time that your task will run
#SBATCH --qos=1week           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/timecalib_stdout.txt     #These are the STDOUT and STDERR files
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/timecalib_stderr.txt

#You selected an array of jobs from 1 to n with n simultaneous jobs
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=max.shafer@gmail.com        #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load R/4.0.0-foss-2018b

#export your required environment variables below
#################################################

export PATH=~/.local/bin:~/.local:$PATH
export LD_LIBRARY_PATH=~/.local/lib:~/.local/lib64:${exec_prefix}/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=~/.local/lib:~/.local/lib64:$LIBRARY_PATH
export CPATH=~/.local/include:~/.local/include/adolc:$CPATH

#add your command lines below
#############################

# Rscript 1_Eutelostomi_tree.R
Rscript 2_Eutelostomi_Hidden_Markov_Models_script.R

