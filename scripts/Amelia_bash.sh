#!/bin/bash

#above is the shebang statement which tells the shell to execute the script via bash shell

#SBATCH --job-name=Artiodactyla_diel_models           #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#SBATCH --tmp=16G
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/Amelia_artiodactyla_models_stdout.txt
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/Amelia_artiodactyla_models_stderr.txt

#You selected an array of jobs from 1 to n with n simultaneous jobs
#SBATCH --mail-type=END,FAIL,TIME_LIMIT,COMPLETE
#SBATCH --mail-user=amelia.mesich@mail.utoronto.ca       #You will be notified via email when your task ends or fails

#This job runs from the current working directory


#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $JOB_ID stores the ID number of your task.


#load your required modules below
#################################

module load R/4.3.1

#export your required environment variables below
#################################################

export PATH=~/.local/bin:~/.local:$PATH
export LD_LIBRARY_PATH=~/.local/lib:~/.local/lib64:${exec_prefix}/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=~/.local/lib:~/.local/lib64:$LIBRARY_PATH
export CPATH=~/.local/include:~/.local/include/adolc:$CPATH

#add your command lines below
#############################

#Rscript Amelia_slurm_script.R max_dinoc artiodactyla ARD bridge_only
Rscript Amelia_max_clade_slurm.R six_state artiodactyla ARD

# For the corHMM script, requires three arguments
# arg1 (states in the model) can be "max_crep" or "max_dinoc" or "six_state" or "four"
# arg2 (phylogenetic tree) can be "cetaceans", "artiodactyla", "artiodactyla_minus_cetaceans"
# arg3 (model) can be "ER", "SYM", "ARD" or "bridge_only"
