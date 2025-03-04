#!/bin/bash

#SBATCH --job-name=TimeCalibrateTrees                   #This is the name of your job
#SBATCH --cpus-per-task=1                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=16G              #This is the memory reserved per core.
#SBATCH --tmp=16G
#Total memory reserved: 16GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/standardmodels_stdout.txt
#SBATCH --error=/scicore/home/schiera/gizevo30/projects/fish_sleep/scripts/logs/standardmodels_stderr.txt

#You selected an array of jobs from 1 to n with n simultaneous jobs
#SBATCH --mail-type=END,FAIL,TIME_LIMIT,COMPLETE
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

# Rscript 1_Fish_phylogeny.R

# Rscript 1_Eutelostomi_tree.R

# For the HMM script, requires two arguments
# name_variable can be all, only_highqual, only_cartilaginous, or only_ingroup, or not_mammals
# dataset_variable can be fish or AllGroups, or mammals, or tetrapods

# The following lines run all the current permutations, one after another

# Rscript 2_General_Hidden_Markov_Models_script.R all fish
Rscript 2_General_Hidden_Markov_Models_script.R only_highqual fish
Rscript 2_General_Hidden_Markov_Models_script.R only_cartilaginous fish
Rscript 2_General_Hidden_Markov_Models_script.R only_ingroup fish

# Rscript 2_General_Hidden_Markov_Models_script.R all mammals

# Rscript 2_General_Hidden_Markov_Models_script.R all tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R not_mammals tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R amphibians tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R sauropsids tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R lepidosauria tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R testudines tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R crocodylia tetrapods
# Rscript 2_General_Hidden_Markov_Models_script.R aves tetrapods

# Rscript 2_General_Hidden_Markov_Models_script.R all AllGroups

