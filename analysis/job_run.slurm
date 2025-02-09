#!/bin/bash
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=20G  # Requested Memory
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o slurm_output/slurm-%j-$1-$2.out  # %j = job ID
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bwrogers@umass.edu

module load r-rocker-ml-verse/4.2.3+apptainer
shopt -s expand_aliases

Rscript --no-restore --quiet --no-save run_fit.r $1 $2
