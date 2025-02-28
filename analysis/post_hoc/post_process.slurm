#!/bin/bash
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --mem=30G  # Requested Memory
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o slurm_output/slurm-%j.out  # %j = job ID
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bwrogers@umass.edu

module load r-rocker-ml-verse/4.2.3+apptainer

Rscript --no-restore --quiet --no-save post_summary_means.R
