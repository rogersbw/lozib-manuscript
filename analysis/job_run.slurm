#!/bin/bash
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=20G  # Requested Memory
#SBATCH -p cpu-long  # Partition
#SBATCH -t 24:00:00  # Job time limit
#SBATCH -o slurm-%j.out  # %j = job ID
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bwrogers@umass.edu

Rscript run_fit.R $1
