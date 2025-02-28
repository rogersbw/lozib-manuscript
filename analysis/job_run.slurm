#!/bin/bash
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --mem=30G  # Requested Memory
#SBATCH -t 48:00:00  # Job time limit
#SBATCH -o slurm_output/slurm-%j.out  # %j = job ID
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=bwrogers@umass.edu

module load r-rocker-ml-verse/4.2.3+apptainer
shopt -s expand_aliases

# Ensure both arguments are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Both first and second arguments are required."
    exit 1
fi

# Create the output directory based on the second argument
OUTPUT_DIR="slurm_output/$2"
mkdir -p "$OUTPUT_DIR"

# Construct the output file name with both user inputs and job ID
OUTPUT_FILE="$OUTPUT_DIR/slurm-${1}-${2}-${SLURM_JOB_ID}.out"

# Redirect stdout and stderr to the new output file
exec > >(tee -a "$OUTPUT_FILE") 2>&1


Rscript --no-restore --quiet --no-save run_fit.r $1 $2
