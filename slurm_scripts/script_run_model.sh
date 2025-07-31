#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=run_model
#SBATCH --error=job.%run_model.err
#SBATCH --output=job.%run_model.out
#SBATCH --partition=gpu

#source /etc/profile.d/modules.sh
# Load necessary modules
#module load meme
#/apps/codes/R-4.3.1/bin/Rscript /home/gaurik/nnp/nnp/test_FIMO.R
#export PATH="/apps/codes/meme-5.5.5/bin:$PATH"
#export PATH="/apps/codes/R-4.3.1/bin:$PATH"
#!/bin/bash

# Load environment modules
. /etc/profile.d/modules.sh

export PATH="/apps/codes/R-4.3.1/bin:$PATH"
export PATH="/apps/codes/meme-5.5.5/bin:$PATH"
# Load required modules
module load bedtools2
module load meme-5.5.5
module load R-4.3.1

# Run the R script
Rscript ../R/run_model.R

