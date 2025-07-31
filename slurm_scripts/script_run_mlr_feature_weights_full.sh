#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=feature_weights
#SBATCH --error=job.%feature_weights.err
#SBATCH --output=job.%feature_weights.out
#SBATCH --partition=gpu

#source /etc/profile.d/modules.sh
# Load necessary modules
#module load meme
#/apps/codes/R-4.3.1/bin/Rscript /home/gaurik/nnp/nnp/test_FIMO.R
#export PATH="/apps/codes/meme-5.5.5/bin:$PATH"
#export PATH="/apps/codes/R-4.3.1/bin:$PATH"
#!/bin/bash

# Load environment modules
#. /etc/profile.d/modules.sh

#export PATH="/apps/codes/R-4.3.1/bin:$PATH"
#export PATH="/apps/codes/meme-5.5.5/bin:$PATH"
# Load required modules
module load anaconda3-2023.3

# Run the R script
python3 ../python/mlr/run_mlr_feature_weights_full.py <motif_activity_file_csv> <labels_csv>

