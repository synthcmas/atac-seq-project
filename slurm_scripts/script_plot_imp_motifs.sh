#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=plot_imp_motifs
#SBATCH --error=job.%plot_imp_motifs.err
#SBATCH --output=job.%plot_imp_motifs.out
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
python3 ../python/mlr/plot_imp_motifs.py <pairs_feature_weights_pkl> <individual_feature_weights_pkl>

