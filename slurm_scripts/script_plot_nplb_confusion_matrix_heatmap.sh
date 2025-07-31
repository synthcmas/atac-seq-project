#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=plot_nplb_confusion_matrix_heatmap
#SBATCH --error=job.%plot_nplb_confusion_matrix_heatmap.err
#SBATCH --output=job.%plot_nplb_confusion_matrix_heatmap.out
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
python3 ../python/nplb/plot_nplb_confusion_matrix_heatmap.py <model_arch_txt> <model_output_txt> <labels_csv> <motifs_csv> <cell_types_csv>

