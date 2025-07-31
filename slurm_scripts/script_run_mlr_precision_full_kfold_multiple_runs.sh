#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=precision
#SBATCH --error=job.%precision.err
#SBATCH --output=job.%precision.out
#SBATCH --partition=gpu

# Load required modules
module load anaconda3-2023.3

# Run the R script
python3 ../python/mlr/run_mlr_precision_full_kfold_multiple_runs.py <motif_activity_file_csv> <labels_csv>

