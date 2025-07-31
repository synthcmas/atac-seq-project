# README

## Project Structure

The project consists of the following directories:

- `R/` - Contains the main R script `run_model.R`.
- `python/` - Contains two subdirectories:
  - `mlr/` - Scripts for running multinomial logistic regression.
  - `nplb/` - Script for plotting confusion matrix and heatmap for NPLB results.
- `jaspar_motifs/` - Contains JASPAR motif files.
- `meme_motifs/` - Contains MEME motif files.
- `data/{dataset_name}/` - Contains dataset-specific files required for running the scripts.
- `slurm_scripts/` - Contains SLURM job scripts for running the analysis on a cluster.

## Prerequisites

Ensure the following dependencies are installed before running the scripts:

### R Packages:

The required R packages are listed in `R/packages.txt`. Install them using the appropriate R package manager.

### Python Dependencies:

Install the required Python libraries using:

```bash
pip install -r python/requirements.txt
```

### Other Dependencies:

- MEME suite
- bedtools

## Running the Main R Script

### Script: `R/run_model.R`

This script runs the main analysis pipeline. The parameters need to be manually set within the script before execution. Key parameters include:

- `cisTopic` (boolean): If `TRUE`, performs cisTopic; otherwise, runs the custom method.
- `same` (boolean): Whether to include same motif pairs.
- `n_seeds_cisTopic` (integer): Number of cisTopic runs (if `cisTopic` or `topic_modeling` is `TRUE`).
- `binary` (boolean): Whether to use raw ATAC-seq count matrix (`FALSE`) or binarized version (`TRUE`).
- `chromVAR_preproc` (boolean): If `TRUE`, applies chromVAR preprocessing.
- `FIMO` (boolean): If `TRUE`, uses FIMO for motif finding.
- `n_cores` (integer): Number of cores for multi-threading.
- `withPairs` (boolean): If `TRUE`, includes motif pair occurrences.
- `chromVAR_norm` (boolean): If `TRUE`, applies chromVAR normalization.
- `bias_correct` (boolean): If `TRUE`, performs bias correction.
- `chromVAR_postproc` (boolean): If `TRUE`, applies chromVAR post-processing.
- `chromVAR_downstream` (boolean): If `TRUE`, runs chromVAR’s downstream clustering pipeline instead of NPLB.
- `topic_modeling` (boolean): If `TRUE`, applies cisTopic to binarized motif activity.
- `data_dir` (string): Path to dataset directory (e.g., `data/furlong`).

### Running on a Cluster

Instead of running the R script manually, use the SLURM job script:

```bash
cd slurm_scripts/
sbatch script_run_model.sh
```

## Required Data Files in `data/{dataset_name}/`

1. `dataset.RData` - A sparse matrix with ATAC-seq count values (columns as cell barcodes, rows as ATAC-seq peaks).
2. `labels.RData` - Vector of true labels for the samples (integer values).
3. `bed_file.bed` - BED file containing ATAC-seq peaks.
4. `regions.fa` - ATAC-seq peaks in FASTA format (for FIMO analysis).
5. `genome_org.txt` - Two-line file specifying genome name in BSgenome format (e.g., `BSgenome.Hsapiens.UCSC.hg38`) and organism name.
6. `lib_size.RData` - Library sizes of samples (required if `chromVAR_preproc` is `TRUE`).

Additionally, the following motif files must be present:

- `jaspar_motifs/{organism_name}.txt` - List of motifs in JASPAR format.
- `meme_motifs/{organism_name}.meme` - List of motifs in MEME format.

## Running Python Scripts

Python scripts are located in `python/`. Instead of running them manually, use the SLURM job scripts. **Ensure the required positional arguments are correctly provided inside the job script for the Python commands.**

### `smlr/`

#### Run Multiple Logistic Regression across multiple runs

```bash
cd slurm_scripts/
sbatch script_run_mlr_precision_full_kfold_multiple_runs.sh
```

Example SLURM script (`script_run_mlr_precision_full_kfold_multiple_runs.sh`):

```bash
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

python3 ../python/mlr/run_mlr_precision_full_kfold_multiple_runs.py <motif_activity_file_csv> <labels_csv>
```

#### Plot Error Bars

```bash
sbatch script_plot_error_bars.sh
```

Example SLURM script (`script_plot_error_bars.sh`):

```bash
#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=plot_error_bars
#SBATCH --error=job.%plot_error_bars.err
#SBATCH --output=job.%plot_error_bars.out
#SBATCH --partition=gpu

# Load required modules
module load anaconda3-2023.3

python3 ../python/mlr/plot_error_bars.py <pairs_precision_results_pkl> <individual_precision_results_pkl> <cell_types_csv>
```

#### Run MLR and Save Feature Weights

```bash
sbatch script_run_mlr_feature_weights_full.sh
```

Example SLURM script (`script_run_mlr_feature_weights_full.sh`):

```bash
#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=feature_weights
#SBATCH --error=job.%feature_weights.err
#SBATCH --output=job.%feature_weights.out
#SBATCH --partition=gpu

# Load required modules
module load anaconda3-2023.3

python3 ../python/mlr/run_mlr_feature_weights_full.py <motif_activity_file_csv> <labels_csv>
```

#### Plot Important Motifs

```bash
sbatch script_plot_imp_motifs.sh
```

Example SLURM script (`script_plot_imp_motifs.sh`):

```bash
#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=plot_imp_motifs
#SBATCH --error=job.%plot_imp_motifs.err
#SBATCH --output=job.%plot_imp_motifs.out
#SBATCH --partition=gpu

# Load required modules
module load anaconda3-2023.3

python3 ../python/mlr/plot_imp_motifs.py <pairs_feature_weights_pkl> <individual_feature_weights_pkl>
```

### `nplb/`

#### Plot Confusion Matrix and Heatmap

```bash
sbatch script_plot_nplb_confusion_matrix_heatmap.sh
```

Example SLURM script (`script_plot_nplb_confusion_matrix_heatmap.sh`):

```bash
#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=plot_nplb_confusion_matrix_heatmap.py
#SBATCH --error=job.%plot_nplb_confusion_matrix_heatmap.err
#SBATCH --output=job.%plot_nplb_confusion_matrix_heatmap.out
#SBATCH --partition=gpu

# Load required modules
module load anaconda3-2023.3

python3 ../python/nplb/plot_nplb_confusion_matrix_heatmap.py <model_arch_txt> <model_output_txt> <labels_csv> <motifs_csv> <cell_types_csv>
```

## Troubleshooting

- **R script fails due to missing packages:** Ensure all R dependencies are installed using Bioconductor and CRAN.
- **Python script errors (e.g., ImportError):** Check that all dependencies are installed via `pip install -r python/requirements.txt`.
- **Incorrect file paths:** Ensure dataset-specific files exist in `data/{dataset_name}/` and paths in scripts are correctly set.
- **Memory issues:** If encountering memory errors, reduce the dataset size or adjust memory allocation parameters.
- **Multithreading issues:** Ensure `n_cores` is set appropriately for available CPU resources.
- **Unexpected results:** Double-check input files, parameter settings, and preprocessing steps.

For any persistent issues, refer to error logs and confirm all input files and dependencies are correctly set up.

