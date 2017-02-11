#!/bin/bash

# Runs smartpca and makes PCA plots using EIGENSRAT package
# Started: Alexey Larionov, Nov 2016
# Last updated: Alexey Larionov, 10Feb2017

# Use:
#sbatch u01_run_smartpca.sh parameters_file populations
#sbatch u01_run_smartpca.sh wecare_only_480_226k_eigenstrat_default_outliers.par UBC:CBC

# Notes:
# Before running the script
# - Make sure that parameters file is consistent with the script
# - It is not necessary to make the output folders: they will be created by the script
# After running the script:
# - Manually REMOVE COLONS ":" from the plots files names !!!
#   Otherwise thiese files may not be opened/copied etc.
# - Manually rename the log (u01_run_smartpca.log).
#   Otherwise it could be overwritten when the sript runs next time. 

# ------------------------------------ #
#         sbatch instructions          #
# ------------------------------------ #

#SBATCH -J smartpca
#SBATCH --time=01:00:00
#SBATCH -A TISCHKOWITZ-SL2
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
#SBATCH -o u01_run_smartpca.log
#SBATCH --qos=INTR

# Minimal required modules (not to be removed)
. /etc/profile.d/modules.sh
module purge
module load default-impi

# Set initial working folder
cd "${SLURM_SUBMIT_DIR}"

## Report settings and run the job
echo "Job name: ${SLURM_JOB_NAME}"
echo "Allocated node: $(hostname)"
echo "Initial working folder:"
echo "${SLURM_SUBMIT_DIR}"
echo ""
echo " ------------------ Output ------------------ "
echo ""

# ------------------------------------ #
#        Start of the script           #
# ------------------------------------ #

# --- Set environment, files and folders --- #

# Stop at errors
set -e 

# Read parameter
parameters_file="${1}"
populations_to_plot="${2}"

# Analysis name (parameters file name w/o extension)
analysis=$(basename "${parameters_file}")
analysis="${parameters_file%.par}"

# Add smartpca bin directory to path
export PATH="/scratch/medgen/tools/eigensoft/EIG-6.1.4/bin/:$PATH"

# Files and folders
working_folder="/scratch/medgen/scripts/wecare_stat_01.17/scripts"
cd "${working_folder}"

data_folder="/scratch/medgen/scripts/wecare_stat_01.17/interim_data"

output_folder="${data_folder}/${analysis}"
mkdir -p "${output_folder}"

evec_file="${output_folder}/${analysis}.evec"
smartpca_log="${output_folder}/${analysis}.log"
gnuplot_file="${output_folder}/${analysis}.xtxt" # must end .xtxt

# --- Report settings to log --- #

echo "Run smartpca and make PCA plots using EIGENSRAT package"
echo ""
echo "Started u01_run_smartpca.sh: $(date +%d%b%Y_%H:%M:%S)"
echo ""
echo "Files and folders"
echo "parameters_file: ${parameters_file}"
echo "populations_to_plot: ${populations_to_plot}"
echo ""
echo "analysis: ${analysis}"
echo "working_folder: ${working_folder}"
echo "data_folder: ${data_folder}"
echo "output_folder: ${output_folder}"
echo " "
echo "evec_file: ${evec_file}"
echo "smartpca_log: ${smartpca_log}"
echo "gnuplot_file: ${gnuplot_file}"
echo ""
echo "Parameters file:"
echo "$(cat ${parameters_file})"
echo ""

# --- Do the analysis --- #

# Run smartpca
smartpca -p "${parameters_file}" > "${smartpca_log}"

# Plot pc1 vs pc2
ploteig \
  -i "${evec_file}" \
  -c 1:2 \
  -p "${populations_to_plot}" \
  -o "${gnuplot_file}" \
  -x

# --- Completion message --- #

echo "Completed: $(date +%d%b%Y_%H:%M:%S)"
echo ""
