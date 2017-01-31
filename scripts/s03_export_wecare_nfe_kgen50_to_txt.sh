#!/bin/bash

# s03_export_wecare_nfe_kgen50_to_txt.sh
# Export wecare_nfe_kgen50 vcf file to plain text tables for downstream analysis in R
# Alexey Larionov, 31Jan2016

# Use:
# sbatch s02_export_wecare_nfe_kgen50_to_txt.sh

# ------------------------------------ #
#         sbatch instructions          #
# ------------------------------------ #

#SBATCH -J export_wecare_nfe_kgen50_to_txt
#SBATCH --time=01:00:00
#SBATCH -A TISCHKOWITZ-SL2
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
#SBATCH -o s03_export_wecare_nfe_kgen50_to_txt.log
#SBATCH --qos=INTR

# Minimal required modules (not to be removed)
. /etc/profile.d/modules.sh
module purge
module load default-impi

# Additional modules for vcfstats
module load python/2.7.5
module load pandoc/1.15.2.1

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

# Stop at errors
set -e

# Start message
echo "Export wecare_nfe_kgen50 vcf file to plain text tables for downstream analysis in R"
echo "Started s02_export_wecare_nfe_kgen50_to_txt.sh: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Files and settings --- #

# Tools and resources
java="/scratch/medgen/tools/java/jre1.8.0_40/bin/java"
gatk="/scratch/medgen/tools/gatk/gatk-3.7-0/GenomeAnalysisTK.jar"
ref_genome="/scratch/medgen/resources/gatk_bundle/b37/decompressed/human_g1k_v37.fasta"

# Files and folders
working_folder="/scratch/medgen/scripts/wecare_stat_01.17/scripts"
data_folder="/scratch/medgen/scripts/wecare_stat_01.17/interim_data"
wecare_nfe_kgen50_vcf="${data_folder}/wecare_nfe_kgen50_filt.vcf"

vcf_txt="${data_folder}/wecare_nfe_kgen50_filt_vcf.txt"
vcf_log="${working_folder}/s03_wecare_nfe_kgen50_filt_vcf_to_txt.log"

gt_txt="${data_folder}/wecare_nfe_kgen50_filt_gt.txt"
gt_log="${working_folder}/s03_wecare_nfe_kgen50_filt_gt_to_txt.log"

gq_txt="${data_folder}/wecare_nfe_kgen50_filt_gq.txt"
gq_log="${working_folder}/s03_wecare_nfe_kgen50_filt_gq_to_txt.log"

dp_txt="${data_folder}/wecare_nfe_kgen50_filt_dp.txt"
dp_log="${working_folder}/s03_wecare_nfe_kgen50_filt_dp_to_txt.log"

# Report to log
echo "Tools and resources:"
echo "java: ${java}"
echo "gatk360: ${gatk}"
echo "ref_genome: ${ref_genome}"
echo ""
echo "Files and folders:"
echo "working_folder: ${working_folder}"
echo "data_folder: ${data_folder}"
echo "wecare_nfe_kgen50_vcf: ${wecare_nfe_kgen50_vcf}"
echo ""
echo "vcf_txt: ${vcf_txt}"
echo "vcf_log: ${vcf_log}"
echo ""
echo "gt_txt: ${gt_txt}"
echo "gt_log: ${gt_log}"
echo ""
echo "gq_txt: ${gq_txt}"
echo "gq_log: ${gq_log}"
echo ""
echo "dp_txt: ${dp_txt}"
echo "dp_log: ${dp_log}"
echo ""

# Set working folder
initial_folder="$(pwd -P)"
cd "${working_folder}"

# --- Export VCF table --- #

# Progress report
echo "Started exporting VCF table"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -V "${wecare_nfe_kgen50_vcf}" \
  -o "${vcf_txt}" \
  -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F kgenVarID -F SplitVarID -F AF -F AC -F AN \
  -AMD &>  "${vcf_log}"
# -AMD : allow missed data

# Progress report
echo "Completed exporting VCF table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export GT table --- #

# Progress report
echo "Started exporting GT table"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -V "${wecare_nfe_kgen50_vcf}" \
  -o "${gt_txt}" \
  -F CHROM -F POS -F kgenVarID -F SplitVarID -GF GT \
  -AMD &>  "${gt_log}"  

# Progress report
echo "Completed exporting GT table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export GQ table --- #

# Progress report
echo "Started exporting GQ table"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -V "${wecare_nfe_kgen50_vcf}" \
  -o "${gq_txt}" \
  -F CHROM -F POS -F kgenVarID -F SplitVarID -GF GQ \
  -AMD &>  "${gq_log}"  

# Progress report
echo "Completed exporting GQ table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Export DP table --- #

# Progress report
echo "Started exporting DP table"

# Export table
"${java}" -Xmx60g -jar "${gatk}" \
  -T VariantsToTable \
  -R "${ref_genome}" \
  -V "${wecare_nfe_kgen50_vcf}" \
  -o "${dp_txt}" \
  -F CHROM -F POS -F kgenVarID -F SplitVarID -GF DP \
  -AMD &>  "${dp_log}"  

# Progress report
echo "Completed exporting DP table: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Final tasks --- #

# Restore working folder
cd "${initial_folder}"

# Completion message
echo "Done all tasks"
echo ""
