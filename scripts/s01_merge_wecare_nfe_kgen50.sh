#!/bin/bash

# s01_merge_wecare_nfe_kgen50.sh
# Merge wecare-nfe and kgen50 vcfs
# Alexey Larionov, 27Jan2016

# Use:
# sbatch s01_merge_wecare_nfe_kgen50.sh

# ------------------------------------ #
#         sbatch instructions          #
# ------------------------------------ #

#SBATCH -J merge_wecare_nfe_kgen50
#SBATCH --time=03:00:00
#SBATCH -A TISCHKOWITZ-SL2
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
#SBATCH -o s01_merge_wecare_nfe_kgen50.log

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
echo "Merge wecare-nfe and kgen50 vcfs"
echo "Started s01_merge_wecare_nfe_kgen50.sh: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Initial settings --- #

# Tools and resources
java8="/scratch/medgen/tools/java/jre1.8.0_40/bin/java"
gatk360="/scratch/medgen/tools/gatk/gatk-3.6-0/GenomeAnalysisTK.jar"
ref_genome="/scratch/medgen/resources/gatk_bundle/b37/decompressed/human_g1k_v37.fasta"
bcftools="/scratch/medgen/tools/bcftools/bcftools-1.2/bin/bcftools"
plot_vcfstats="/scratch/medgen/tools/bcftools/bcftools-1.2/bin/plot-vcfstats"

# Files and folders
working_folder="/scratch/medgen/scripts/wecare_stat_01.17/scripts"
target_folder="/scratch/medgen/scripts/wecare_stat_01.17/interim_data"

wecare_nfe_vcf="/scratch/medgen/users/alexey/wecare_6_jan2017/wecare_nfe_nov2016_vqsr_shf_sma_ann/wecare_nfe_nov2016_vqsr_shf_sma_ann.vcf"
kgen50_vcf="/scratch/medgen/resources/phase3_1k_release20130502/selected_samples/kgen_selected_50_cat_sma_cln_id.vcf"

wecare_nfe_kgen50_vcf="${target_folder}/wecare_nfe_kgen50.vcf"
combining_log="${working_folder}/s01_combining.log"

vcfstats_folder="${target_folder}/vcfstats"
mkdir -p "${vcfstats_folder}"

vcf_stats="${vcfstats_folder}/wecare_nfe_kgen50.vchk"

# Report to log
echo "Tools and resources:"
echo "java8: ${java8}"
echo "gatk360: ${gatk360}"
echo "ref_genome: ${ref_genome}"
echo "bcftools: ${bcftools}"
echo "plot_vcfstats: ${plot_vcfstats}"
echo ""
echo "Files and folders:"
echo "working_folder: ${working_folder}"
echo "target_folder: ${target_folder}"
echo "wecare_nfe_vcf: ${wecare_nfe_vcf}"
echo "kgen50_vcf: ${kgen50_vcf}"
echo "wecare_nfe_kgen50_vcf: ${wecare_nfe_kgen50_vcf}"
echo "combining_log: ${combining_log}"
echo "vcfstats_folder: ${vcfstats_folder}"
echo "vcf_stats: ${vcf_stats}"
echo ""

initial_folder="$(pwd -P)"
cd "${working_folder}"

# --- Combine vcf-s --- #

echo "Num of variants in wecare-nfe vcf:"
grep -v "^#" "${wecare_nfe_vcf}" | wc -l
echo ""

echo "Num of variants in kgen-50 vcf:"
grep -v "^#" "${kgen50_vcf}" | wc -l
echo ""

# Progress report
echo "Started combining vcf-s"

# Combine vcf-s
"${java8}" -Xmx60g -jar "${gatk360}" \
  -T CombineVariants \
  -R "${ref_genome}" \
  -V "${wecare_nfe_vcf}" \
  -V "${kgen50_vcf}" \
  -o "${wecare_nfe_kgen50_vcf}" \
  -genotypeMergeOptions UNIQUIFY \
  -env \
  -nt 12 &>> "${combining_log}"

# UNIQUIFY: Make all sample genotypes unique by file. 
#           Each sample shared across VCFs gets named sample.VCF
# env: Exclude sites where no variation is present after merging

# Progress report
echo "Num of combined variants:"
grep -v "^#" "${wecare_nfe_kgen50_vcf}" | wc -l
echo ""

# --- Run vcf-stats --- #

echo "Started calculating vcfstats"
echo ""

# Calculate vcf stats
"${bcftools}" stats -F "${ref_genome}" "${wecare_nfe_kgen50_vcf}" > "${vcf_stats}" 

# Plot the stats
"${plot_vcfstats}" "${vcf_stats}" -p "${vcfstats_folder}/"
echo ""

# Progress report
echo "Completed calculating vcfstats: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Final tasks --- #

# Restore working folder
cd "${initial_folder}"

# Completion message
echo "Done all tasks"
echo ""
