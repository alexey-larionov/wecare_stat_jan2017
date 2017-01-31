#!/bin/bash

# s02_filter_wecare_nfe_kgen50.sh
# Filter wecare-NFE-kgen50: keep only variants present in both datasets (wecare-NFE and kgen50)
# Alexey Larionov, 31Jan2016

# Use:
# sbatch s02_filter_wecare_nfe_kgen50.sh

# ------------------------------------ #
#         sbatch instructions          #
# ------------------------------------ #

#SBATCH -J filter_wecare_nfe_kgen50
#SBATCH --time=01:00:00
#SBATCH -A TISCHKOWITZ-SL2
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH -p sandybridge
#SBATCH -o s02_filter_wecare_nfe_kgen50.log
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
echo "Filter wecare-NFE-kgen50: keep only variants present in both datasets (wecare-NFE and kgen50)"
echo "Started s02_filter_wecare_nfe_kgen50.sh: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Tools and resources
java="/scratch/medgen/tools/java/jre1.8.0_40/bin/java"
gatk="/scratch/medgen/tools/gatk/gatk-3.7-0/GenomeAnalysisTK.jar"
ref_genome="/scratch/medgen/resources/gatk_bundle/b37/decompressed/human_g1k_v37.fasta"

# Files and folders
working_folder="/scratch/medgen/scripts/wecare_stat_01.17/scripts"
data_folder="/scratch/medgen/scripts/wecare_stat_01.17/interim_data"

wecare_nfe_kgen50_raw_vcf="${data_folder}/wecare_nfe_kgen50.vcf"
wecare_nfe_kgen50_filt_vcf="${data_folder}/wecare_nfe_kgen50_filt.vcf"
filtering_log="${working_folder}/s02_filtering.log"

cd "${working_folder}"

# Report to log
echo "Tools and resources:"
echo "java: ${java}"
echo "gatk: ${gatk}"
echo "ref_genome: ${ref_genome}"
echo ""
echo "Files and folders:"
echo "working_folder: ${working_folder}"
echo "data_folder: ${data_folder}"
echo ""
echo "wecare_nfe_kgen50_raw_vcf: ${wecare_nfe_kgen50_raw_vcf}"
echo "wecare_nfe_kgen50_filt_vcf: ${wecare_nfe_kgen50_filt_vcf}"
echo "filtering_log: ${filtering_log}"
echo ""
echo "Started filtering vcf"

# Combine vcf-s
"${java}" -Xmx60g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -V "${wecare_nfe_kgen50_raw_vcf}" \
  -o "${wecare_nfe_kgen50_filt_vcf}" \
  -select "!empty(kgenVarID) and !empty(SplitVarID)" \
  -nt 12 &>> "${filtering_log}"

# Progress report
echo "Completed filtering vcf: $(date +%d%b%Y_%H:%M:%S)"
echo ""
echo "Num of source variants:"
grep -v "^#" "${wecare_nfe_kgen50_raw_vcf}" | wc -l
echo ""
echo "Num of filtered variants:"
grep -v "^#" "${wecare_nfe_kgen50_filt_vcf}" | wc -l
echo ""

# Completion message
echo "Done all tasks"
echo ""
