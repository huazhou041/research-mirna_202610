#!/bin/bash
#SBATCH --partition=cpu_long,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fork3corr
#SBATCH --output=logs/fork3corr_%A_%a.log
#SBATCH --array=1-96

# CV-honest ComBat + SUPERVISED correlation feature selection (train-fold only).
# One array task = one (variant, learner, varcut) cell. 48h walltime because raw
# GLM at high varcut is slow (see fork3). Params from codes/fork3_jobs.txt.
set +u
module load miniconda3/4.6.14
source activate r_env
set -u

cd /gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork3_corr_fs
mkdir -p res logs raw_pred

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" codes/fork3_jobs.txt)
variant=$(echo "$line" | awk '{print $1}')
learner=$(echo "$line" | awk '{print $2}')
varcut=$(echo  "$line" | awk '{print $3}')
echo "task ${SLURM_ARRAY_TASK_ID}: variant=$variant learner=$learner varcut=$varcut fs=corr"

Rscript codes/tune_combined_models_cvcombat.R \
  rdata/combined_plain_miRNA.RDS \
  rdata/combined_variants_metadata.csv \
  "res/combined_${variant}_miRNA" \
  --variant "$variant" --learner "$learner" --varcut "$varcut" \
  --fs_method corr --nboot 100 --kfold 10
