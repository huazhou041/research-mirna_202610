#!/bin/bash
#SBATCH --partition=cpu_short,cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=fork3_cvcombat
#SBATCH --output=logs/fork3_%A_%a.log
#SBATCH --array=1-96

# One array task = one (variant, learner, varcut) cell of the CV-honest
# 4-way batch-handling comparison. Params read from codes/fork3_jobs.txt.
set +u
module load miniconda3/4.6.14
source activate r_env
set -u

cd /gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork3
mkdir -p res logs raw_pred

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" codes/fork3_jobs.txt)
variant=$(echo "$line" | awk '{print $1}')
learner=$(echo "$line" | awk '{print $2}')
varcut=$(echo  "$line" | awk '{print $3}')

echo "task ${SLURM_ARRAY_TASK_ID}: variant=$variant learner=$learner varcut=$varcut"

Rscript codes/tune_combined_models_cvcombat.R \
  rdata/combined_plain_miRNA.RDS \
  rdata/combined_variants_metadata.csv \
  "res/combined_${variant}_miRNA" \
  --variant "$variant" --learner "$learner" --varcut "$varcut" \
  --nboot 100 --kfold 10
