#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=corrxfer
#SBATCH --output=logs/corrxfer_%A_%a.log
#SBATCH --array=1-24
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/corr_transfer
mkdir -p cells logs
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" codes/corr_transfer_jobs.txt)
s=$(echo "$line"|awk '{print $1}'); v=$(echo "$line"|awk '{print $2}')
echo "cell: $s v$v (ntasks=$SLURM_NTASKS)"
Rscript codes/run_corr_transfer_cell.R "$s" "$v"
