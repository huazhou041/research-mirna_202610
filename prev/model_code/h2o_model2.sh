#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load miniconda3/4.6.14
source activate r_env 

f1=/gpfs/scratch/zhouh05/machinelearning/research/passharvey/rnaseq_igor/fpkm_rm_mismatch/original/control/dataR/${1}.RDS
f3=huazhou_${1}_v${2}
[ ! -f "huazhou_${1}_rf_v${2}_cvSummary.csv" ] && ./deep-drug-response/tune_drug.R recurrence $f1 '/gpfs/scratch/zhouh05/machinelearning/research/passharvey/rnaseq_igor/fpkm_rm_mismatch/original/control/deep-drug-response/hyperparams/grid_rf.json' huazhou_${1}_rf_v${2} -v $2 --learner 'RF' -k 10 
