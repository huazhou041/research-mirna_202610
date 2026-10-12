#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Extended modeling script: runs GLM, RF, and XGB on same data split
# Usage: sbatch extended_model.sh <datapath> <dataname> <varcut>

datapath=$1
dataname=$2
varcut=$3

[ ! -d "res" ] && mkdir -p res
[ ! -d "logs" ] && mkdir -p logs

module load miniconda3/4.6.14
source activate r_env

Rscript codes/tune_drug_extended.R \
  recurrence \
  ${datapath}/${dataname}.RDS \
  codes/hyperparams.json \
  res/${dataname}_v${varcut} \
  --varcut ${varcut} \
  --kfold 10 \
  --nboot 100
