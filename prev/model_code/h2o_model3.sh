#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00

module load miniconda3/4.6.14
source activate r_env 

f1=./${1}/${2}.RDS
[ ! -d "res" ] && mkdir -p res
[ ! -f "$res/${2}_glm_v${3}_cvSummary.csv" ] && ./deep-drug-response/tune_drug.R recurrence $f1 './deep-drug-response/hyperparams/grid_glm.json' ./res/${2}_glm_v${3} -v $3 --learner 'GLM' -k 10 --nboot 100 
