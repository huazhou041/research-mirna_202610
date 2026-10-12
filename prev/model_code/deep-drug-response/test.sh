#!/bin/bash
for i in 1 2 5 10 20 50 100
do
./tune_drug.R recurrence '/Volumes/hua_mac/huaback/research/aris/passharvey/proj_202001/Lung RNA-seq data and clinical information.Cohort.Map.b6/theo/dataR/control_fpkm.RDS' '/Volumes/hua_mac/huaback/research/aris/passharvey/proj_202001/Lung RNA-seq data and clinical information.Cohort.Map.b6/theo/deep-drug-response/hyperparams/grid_dl.json' 'huazhou_v${i}' -v $i -c 'system' &
done
