#!/bin/bash
curname='pred_mod'
[ ! -d "logs" ] && mkdir -p logs
datapath=rdata/
dataname=EV_RNA
for i in 5 10 20 50 100 200 500 1000
#for i in 5
	do
		for subset in Clinical miRNA Both
			do
			data=${dataname}_${subset}
			[ ! -f "./res/${data}_glm_v${i}_cvSummary.csv" ] && sbatch --nice=1000 -J ${curname} -o logs/${data}_v${i}.o_glm h2o_model3.sh ${datapath} ${data} $i
			done
	done

