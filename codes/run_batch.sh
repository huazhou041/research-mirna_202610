#!/bin/bash
curname='pred_mod'
[ ! -d "logs" ] && mkdir -p logs
[ ! -d "res"  ] && mkdir -p res
datapath=rdata

for i in 5 10 20 50 100 200 500 1000
do
    # 166-sample set: Clinical, miRNA, Both
    for subset in Clinical miRNA Both
    do
        data=nyu_166_${subset}
        [ ! -f "./res/${data}_glm_v${i}_cvSummary.csv" ] && \
            sbatch --nice=1000 -J ${curname} -o logs/${data}_v${i}.o_glm \
            h2o_model3.sh ${datapath} ${data} ${i}
    done

    # 82-sample set: miRNA only
    data=nyu_82_miRNA
    [ ! -f "./res/${data}_glm_v${i}_cvSummary.csv" ] && \
        sbatch --nice=1000 -J ${curname} -o logs/${data}_v${i}.o_glm \
        h2o_model3.sh ${datapath} ${data} ${i}
done
