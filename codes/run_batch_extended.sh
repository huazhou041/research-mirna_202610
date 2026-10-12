#!/bin/bash
# Run extended models: GLM, RF, XGB in a single job per dataset/varcut combo
# Each job runs all 3 learners on the same data split (32 total jobs)

curname='pred_mod_ext_all'
[ ! -d "logs" ] && mkdir -p logs
[ ! -d "res"  ] && mkdir -p res
datapath=rdata

for i in 5 10 20 50 100 200 500 1000
do
    # 166-sample set: Clinical, miRNA, Both
    for subset in Clinical miRNA Both
    do
        data=nyu_166_${subset}

        # Check if any learner output already exists
        if [ ! -f "./res/${data}_v${i}_GLM_cvSummary.csv" ]; then
            echo "Submitting: ${data}_v${i} (all learners)"
            sbatch --nice=1000 -J ${curname} -o logs/${data}_v${i}.log \
              codes/extended_model.sh ${datapath} ${data} ${i}
        fi
    done

    # 82-sample set: miRNA only
    data=nyu_82_miRNA

    if [ ! -f "./res/${data}_v${i}_GLM_cvSummary.csv" ]; then
        echo "Submitting: ${data}_v${i} (all learners)"
        sbatch --nice=1000 -J ${curname} -o logs/${data}_v${i}.log \
          codes/extended_model.sh ${datapath} ${data} ${i}
    fi
done

echo "Extended model jobs submitted (32 jobs: GLM/RF/XGB combined per dataset/varcut)"
