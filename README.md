# research-mirna_202610

Code for the plasma miRNA cancer-detection modeling project (Harvey/Aris lab).
Only code, configs and notes are in this repo. Data, model outputs and
rendered dashboards are not included.

## Layout (roughly in chronological order)
- `prev/`: pipeline reused from the earlier EDRN project (data prep, `tune_drug.R` elastic-net/RF/DL tuning, SLURM batch scripts, AUC summaries).
- `codes/`: first pass on the new cohorts (166 / 82 samples): data prep, extended tuning, unified 3-cohort GLM, feature-importance consensus.
- `results/codes*`, `results/summary/`: analysis scripts and notes for the June results.
- `20260626_full/`: full analysis forks:
  - `fork1_results/`: archive and source→target transfer plots
  - `fork2_*`: combined models and source→target transfer, including direction-fixed rescoring
  - `fork3_cluster/`, `fork3_corr_fs_cluster/`: pooled 3-cohort models with CV-internal ComBat, MAD vs correlation feature selection (SLURM arrays)
  - `mirna_20260714/`: correlation-FS transfer runs
- `20260715/`: gathers AUC/ROC tables and builds the interactive dashboard (`codes/build_dashboard.R` + `codes/dashboard_template.html`). The dashboard build needs `plotly-2.32.0.min.js` in `20260715/codes/`, which is not included here. Download it from https://cdn.plot.ly/plotly-2.32.0.min.js.
