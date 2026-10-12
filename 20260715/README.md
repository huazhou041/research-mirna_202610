# Plasma miRNA — Transfer & Combined-Cohort AUC Dashboard (2026-07-15)

Self-contained interactive dashboard (inlined Plotly, no internet needed), modeled
on the erdn_202607 `auc_dashboard_..._fsmethods.html`, adapted to this data's schema.

## Open
`htmls/mirna_auc_dashboard_20260715.html` — open in any browser.

## Contents (both feature-selection methods: MAD vs correlation)
- **Transfer AUC** — each cohort's whole-dataset GLM applied to the other two
  (heatmap by source→target × cutoff; line view compares MAD vs corr).
- **Combined AUC** — pooled 3-cohort model, 4 batch-handling variants
  (plain/cov/combat/both) × learner, mean CV-AUC ± 95% CI vs cutoff.
- **ROC curves** — transfer curves (source→target predictions) + combined curves
  (pooled across all fold×boot predictions); multi-select overlay.
- **Data explorer** — all AUC rows, sortable + CSV export.

## Data provenance (data/)
- `combined_{auc,roc}_long.csv` — from fork3 (MAD) + fork3_corr_fs (corr) raw_pred
  cvSummary/pred_prob (per-fold mean AUC + pooled ROC). Built by
  `codes/gather_combined_cluster.R` (run on the cluster).
- `transfer_{auc,roc}_long.csv` — whole-source GLM refit applied to other cohorts,
  both FS methods. Built locally by `codes/gather_transfer_local.R`.
- AUC uses pROC direction="<" throughout; correlation FS is fold-internal (leakage-proof).

Build: `Rscript codes/build_dashboard.R` (injects the 4 CSVs + Plotly into
`codes/dashboard_template.html`).
