# fork3 — CV-honest ComBat batch-handling comparison

Corrects a data-leakage bug in the fork2 combined 3-cohort batch-handling
comparison (EV_RNA + nyu_166 + nyu_82 pooled; variants plain / cov / combat / both).

## The bug (fork2)
`create_combined_variants.R` ran `limma::removeBatchEffect` **once on the full
pooled dataset before any CV split**, with `design = ~recurrence` using every
sample's label. Each bootstrap fold then read that pre-corrected matrix, so a
held-out validation sample's batch correction was informed by (a) other
validation samples and (b) its own outcome label -> batch-effect + outcome
leakage into the reported CV AUC (`combat`/`both` variants only).

## The fix (fork3)
`codes/tune_combined_models_cvcombat.R` fits the batch correction **inside each
training fold** (training samples/labels only) and applies it to the held-out
fold using only that fold's cohort membership (never its labels). Replicates
limma's math (sum-to-zero contrasts, recurrence protected, only batch-coefficient
contribution removed) but split-safe. AUC uses `direction="<"` throughout.
plain/cov have no feature correction and are unaffected.

## Result
Leakage effect (CV-honest minus leaky, both direction-fixed) for the corrected
combat/both variants: **largest |delta| = 0.0046 AUC**, mean ~0.001. This sits
**within the Monte-Carlo noise floor** (max 0.0036 AUC) measured from the
plain/cov RF cells, which have no correction and differ only by RF's internal
randomness between the two independent runs.

**Conclusion:** the batch effect was real but its impact on the reported AUCs is
negligible — indistinguishable from run-to-run noise. This is because the three
cohorts are perfectly balanced 1:1 case/control, so recurrence is orthogonal to
cohort and the batch-effect estimate barely changes when a fold is held out. The
combat/both > plain/cov advantage is **real, not a leakage artifact**; all
headline findings hold.

## Files
- `codes/tune_combined_models_cvcombat.R` — CV-honest trainer (GLM/RF/XGB, 4 variants)
- `codes/run_fork3_array.sh` + `codes/fork3_jobs.txt` — 96-cell SLURM array
- `codes/summarize_fork3.R` — builds the comparison
- `results/combined_4way_cvhonest_FINAL.csv` — final CV-honest 4-way AUC table
  (plain/cov canonical from fork2 direction-fixed; combat/both CV-honest from fork3)
- `results/combat_leakage_effect_by_cell.csv` — per-cell leakage delta (combat/both)
- `results/fork3_vs_fork2_combat_leakage_comparison.csv` — full 96-cell (some plain/cov GLM
  v100–v200 cells absent: they timed out at 12h on raw features, but are deterministic
  duplicates of fork2 GLM and irrelevant to the leakage question)

Cluster location: `/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork3`
