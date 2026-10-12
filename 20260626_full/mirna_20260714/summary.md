# Combined 3-Cohort miRNA Model — Feature-Selection Summary

**Date:** 2026-07-14
**Scope:** Combined EV_RNA + NYU-166 + NYU-82 plasma-miRNA recurrence classifier
(pooled 548 samples × 798 common miRNA features), batch-handling comparison
(plain / cov / combat / both) across GLM, RF, and XGB at feature cutoffs
5–1000. Root: `20260626_full/`.

---

## Two feature-selection versions

We evaluated two ways of choosing which miRNAs enter each model. In **both**,
selection is performed **inside each cross-validation training fold only**
(never on the held-out fold or the full dataset), so both are leakage-proof.

### 1. MAD selection (unsupervised) — **the version we use**
Rank miRNAs by **median absolute deviation** (variance) on the training fold and
keep the top-N. Ignores the outcome label, so it selects stable, high-variance
features that generalize across folds.

- **Canonical folder:** `fork3/`
  (CV-honest ComBat batch correction fit inside each fold; ROC `direction="<"`.)
  - Final table: `fork3/results/combined_4way_cvhonest_FINAL.csv`
- **Earlier MAD folders (superseded, kept for provenance):**
  - `fork2_new/` — original MAD run. Two issues later fixed: ComBat was fit once
    on the whole pooled dataset before CV (batch/outcome leakage), and AUC used
    pROC `direction="auto"`. Holds the original combined 4-way
    (`combined_4way_comparison.csv`) and the 3×3 cross-cohort transfer table
    (`all_sources_fixed_models_on_other2_summary.csv`).
  - `fork2_fixed_direction/` — same MAD numbers with only the ROC-direction fix
    applied (leakage not yet corrected).
- **fork3 vs fork2_new difference:** the ComBat leakage fix moved AUCs by
  ≤0.005 (within Monte-Carlo noise) because the three cohorts are balanced 1:1
  case/control, so the batch estimate barely changes when a fold is held out.
  The batch-correction advantage is real, not a leakage artifact.

**Best MAD results (fork3, CV-honest):** XGB and RF with ComBat correction peak
at **AUC ≈ 0.845–0.850 at cutoffs v20–v50**; ComBat (`combat`/`both`) beats no
correction (`plain`/`cov`) across learners.

### 2. Correlation selection (supervised) — tested, **not used**
Rank miRNAs by **|Pearson correlation with the recurrence label|** on the
training fold and keep the top-N (point-biserial; training labels only, so still
leakage-proof).

- **Folder:** `fork3_corr_fs/`
  (Identical CV-honest ComBat pipeline; only the selection criterion differs.)
  - Table: `fork3_corr_fs/results/fork3_corr_fs_4way.csv`
  - Head-to-head vs MAD: `fork3_corr_fs/results/corr_fs_vs_mad_fs_comparison.csv`
- **Result: uniformly worse than MAD**, by **−0.05 to −0.11 AUC** across every
  learner/variant, worst at v20 (≈ −0.10). The gap closes to ≈0 only at v1000
  (where "top-1000" keeps essentially all features, so selection is moot).
- **Why:** selection instability — ranking 798 features by outcome-correlation on
  ~490 training samples captures fold-specific spurious correlations that do not
  generalize to the held-out fold. MAD's variance-based, outcome-agnostic ranking
  is stable across folds. Note this is *not* leakage inflating MAD; the
  correlation selection was fully train-fold-only and simply generalizes worse.

**Conclusion:** keep **MAD selection (fork3)** as the model. Supervised
correlation filtering hurts performance here and is not used.

---

## Chromosome / sex / mitochondrial filtering — not done

We considered excluding sex-linked (chrX/chrY) and mitochondrial miRNAs from
selection but **did not implement it**, because a chromosome-based filter is not
well-defined for these miRNA features:

- miRNA probes are named by mature sequence, and the same mature miRNA is often
  transcribed from **multiple genomic loci on different chromosomes** (e.g.
  let-7f, the miR-17-92 cluster and its chrX paralog miR-106a-363). A probe
  measures a mix, so a feature cannot be cleanly assigned to one chromosome.
- Human mitochondrial miRNAs are **not annotated in miRBase** (no chrMT miRNAs),
  so a mitochondrial filter would remove nothing.

Because a feature's "chromosome" is ambiguous for miRNA, chromosome/sex filtering
was dropped rather than applied on an unreliable mapping.

---

## Folder quick-reference

| Folder | Selection | ComBat | ROC dir | Use |
|---|---|---|---|---|
| `fork2_new/` | MAD | leaky (whole-data) | auto | original (superseded) |
| `fork2_fixed_direction/` | MAD | leaky | `"<"` | direction fix only |
| **`fork3/`** | **MAD** | **CV-honest** | `"<"` | **canonical — use this** |
| `fork3_corr_fs/` | correlation | CV-honest | `"<"` | tested; worse; not used |

**Bottom line:** report from **`fork3`** (MAD, CV-honest ComBat). `fork2_new` is
the older MAD version it supersedes; `fork3_corr_fs` documents that supervised
correlation selection performs worse; chromosome filtering was intentionally not
done.

---

## Assembled results: `results/{mad,correlation}/{transfer,combined}/`

Both selection methods are provided for the two analyses, mirroring the layout
of `fork2_fixed_direction/`:

| | `transfer/` (3×3 cross-cohort) | `combined/` (4-way batch handling) |
|---|---|---|
| **`mad/`** | `all_sources_fixed_models_on_other2_summary.csv` (from fork2_fixed_direction; MAD FS, direction-fixed incl. `Source_CV_AUC`) + `source_transfer_archives/` (venn + coefficient plots) | `combined_4way_auc.csv` (fork3, CV-honest ComBat, MAD FS) |
| **`correlation/`** | `all_sources_corr_transfer_summary.csv` (correlation FS; whole-source refit applied to the other 2 cohorts; `Source_CV_AUC` = 50×10 bootstrap with fold-internal correlation FS) | `combined_4way_auc.csv` (fork3_corr_fs, CV-honest ComBat, correlation FS) |

- **Transfer** = per-cohort GLM refit on one source, applied to the other two
  (no ComBat, no pooling); `Source_CV_AUC` is the honest within-cohort CV.
- **Combined** = pooled 3-cohort model with the 4 batch-handling variants.

**Transfer finding (MAD vs correlation):** correlation-selected source features
transfer **worse** to other cohorts, especially at small cutoffs (e.g.
EV_RNA→NYU-82 is 0.25 at v5 with correlation vs 0.85 with MAD), even though the
within-cohort `Source_CV_AUC` is comparable between the two. Same conclusion as
the combined analysis: supervised correlation selection overfits the source and
does not generalize. **MAD remains the method of choice.**

Provenance code for the correlation transfer: `codes/run_corr_transfer_cell.R`
(+ cluster array), run on `mirna_202606/corr_transfer/`.
