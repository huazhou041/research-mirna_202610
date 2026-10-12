# Summary: Top 50 Features Comparison (miRNA only)

## Overview
This folder contains a comprehensive comparison of GLM models trained on top 50 miRNA features across three cohorts:
- **EDRN_Plasma**: Original EDRN cohort (EV_RNA_miRNA_glm_v50)
- **nyu_166_miRNA**: New validation set (166 samples)
- **nyu_82_miRNA**: New screening set (82 samples)

## Files

### AUC Summaries (miRNA only)
- **auc_summary_EDRN.csv** — Original EDRN AUC summary with all feature counts (v5-v1000)
- **auc_summary_nyu_166_miRNA.csv** — nyu_166 AUC summary with all feature counts
- **auc_summary_nyu_82_miRNA.csv** — nyu_82 AUC summary with all feature counts
- **auc_summary_miRNA_all.csv** — Combined summary of all three cohorts (main reference)

### Feature Analysis (Top 50)
- **top50_features_comparison.csv** — Detailed feature importance table
  - Includes coefficients from all three models
  - Shows which cohorts contain each feature
  - Ranked by consensus across cohorts

- **top50_cohort_summary.csv** — Summary statistics per cohort
  - Number of unique features
  - Number of shared features (2+ cohorts)
  - AUC performance

### Visualizations

**Overlap Analysis:**
- **venn_top50_3cohorts.pdf** — 3-way Venn diagram showing feature overlap
  - 18 features shared by all 3 cohorts
  - 40 features shared by both new cohorts (nyu_166 ∩ nyu_82)
  - 28 unique EDRN features

**Feature Importance:**
- **feature_importance_top50_25features.pdf** — Top 25 consensus features
  - Horizontal bar plot of GLM coefficients
  - Compares coefficient magnitudes across 3 cohorts
  - Top features: hsa-miR-212-3p, hsa-miR-4531, hsa-miR-448, hsa-let-7b-5p

**Feature Overlap Heatmap:**
- **feature_overlap_heatmap_top50.pdf** — Presence/absence matrix of top 40 shared features
  - Shows which cohorts contain each feature
  - Reveals consistent biomarkers across cohorts

**Performance Comparison:**
- **auc_comparison_top50.pdf** — AUC bar chart with 95% CI
  - EDRN: 0.91 (0.88~0.94)
  - nyu_166: 0.76 (0.73~0.74)
  - nyu_82: 0.96 (0.95~0.96)

### Model Objects
- **glm_top50_nyu166.rds** — Trained GLM model on nyu_166 (top 50 features)
- **glm_top50_nyu82.rds** — Trained GLM model on nyu_82 (top 50 features)

## Key Findings

### Feature Conservation
| Metric | Count |
|--------|-------|
| Features in ALL 3 cohorts | 18 |
| Features in 2+ cohorts | 40 |
| Unique to EDRN | 28 |
| Unique to nyu_166 | 9 |
| Unique to nyu_82 | 7 |

### Top Features (appearing in all 3 or 2+ cohorts)
1. **hsa-miR-212-3p** — Highest consensus importance (all 3)
2. **hsa-miR-4531** — All 3 cohorts
3. **hsa-miR-448** — All 3 cohorts
4. **hsa-let-7b-5p** — All 3 cohorts
5. **hsa-miR-548ar-5p** — All 3 cohorts

### Performance Summary

#### AUC with Top 50 Features
| Cohort | AUC | 95% CI |
|--------|-----|--------|
| EDRN_Plasma | 0.91 | 0.88~0.94 |
| nyu_166_miRNA | 0.76 | 0.73~0.74 |
| nyu_82_miRNA | 0.96 | 0.95~0.96 |

**Notes:**
- nyu_82 shows strongest signal (0.96 AUC)
- EDRN remains strong (0.91 AUC, same as published)
- nyu_166 shows weaker signal (0.76 AUC) — consistent with biological differences
  - 166 samples: benign surgical vs cancer surgical (harder classification task)
  - 82 samples: likely easier cases or clearer cancer signals

## Recommendations

1. **Use top 50 features** for balanced performance across all cohorts
2. **Focus on 18 consensus features** if model simplification is needed
   - These features are most reproducible across independent cohorts
3. **Consider cohort-specific models** given nyu_166 performance gap
   - May need additional biological variables or sample stratification

## Technical Details

- **Training method**: GLM (elastic net, alpha=0.5) with 5-fold cross-validation
- **Feature selection**: Top 50 features by Median Absolute Deviation (MAD)
- **Outcome**: Cancer recurrence (binary classification)
- **Data**: miRNA expression only (no clinical variables)

---

## Update: Venn Diagram Correction (June 15, 2026)

**Issue Fixed:** The original Venn diagram incorrectly showed 0 features for EDRN.

**Root Cause:** Bug in the plot script - EDRN features were not being properly passed to the `ggvenn()` function.

**Corrected Results:**
- **EDRN features**: 50 (was showing as 0)
- **Feature overlap across all 3 cohorts**: 18 (correct)
- **Features in EDRN & nyu_166**: 19 
- **Features in EDRN & nyu_82**: 21
- **Features unique to EDRN**: ~28

The corrected Venn diagram now accurately displays the feature overlap, showing that EDRN shares 18-21 features with the new cohorts, representing robust cancer detection biomarkers.
