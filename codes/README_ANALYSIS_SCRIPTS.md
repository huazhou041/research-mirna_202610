# Analysis Scripts for EDRN Model Comparison

## Overview
These scripts generate visualizations and tables comparing EDRN models with nyu_166 and nyu_82 datasets.

---

## Scripts

### 1. `edrn_all_models_analysis.R`
**Purpose:** Apply all EDRN models (v5-v1000) directly to new datasets and generate summary table

**Output:**
- `results/summary/edrn_all_models_summary.csv` — Main results table
- `results/summary/edrn_v5_features.csv` through `edrn_v1000_features.csv` — Individual model features & coefficients
- `results/summary/edrn_all_features_combined.csv` — All unique features across models
- `results/summary/edrn_consensus_features_3plus.csv` — Consensus features (appearing in 3+ models)

**Key Output: edrn_all_models_summary.csv columns:**
- Version: EDRN model version (v5, v10, v20, v50, v100, v200, v500, v1000)
- N_Features: Number of features in that model
- EDRN_AUC: Original EDRN model AUC
- AUC_166: When applied to nyu_166
- Delta_166: Performance change for nyu_166
- AUC_82: When applied to nyu_82
- Delta_82: Performance change for nyu_82

---

### 2. `regenerate_venn_importance.R`
**Purpose:** Create Venn diagrams showing feature overlap across 3 models

**Output:**
- `results/plots/venn_top20_3models.pdf` — 3-way Venn for top 20 features
- `results/plots/venn_top50_3models.pdf` — 3-way Venn for top 50 features

**Compares:**
- nyu_166 top 20/50 (by MAD)
- nyu_82 top 20/50 (by MAD)
- EDRN v20/v50 (20/50 features from trained models)

---

### 3. `feature_importance_consensus.R`
**Purpose:** Compare feature coefficients across models for consensus features

**Output:**
- `results/plots/feature_importance_top20_3models.pdf` — Consensus features (7 features, 3 bars each)
- `results/plots/feature_importance_top50_3models.pdf` — Consensus features (18 features, 3 bars each)

**Format:** 
- Each feature has 3 bars showing coefficients from:
  - Red bar = nyu_166 GLM coefficient
  - Blue bar = nyu_82 GLM coefficient
  - Green bar = EDRN GLM coefficient
- Only shows features appearing in ALL 3 models

---

## How to Run

Run each script individually:

```bash
# Generate EDRN model summary table & individual feature files
Rscript codes/edrn_all_models_analysis.R

# Generate Venn diagrams
Rscript codes/regenerate_venn_importance.R

# Generate feature importance comparison plots
Rscript codes/feature_importance_consensus.R
```

Or run all three:
```bash
Rscript codes/edrn_all_models_analysis.R && \
Rscript codes/regenerate_venn_importance.R && \
Rscript codes/feature_importance_consensus.R
```

---

## Key Findings Summary

### EDRN Model Transfer Performance
- **Best transfer:** EDRN v5 (0.827 AUC on nyu_82)
- **Worst transfer:** EDRN v500-v1000 (severe overfitting)
- **Pattern:** Smaller models generalize better

### Consensus Features
- **Top 20 analysis:** 7 features appear in all 3 models
- **Top 50 analysis:** 18 features appear in all 3 models
- These consensus features represent the most robust biomarkers

### Model Comparisons
- EDRN transfers well to nyu_82 (similar population)
- EDRN transfers poorly to nyu_166 (different disease pattern)
- nyu_82 (screening) more similar to EDRN (healthy vs cancer)
- nyu_166 (benign surgical vs cancer surgical) is harder task

---

## File Locations

All outputs saved to:
- **PDFs:** `/results/plots/`
- **CSVs:** `/results/summary/`

Main reference files:
- `edrn_all_models_summary.csv` — Start here for performance comparison
- `venn_top50_3models.pdf` — Feature overlap visualization
- `feature_importance_top50_3models.pdf` — Coefficient comparison for consensus features

