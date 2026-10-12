library(tidyverse)
library(glmnet)
library(pROC)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

cat('=== ALL EDRN MODELS: COMPLETE FEATURE & COEFFICIENT ANALYSIS ===\n\n')

out_dir = '20260626_full/results'
summary_dir = file.path(out_dir, 'summary')
model_dir = '20260626_full/models'
dir.create(summary_dir, recursive=TRUE, showWarnings=FALSE)

# Load new data
d_166 = readRDS('data/model_data/nyu_166_miRNA.RDS')
d_82 = readRDS('data/model_data/nyu_82_miRNA.RDS')

d_166_clean = d_166[!is.na(d_166$recurrence), ]
d_82_clean = d_82[!is.na(d_82$recurrence), ]

varlist = c(5, 10, 20, 50, 100, 200, 500, 1000)

# Load EDRN AUC summary
auc_edrn = read.csv('results/summary/auc_summary_EDRN.csv') %>%
  filter(dataname == 'miRNA')

# Store all model info
all_models_info = list()

for (v in varlist) {
  cat(sprintf('Processing EDRN v%d...\n', v))

  model_file = file.path(model_dir, sprintf('EV_RNA_miRNA_glm_v%d.rds', v))

  if (!file.exists(model_file)) {
    cat(sprintf('  Model file not found: %s\n', model_file))
    next
  }

  prev_model = readRDS(model_file)
  coef_all = coef(prev_model, s='lambda.min')

  # Extract features and coefficients
  model_features = rownames(coef_all)[-1]  # exclude intercept
  coefs = as.vector(coef_all[-1])
  names(coefs) = model_features
  active_features = names(coefs)[coefs != 0]

  # Get EDRN original AUC
  col_name = paste0('var_', v)
  edrn_auc_str = auc_edrn[[col_name]]
  edrn_auc = as.numeric(sub(' .*$', '', edrn_auc_str))

  # Apply to nyu_166
  X_166 = as.matrix(d_166_clean[, model_features, drop=FALSE])
  pred_166 = predict(prev_model, newx=X_166, s='lambda.min', type='response')
  auc_166 = roc(d_166_clean$recurrence, as.vector(pred_166), quiet=TRUE)$auc

  # Apply to nyu_82
  X_82 = as.matrix(d_82_clean[, model_features, drop=FALSE])
  pred_82 = predict(prev_model, newx=X_82, s='lambda.min', type='response')
  auc_82 = roc(d_82_clean$recurrence, as.vector(pred_82), quiet=TRUE)$auc

  all_models_info[[as.character(v)]] = list(
    v = v,
    model_features = model_features,
    active_features = active_features,
    coefs = coefs,
    edrn_auc = edrn_auc,
    auc_166 = auc_166,
    auc_82 = auc_82,
    delta_166 = auc_166 - edrn_auc,
    delta_82 = auc_82 - edrn_auc,
    n_features = length(model_features),
    n_nonzero = length(active_features)
  )
}

# ===== CREATE SUMMARY TABLE =====
cat('\n\n=== SUMMARY TABLE: ALL EDRN MODELS ===\n\n')

summary_df = data.frame(
  Version = as.numeric(names(all_models_info)),
  N_Features = sapply(all_models_info, function(x) x$n_features),
  N_Nonzero = sapply(all_models_info, function(x) x$n_nonzero),
  EDRN_AUC = sapply(all_models_info, function(x) x$edrn_auc),
  AUC_166 = sapply(all_models_info, function(x) x$auc_166),
  Delta_166 = sapply(all_models_info, function(x) x$delta_166),
  AUC_82 = sapply(all_models_info, function(x) x$auc_82),
  Delta_82 = sapply(all_models_info, function(x) x$delta_82)
)

print(summary_df)

write.csv(summary_df, file.path(summary_dir, 'edrn_all_models_summary.csv'),
          row.names=FALSE, quote=FALSE)
cat(sprintf('\nSaved: %s\n', file.path(summary_dir, 'edrn_all_models_summary.csv')))

# ===== CREATE DETAILED FEATURE TABLES FOR EACH MODEL =====
cat('\n\n=== DETAILED FEATURE TABLES ===\n\n')

for (v in varlist) {
  v_str = as.character(v)
  if (!v_str %in% names(all_models_info)) next

  info = all_models_info[[v_str]]

  # Create feature table
  feat_df = data.frame(
    Version = v,
    Feature = info$model_features,
    Coefficient = info$coefs,
    AbsCoef = abs(info$coefs),
    Nonzero = info$coefs != 0,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(AbsCoef))

  # Save individual feature file
  filename = file.path(summary_dir, sprintf('edrn_v%d_features.csv', v))
  write.csv(feat_df, filename, row.names=FALSE, quote=FALSE)

  cat(sprintf('EDRN v%d: %d input features, %d nonzero coefficients\n',
              v, length(info$model_features), length(info$active_features)))
  cat(sprintf('  EDRN AUC: %.3f\n', info$edrn_auc))
  cat(sprintf('  Applied to 166: %.3f (Delta %+.3f)\n', info$auc_166, info$delta_166))
  cat(sprintf('  Applied to 82:  %.3f (Delta %+.3f)\n', info$auc_82, info$delta_82))
  cat(sprintf('  Top 5 features by coefficient:\n'))

  for (i in 1:min(5, nrow(feat_df))) {
    f = feat_df$Feature[i]
    c = feat_df$Coefficient[i]
    cat(sprintf('    %2d. %-50s coef: %+.6f\n', i, f, c))
  }
  cat('\n')
}

# ===== COMBINED FEATURE TABLE (all features from all models) =====
cat('\n=== COMBINED FEATURE ANALYSIS ===\n')

all_features = unique(unlist(lapply(all_models_info, function(x) x$model_features)))

combined_feat_df = data.frame(
  Feature = all_features,
  stringsAsFactors = FALSE
)

# Add columns for each model
for (v in varlist) {
  v_str = as.character(v)
  if (!v_str %in% names(all_models_info)) next

  info = all_models_info[[v_str]]
  coef_vec = info$coefs

  # Create column with coefficient or NA
  col_name = sprintf('v%d_coef', v)
  combined_feat_df[[col_name]] = sapply(all_features, function(f) {
    if (f %in% names(coef_vec)) coef_vec[f] else NA
  })
}

# Count how many models contain each feature
combined_feat_df$N_Models = rowSums(!is.na(combined_feat_df[, -1]))
combined_feat_df$N_Nonzero_Models = rowSums(combined_feat_df[, -c(1, ncol(combined_feat_df))] != 0, na.rm=TRUE)
combined_feat_df$Mean_AbsCoef = rowMeans(abs(combined_feat_df[, grep('_coef$', names(combined_feat_df))]), na.rm=TRUE)

combined_feat_df = combined_feat_df %>%
  arrange(desc(N_Nonzero_Models), desc(N_Models), desc(Mean_AbsCoef))

# Save
write.csv(combined_feat_df, file.path(summary_dir, 'edrn_all_features_combined.csv'),
          row.names=FALSE, quote=FALSE)

cat('Saved: edrn_all_features_combined.csv\n')
cat('Total unique features across all EDRN models:', nrow(combined_feat_df), '\n')
cat('Features appearing in multiple models:\n')

multi_model_feat = combined_feat_df %>%
  filter(N_Nonzero_Models >= 2) %>%
  select(Feature, N_Models, N_Nonzero_Models, Mean_AbsCoef) %>%
  head(20)

print(multi_model_feat)

# ===== TOP CONSENSUS FEATURES (appearing in 3+ models) =====
cat('\n\n=== TOP CONSENSUS EDRN FEATURES (3+ models) ===\n\n')

consensus_edrn = combined_feat_df %>%
  filter(N_Nonzero_Models >= 3) %>%
  select(Feature, N_Models, N_Nonzero_Models, Mean_AbsCoef) %>%
  arrange(desc(Mean_AbsCoef))

print(consensus_edrn)

write.csv(consensus_edrn, file.path(summary_dir, 'edrn_consensus_features_3plus.csv'),
          row.names=FALSE, quote=FALSE)

cat('\nSaved: edrn_consensus_features_3plus.csv\n')
