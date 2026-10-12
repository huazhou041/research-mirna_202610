#!/usr/bin/env Rscript
# Comprehensive analysis: compare GLM, RF, XGB across datasets and varcounts
library(tidyverse)
library(pROC)
library(ggplot2)
library(gridExtra)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

# Load all cvSummary files
load_all_results = function() {
  res_dir = 'results'
  datasets = c('nyu_166_Clinical', 'nyu_166_miRNA', 'nyu_166_Both', 'nyu_82_miRNA')
  varcounts = c(5, 10, 20, 50, 100, 200, 500, 1000)
  learners = c('GLM', 'RF', 'XGB')

  all_results = data.frame()

  # Try both old format (dataset_learner_varcount) and new format (dataset_varcount_learner)
  for (dataset in datasets) {
    for (varcut in varcounts) {
      for (learner in learners) {
        # Old format: nyu_166_Both_GLM_v1000_cvSummary.csv (in raw_pred subfolder)
        file_old = file.path(res_dir, 'raw_pred', paste0(dataset, '_', learner, '_v', varcut, '_cvSummary.csv'))
        # New format: nyu_166_v1000_GLM_cvSummary.csv (in raw_pred subfolder)
        file_new = file.path(res_dir, 'raw_pred', paste0(dataset, '_v', varcut, '_', learner, '_cvSummary.csv'))

        file = if (file.exists(file_old)) file_old else if (file.exists(file_new)) file_new else NA

        if (!is.na(file)) {
          df = read.csv(file)
          df$dataset = dataset
          df$varcut = varcut
          df$learner = learner
          all_results = rbind(all_results, df)
        }
      }
    }
  }

  return(all_results)
}

cat('Loading all results...\n')
results = load_all_results()

cat('Total rows loaded:', nrow(results), '\n')
cat('Datasets:', n_distinct(results$dataset), '\n')
cat('Learners:', n_distinct(results$learner), '\n')
cat('Varcounts:', n_distinct(results$varcut), '\n')

# ============================================================
# SUMMARY TABLE: Mean AUC by dataset, learner, varcut
# ============================================================
cat('\n=== SUMMARY: Mean AUC by Dataset and Learner ===\n')

summary_table = results %>%
  group_by(dataset, learner, varcut) %>%
  summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    sd_auc = sd(AUC, na.rm = TRUE),
    min_auc = min(AUC, na.rm = TRUE),
    max_auc = max(AUC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(dataset, varcut, desc(mean_auc))

print(summary_table %>% head(30))

# ============================================================
# BEST MODEL PER DATASET
# ============================================================
cat('\n=== BEST MODEL PER DATASET (across all varcounts) ===\n')

best_models = results %>%
  group_by(dataset, learner, varcut) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = 'drop') %>%
  group_by(dataset) %>%
  slice_max(mean_auc, n = 3)

print(best_models)

# ============================================================
# LEARNER COMPARISON BY DATASET
# ============================================================
cat('\n=== LEARNER COMPARISON (aggregated across all varcounts) ===\n')

learner_agg = results %>%
  group_by(dataset, learner) %>%
  summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    sd_auc = sd(AUC, na.rm = TRUE),
    n_points = n(),
    .groups = 'drop'
  ) %>%
  arrange(dataset, desc(mean_auc))

print(learner_agg)

# ============================================================
# PLOT 1: AUC by varcut for each dataset/learner combo
# ============================================================
cat('\nGenerating plots...\n')

p1 = results %>%
  group_by(dataset, learner, varcut) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = varcut, y = mean_auc, color = learner, shape = learner)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.6) +
  facet_wrap(~dataset, scales = 'free_y', nrow = 2) +
  labs(
    title = 'Mean AUC by Feature Count and Learner',
    x = 'Number of Features (MAD-selected)',
    y = 'Mean AUC',
    color = 'Learner',
    shape = 'Learner'
  ) +
  scale_x_log10() +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    strip.text = element_text(face = 'bold')
  )

ggsave(file.path('results/plots', 'learner_comparison_by_varcut.pdf'), p1, width = 12, height = 8)
cat('Saved: learner_comparison_by_varcut.pdf\n')

# ============================================================
# PLOT 2: Distribution of AUC across bootstrap iterations
# ============================================================
p2 = results %>%
  ggplot(aes(x = learner, y = AUC, fill = learner)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~dataset, nrow = 2) +
  labs(
    title = 'Distribution of AUC Across Bootstrap Iterations',
    x = 'Learner',
    y = 'AUC',
    fill = 'Learner'
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    strip.text = element_text(face = 'bold')
  )

ggsave(file.path('results/plots', 'auc_distribution_by_learner.pdf'), p2, width = 12, height = 8)
cat('Saved: auc_distribution_by_learner.pdf\n')

# ============================================================
# PLOT 3: Heatmap of mean AUC
# ============================================================
heatmap_data = results %>%
  group_by(dataset, learner, varcut) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = learner, values_from = mean_auc)

p3 = results %>%
  group_by(dataset, learner, varcut) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = factor(varcut), y = fct_reorder(dataset, mean_auc, .desc = TRUE), fill = mean_auc)) +
  geom_tile(color = 'white', size = 0.5) +
  facet_wrap(~learner, nrow = 1) +
  scale_fill_viridis_c(name = 'Mean AUC') +
  labs(
    title = 'Mean AUC Heatmap: Dataset × Varcount × Learner',
    x = 'Feature Count',
    y = 'Dataset'
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = 'bold')
  )

ggsave(file.path('results/plots', 'auc_heatmap.pdf'), p3, width = 14, height = 6)
cat('Saved: auc_heatmap.pdf\n')

# ============================================================
# STATISTICS: Test if learners differ significantly
# ============================================================
cat('\n=== STATISTICAL COMPARISON ===\n')

for (ds in unique(results$dataset)) {
  cat('\nDataset:', ds, '\n')
  ds_data = filter(results, dataset == ds)

  # ANOVA across learners
  if (length(unique(ds_data$learner)) == 3) {
    aov_result = aov(AUC ~ learner, data = ds_data)
    f_stat = summary(aov_result)[[1]]$F[1]
    p_val = summary(aov_result)[[1]]$Pr[1]
    cat(sprintf('  ANOVA: F=%5.2f, p=%0.4f\n', f_stat, p_val))
  }

  # Pairwise learner means
  learner_means = ds_data %>%
    group_by(learner) %>%
    summarise(mean_auc = mean(AUC, na.rm = TRUE), .groups = 'drop')
  print(learner_means)
}

# ============================================================
# OUTPUT SUMMARY CSV
# ============================================================
write.csv(summary_table, 'results/all_learners_summary.csv', row.names = FALSE)
write.csv(learner_agg, 'results/learner_aggregated.csv', row.names = FALSE)

cat('\n=== Analysis complete ===\n')
cat('Summary saved to:\n')
cat('  - results/all_learners_summary.csv\n')
cat('  - results/learner_aggregated.csv\n')
cat('  - results/plots/*.pdf\n')
