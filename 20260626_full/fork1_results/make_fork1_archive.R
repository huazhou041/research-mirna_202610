suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggvenn)
})

root <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605"
base_dir <- file.path(root, "20260626_full/fork1_results")
raw_dir <- file.path(base_dir, "raw")
archive_dir <- file.path(base_dir, "Archive")
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)

var_levels <- c(5, 10, 20, 50, 100, 200, 500, 1000)
learner_levels <- c("GLM", "RF", "XGB")

theme_archive <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

save_pdf <- function(path, plot, width = 7.5, height = 5) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
}

# Transfer-back: NYU 166/82 GLM models applied to EV_RNA.
transfer <- read_csv(
  file.path(raw_dir, "transfer_results/nyu_glm_models_on_EV_RNA_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    varcut = factor(varcut, levels = var_levels),
    source_model = recode(source_model,
      "nyu_166_miRNA" = "NYU 166 GLM",
      "nyu_82_miRNA" = "NYU 82 GLM"
    )
  ) %>%
  arrange(source_model, as.integer(as.character(varcut)))

transfer_best <- transfer %>%
  group_by(source_model) %>%
  slice_max(auc, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(transfer, file.path(archive_dir, "transfer_nyu_models_on_EVRNA_summary.csv"))
write_csv(transfer_best, file.path(archive_dir, "transfer_nyu_models_on_EVRNA_best.csv"))

p_transfer_auc <- ggplot(transfer, aes(varcut, auc, color = source_model, group = source_model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_point(data = transfer_best, aes(varcut, auc), size = 3.8, shape = 21, stroke = 1.1, fill = "white") +
  scale_y_continuous(limits = c(0.45, max(transfer$auc, na.rm = TRUE) + 0.04)) +
  labs(x = "Feature cutoff", y = "AUC on EV_RNA", color = "Source model") +
  theme_archive()
save_pdf(file.path(archive_dir, "transfer_auc_by_varcut.pdf"), p_transfer_auc)

p_transfer_nonzero <- ggplot(transfer, aes(varcut, n_nonzero, color = source_model, group = source_model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  labs(x = "Feature cutoff", y = "Nonzero GLM coefficients", color = "Source model") +
  theme_archive()
save_pdf(file.path(archive_dir, "transfer_nonzero_features_by_varcut.pdf"), p_transfer_nonzero)

# Combined 3-cohort retraining CV summaries.
cv_files <- list.files(
  file.path(raw_dir, "raw_pred"),
  pattern = "^combined_3cohort_miRNA_(GLM|RF|XGB)_v[0-9]+_cvSummary\\.csv$",
  full.names = TRUE
)

combined_cv <- map_dfr(cv_files, function(path) {
  b <- basename(path)
  tibble(
    learner = str_match(b, "^combined_3cohort_miRNA_([A-Z]+)_v")[, 2],
    varcut = as.integer(str_match(b, "_v([0-9]+)_cvSummary\\.csv$")[, 2])
  ) %>%
    bind_cols(read_csv(path, show_col_types = FALSE))
}) %>%
  mutate(
    learner = factor(learner, levels = learner_levels),
    varcut = factor(varcut, levels = var_levels)
  )

combined_summary <- combined_cv %>%
  group_by(learner, varcut) %>%
  summarise(
    n = n(),
    mean_auc = mean(AUC, na.rm = TRUE),
    sd_auc = sd(AUC, na.rm = TRUE),
    median_auc = median(AUC, na.rm = TRUE),
    q25_auc = quantile(AUC, 0.25, na.rm = TRUE),
    q75_auc = quantile(AUC, 0.75, na.rm = TRUE),
    min_auc = min(AUC, na.rm = TRUE),
    max_auc = max(AUC, na.rm = TRUE),
    .groups = "drop"
  )

metadata <- read_csv(
  file.path(raw_dir, "analysis_outputs/combined_3cohort_model_metadata.csv"),
  show_col_types = FALSE
) %>%
  mutate(learner = factor(learner, levels = learner_levels), varcut = factor(varcut, levels = var_levels))

combined_summary <- combined_summary %>%
  left_join(metadata, by = c("learner", "varcut")) %>%
  arrange(learner, as.integer(as.character(varcut)))

combined_best <- combined_summary %>%
  group_by(learner) %>%
  slice_max(mean_auc, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(combined_cv, file.path(archive_dir, "combined_3cohort_model_cv_auc_long.csv"))
write_csv(combined_summary, file.path(archive_dir, "combined_3cohort_model_summary.csv"))
write_csv(combined_best, file.path(archive_dir, "combined_3cohort_model_best_by_learner.csv"))

p_combined_line <- ggplot(combined_summary, aes(varcut, mean_auc, color = learner, group = learner)) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), width = 0.15, alpha = 0.65) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_point(data = combined_best, aes(varcut, mean_auc), size = 3.8, shape = 21, stroke = 1.1, fill = "white") +
  labs(x = "Feature cutoff", y = "Mean CV AUC", color = "Learner") +
  theme_archive()
save_pdf(file.path(archive_dir, "combined_cv_auc_by_varcut.pdf"), p_combined_line)

p_combined_box <- ggplot(combined_cv, aes(varcut, AUC, fill = learner)) +
  geom_boxplot(outlier.alpha = 0.35, width = 0.72) +
  facet_wrap(~ learner, nrow = 1) +
  guides(fill = "none") +
  labs(x = "Feature cutoff", y = "CV fold AUC") +
  theme_archive()
save_pdf(file.path(archive_dir, "combined_cv_auc_boxplot.pdf"), p_combined_box, width = 9, height = 4.5)

p_combined_heat <- ggplot(combined_summary, aes(varcut, learner, fill = mean_auc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean_auc)), size = 3.2) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(x = "Feature cutoff", y = NULL, fill = "Mean CV AUC") +
  theme_archive() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(file.path(archive_dir, "combined_cv_auc_heatmap.pdf"), p_combined_heat, width = 7.5, height = 3.8)

# Feature importance and Venn plots for the best varcut per learner.
importance <- read_csv(
  file.path(raw_dir, "analysis_outputs/combined_3cohort_model_feature_importance.csv"),
  show_col_types = FALSE
) %>%
  mutate(learner = factor(learner, levels = learner_levels), varcut = factor(varcut, levels = var_levels))

best_keys <- combined_best %>% select(learner, varcut)
importance_best <- importance %>%
  inner_join(best_keys, by = c("learner", "varcut")) %>%
  filter(is.finite(importance), importance > 0) %>%
  group_by(learner) %>%
  arrange(desc(importance), feature, .by_group = TRUE) %>%
  mutate(best_rank = row_number()) %>%
  ungroup()

write_csv(importance, file.path(archive_dir, "combined_3cohort_feature_importance_all.csv"))
write_csv(importance_best, file.path(archive_dir, "combined_3cohort_feature_importance_best_models.csv"))

make_importance_plot <- function(top_n) {
  top_tbl <- importance_best %>%
    group_by(learner) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(feature = fct_reorder(feature, importance))

  p <- ggplot(top_tbl, aes(importance, feature, fill = learner)) +
    geom_col(width = 0.72) +
    facet_wrap(~ learner, scales = "free_y", nrow = 1) +
    guides(fill = "none") +
    labs(x = "Model-specific importance", y = NULL) +
    theme_archive() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  write_csv(top_tbl, file.path(archive_dir, sprintf("feature_importance_top%d_3models.csv", top_n)))
  save_pdf(file.path(archive_dir, sprintf("feature_importance_top%d_3models.pdf", top_n)), p, width = 12, height = ifelse(top_n == 20, 7, 11))
}

make_venn <- function(top_n) {
  sets <- importance_best %>%
    group_by(learner) %>%
    slice_head(n = top_n) %>%
    summarise(features = list(feature), .groups = "drop") %>%
    deframe()

  p <- ggvenn(
    sets,
    fill_color = c("#3B7EA1", "#C76F41", "#5B8F58"),
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 3.5
  ) +
    labs(title = sprintf("Top %d Features From Best Combined Models", top_n)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  overlap_tbl <- tibble(feature = sort(unique(unlist(sets)))) %>%
    mutate(
      in_GLM = feature %in% sets[["GLM"]],
      in_RF = feature %in% sets[["RF"]],
      in_XGB = feature %in% sets[["XGB"]],
      n_models = in_GLM + in_RF + in_XGB
    ) %>%
    arrange(desc(n_models), feature)

  write_csv(overlap_tbl, file.path(archive_dir, sprintf("venn_top%d_3models_features.csv", top_n)))
  save_pdf(file.path(archive_dir, sprintf("venn_top%d_3models.pdf", top_n)), p, width = 6.5, height = 5.5)
}

walk(c(20, 50), make_importance_plot)
walk(c(20, 50), make_venn)

cat("Wrote fork1 Archive outputs to:", archive_dir, "\n")
cat("Transfer summaries:", nrow(transfer), "rows\n")
cat("Combined CV rows:", nrow(combined_cv), "\n")
cat("Combined model summaries:", nrow(combined_summary), "rows\n")
