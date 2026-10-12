suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggvenn)
  library(glmnet)
})

root <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605"
base_dir <- file.path(root, "20260626_full/fork2_results")
archive_root <- file.path(base_dir, "source_transfer_archives")
model_dir <- file.path(base_dir, "models")
sources <- c("EV_RNA", "nyu_166", "nyu_82")
varcuts <- c(5, 10, 20, 50, 100, 200, 500, 1000)

pretty_name <- function(x) {
  recode(x, EV_RNA = "EV_RNA", nyu_166 = "NYU 166", nyu_82 = "NYU 82")
}

model_path <- function(cohort, v) {
  file.path(model_dir, sprintf("%s_miRNA_glm_v%d_fixed_whole_features.rds", cohort, v))
}

extract_coefs <- function(model) {
  cm <- coef(model, s = "lambda.min")
  coefs <- as.numeric(cm[-1, 1])
  names(coefs) <- rownames(cm)[-1]
  coefs
}

theme_archive <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

save_pdf <- function(path, plot, width = 8, height = 5) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
}

make_one_source <- function(source) {
  targets <- setdiff(sources, source)
  archive_dir <- file.path(archive_root, paste0(source, "_source"), "Archive")
  summary_dir <- file.path(archive_root, paste0(source, "_source"), "summary")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  summary_file <- file.path(archive_dir, sprintf("%s_models_on_other2_summary.csv", source))
  summary_df <- read_csv(summary_file, show_col_types = FALSE)

  long_auc <- summary_df %>%
    select(Version, Source_CV_AUC, starts_with("AUC_")) %>%
    pivot_longer(-Version, names_to = "Cohort", values_to = "AUC") %>%
    mutate(
      Cohort = case_when(
        Cohort == "Source_CV_AUC" ~ paste0(pretty_name(source), " source CV"),
        TRUE ~ paste0("Applied to ", pretty_name(sub("^AUC_", "", Cohort)))
      ),
      Version = factor(Version, levels = varcuts)
    )

  p_auc <- ggplot(long_auc, aes(Version, AUC, color = Cohort, group = Cohort)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    labs(x = "Feature cutoff", y = "AUC", color = NULL) +
    theme_archive()
  save_pdf(file.path(archive_dir, "auc_by_varcut.pdf"), p_auc)

  make_venn_importance <- function(v, top_n, height) {
    cohorts <- c(source, targets)
    coefs_by_cohort <- setNames(lapply(cohorts, function(cohort) extract_coefs(readRDS(model_path(cohort, v)))), cohorts)
    sets <- lapply(coefs_by_cohort, names)
    names(sets) <- paste0(pretty_name(names(sets)), "\nv", v)

    p_venn <- ggvenn(
      sets,
      fill_color = c("#377EB8", "#E41A1C", "#4DAF4A"),
      stroke_size = 2.2,
      set_name_size = 5,
      text_size = 4
    ) +
      labs(
        title = sprintf("Feature Overlap: v%d Fixed GLM Models", v),
        subtitle = paste(pretty_name(source), "fixed source model and the other two cohorts")
      ) +
      theme_bw() +
      theme(plot.title = element_text(face = "bold", size = 13))
    save_pdf(file.path(archive_dir, sprintf("venn_top%d_3models.pdf", top_n)), p_venn, width = 9, height = 8)

    overlap_tbl <- tibble(Feature = sort(unique(unlist(sets)))) %>%
      mutate(
        in_source = Feature %in% names(coefs_by_cohort[[source]]),
        in_target1 = Feature %in% names(coefs_by_cohort[[targets[1]]]),
        in_target2 = Feature %in% names(coefs_by_cohort[[targets[2]]]),
        n_models = in_source + in_target1 + in_target2
      ) %>%
      arrange(desc(n_models), Feature)
    write_csv(overlap_tbl, file.path(archive_dir, sprintf("venn_top%d_3models_features.csv", top_n)))

    consensus <- Reduce(intersect, lapply(coefs_by_cohort, names))
    wide <- tibble(Feature = consensus)
    for (cohort in cohorts) {
      wide[[paste0("coef_", cohort)]] <- abs(coefs_by_cohort[[cohort]][consensus])
    }
    imp_tbl <- wide %>%
      mutate(Mean_AbsCoef = rowMeans(across(starts_with("coef_")), na.rm = TRUE)) %>%
      arrange(desc(Mean_AbsCoef)) %>%
      slice_head(n = ifelse(top_n == 20, 15, 20)) %>%
      pivot_longer(starts_with("coef_"), names_to = "Model", values_to = "AbsCoef") %>%
      mutate(Model = sub("^coef_", "", Model), Model = pretty_name(Model))
    write_csv(imp_tbl, file.path(archive_dir, sprintf("feature_importance_top%d_3models.csv", top_n)))

    p_imp <- ggplot(imp_tbl, aes(reorder(Feature, Mean_AbsCoef), AbsCoef, fill = Model)) +
      geom_col(position = "dodge", width = 0.72) +
      coord_flip() +
      labs(
        title = sprintf("Consensus Feature Coefficients: v%d Fixed Models", v),
        subtitle = "Features present in all three fixed cohort models",
        x = "Feature",
        y = "|Coefficient|",
        fill = "Model"
      ) +
      theme_archive() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    save_pdf(file.path(archive_dir, sprintf("feature_importance_top%d_3models.pdf", top_n)), p_imp, width = 11, height = height)

    file.copy(file.path(archive_dir, sprintf("venn_top%d_3models.pdf", top_n)),
              file.path(summary_dir, sprintf("venn_top%d.pdf", top_n)), overwrite = TRUE)
    file.copy(file.path(archive_dir, sprintf("feature_importance_top%d_3models.pdf", top_n)),
              file.path(summary_dir, sprintf("feature_importance_top%d.pdf", top_n)), overwrite = TRUE)
  }

  make_venn_importance(20, 20, 7.5)
  make_venn_importance(50, 50, 8.5)
  file.copy(summary_file, file.path(summary_dir, "models_on_other2_summary.csv"), overwrite = TRUE)
}

walk(sources, make_one_source)

cat("Wrote fork2 source-transfer plots under:", archive_root, "\n")
