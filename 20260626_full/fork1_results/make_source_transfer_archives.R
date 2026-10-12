suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggvenn)
  library(glmnet)
  library(pROC)
})

root <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605"
base_dir <- file.path(root, "20260626_full/fork1_results/source_transfer_archives")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

varlist <- c(5, 10, 20, 50, 100, 200, 500, 1000)

datasets <- list(
  EV_RNA = readRDS(file.path(root, "prev/data/model_data/EV_RNA_miRNA.RDS")),
  nyu_166 = readRDS(file.path(root, "data/model_data/nyu_166_miRNA.RDS")),
  nyu_82 = readRDS(file.path(root, "data/model_data/nyu_82_miRNA.RDS"))
) %>%
  lapply(function(d) d[!is.na(d$recurrence), , drop = FALSE])

model_path <- function(cohort, v) {
  if (cohort == "EV_RNA") {
    file.path(root, sprintf("20260626_full/models/EV_RNA_miRNA_glm_v%d.rds", v))
  } else {
    file.path(root, sprintf("results/raw_pred/%s_miRNA_glm_v%d.rds", cohort, v))
  }
}

cv_path <- function(cohort, v) {
  if (cohort == "EV_RNA") {
    NA_character_
  } else {
    file.path(root, sprintf("results/raw_pred/%s_miRNA_glm_v%d_cvSummary.csv", cohort, v))
  }
}

pretty_name <- function(cohort) {
  recode(cohort, EV_RNA = "EV_RNA", nyu_166 = "NYU 166", nyu_82 = "NYU 82")
}

extract_coefs <- function(model) {
  cm <- coef(model, s = "lambda.min")
  features <- rownames(cm)[-1]
  coefs <- as.numeric(cm[-1, 1])
  names(coefs) <- features
  coefs
}

predict_auc <- function(model, data, features) {
  missing_features <- setdiff(features, colnames(data))
  if (length(missing_features) > 0) {
    stop("Missing features in target data: ", paste(missing_features, collapse = ", "))
  }
  pred <- as.numeric(predict(model, newx = as.matrix(data[, features, drop = FALSE]), s = "lambda.min", type = "response"))
  auc <- as.numeric(pROC::roc(data$recurrence, pred, quiet = TRUE)$auc)
  list(pred = pred, auc = auc)
}

source_cv_auc <- function(cohort, v, model, source_data, features) {
  path <- cv_path(cohort, v)
  if (!is.na(path) && file.exists(path)) {
    cv <- read_csv(path, show_col_types = FALSE)
    return(list(auc = mean(cv$AUC, na.rm = TRUE), source = "mean_cv_auc"))
  }
  apparent <- predict_auc(model, source_data, features)$auc
  list(auc = apparent, source = "apparent_auc")
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

make_one_source <- function(source, targets) {
  out_dir <- file.path(base_dir, paste0(source, "_source"))
  archive_dir <- file.path(out_dir, "Archive")
  summary_dir <- file.path(out_dir, "summary")
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  all_info <- list()

  for (v in varlist) {
    model <- readRDS(model_path(source, v))
    coefs <- extract_coefs(model)
    features <- names(coefs)
    source_auc <- source_cv_auc(source, v, model, datasets[[source]], features)

    row <- tibble(
      Version = v,
      Source = source,
      N_Features = length(features),
      N_Nonzero = sum(coefs != 0),
      Source_AUC = source_auc$auc,
      Source_AUC_Type = source_auc$source
    )

    for (target in targets) {
      res <- predict_auc(model, datasets[[target]], features)
      row[[paste0("AUC_", target)]] <- res$auc
    }

    all_info[[as.character(v)]] <- list(
      summary = row,
      coefs = coefs,
      features = features
    )
  }

  summary_df <- bind_rows(lapply(all_info, `[[`, "summary")) %>%
    arrange(Version)

  write_csv(summary_df, file.path(archive_dir, sprintf("%s_models_on_other2_summary.csv", source)))
  write_csv(summary_df, file.path(summary_dir, "models_on_other2_summary.csv"))
  write_csv(
    summary_df %>% slice_max(Source_AUC, n = 1, with_ties = FALSE),
    file.path(archive_dir, sprintf("%s_source_best_cv_model.csv", source))
  )

  long_auc <- summary_df %>%
    select(Version, Source_AUC, starts_with("AUC_")) %>%
    pivot_longer(-Version, names_to = "Cohort", values_to = "AUC") %>%
    mutate(
      Cohort = case_when(
        Cohort == "Source_AUC" ~ paste0(pretty_name(source), " source CV"),
        TRUE ~ paste0("Applied to ", pretty_name(sub("^AUC_", "", Cohort)))
      ),
      Version = factor(Version, levels = varlist)
    )

  p_auc <- ggplot(long_auc, aes(Version, AUC, color = Cohort, group = Cohort)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    labs(x = "Feature cutoff", y = "AUC", color = NULL) +
    theme_archive()
  save_pdf(file.path(archive_dir, "auc_by_varcut.pdf"), p_auc)

  feature_tables <- map_dfr(varlist, function(v) {
    coefs <- all_info[[as.character(v)]]$coefs
    tibble(
      Version = v,
      Feature = names(coefs),
      Coefficient = as.numeric(coefs),
      AbsCoef = abs(as.numeric(coefs)),
      Nonzero = as.numeric(coefs) != 0
    ) %>%
      arrange(desc(AbsCoef))
  })

  write_csv(feature_tables, file.path(archive_dir, sprintf("%s_all_model_features.csv", source)))
  for (v in varlist) {
    write_csv(
      feature_tables %>% filter(Version == v),
      file.path(archive_dir, sprintf("%s_v%d_features.csv", source, v))
    )
  }

  make_venn_and_importance <- function(v, top_n, height) {
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
        title = sprintf("Feature Overlap: v%d GLM Models", v),
        subtitle = paste(pretty_name(source), "source model and the other two cohorts")
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
    if (length(consensus) == 0) {
      imp_tbl <- tibble(
        Feature = character(),
        Model = character(),
        AbsCoef = numeric(),
        Mean_AbsCoef = numeric()
      )
    } else {
      wide <- tibble(Feature = consensus)
      for (cohort in cohorts) {
        wide[[paste0("coef_", cohort)]] <- abs(coefs_by_cohort[[cohort]][consensus])
      }
      wide <- wide %>%
        mutate(Mean_AbsCoef = rowMeans(across(starts_with("coef_")), na.rm = TRUE)) %>%
        arrange(desc(Mean_AbsCoef)) %>%
        slice_head(n = ifelse(top_n == 20, 15, 20))
      imp_tbl <- wide %>%
        pivot_longer(starts_with("coef_"), names_to = "Model", values_to = "AbsCoef") %>%
        mutate(Model = sub("^coef_", "", Model), Model = pretty_name(Model))
    }

    write_csv(imp_tbl, file.path(archive_dir, sprintf("feature_importance_top%d_3models.csv", top_n)))

    if (nrow(imp_tbl) == 0) {
      p_imp <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No features overlap across all 3 models") +
        theme_void()
    } else {
      p_imp <- ggplot(imp_tbl, aes(reorder(Feature, Mean_AbsCoef), AbsCoef, fill = Model)) +
        geom_col(position = "dodge", width = 0.72) +
        coord_flip() +
        labs(
          title = sprintf("Consensus Feature Coefficients: v%d Models", v),
          subtitle = "Features present in all three cohort models",
          x = "Feature",
          y = "|Coefficient|",
          fill = "Model"
        ) +
        theme_archive() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    }

    save_pdf(file.path(archive_dir, sprintf("feature_importance_top%d_3models.pdf", top_n)), p_imp, width = 11, height = height)
    file.copy(
      file.path(archive_dir, sprintf("feature_importance_top%d_3models.pdf", top_n)),
      file.path(summary_dir, sprintf("feature_importance_top%d.pdf", top_n)),
      overwrite = TRUE
    )
    file.copy(
      file.path(archive_dir, sprintf("venn_top%d_3models.pdf", top_n)),
      file.path(summary_dir, sprintf("venn_top%d.pdf", top_n)),
      overwrite = TRUE
    )
  }

  make_venn_and_importance(20, 20, 7.5)
  make_venn_and_importance(50, 50, 8.5)

  invisible(summary_df)
}

summary_82 <- make_one_source("nyu_82", c("EV_RNA", "nyu_166"))
summary_166 <- make_one_source("nyu_166", c("EV_RNA", "nyu_82"))

cat("Wrote source-transfer archives under:", base_dir, "\n")
cat("NYU82 source rows:", nrow(summary_82), "\n")
cat("NYU166 source rows:", nrow(summary_166), "\n")
