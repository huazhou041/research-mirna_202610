#!/usr/bin/env Rscript
# Assemble the self-contained miRNA AUC dashboard: inject the 4 long CSVs (as JSON)
# and Plotly.js into dashboard_template.html.
suppressPackageStartupMessages(library(jsonlite))
setwd("/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260715")

combined_auc <- read.csv("data/combined_auc_long.csv", stringsAsFactors = FALSE)
combined_roc <- read.csv("data/combined_roc_long.csv", stringsAsFactors = FALSE)
transfer_auc <- read.csv("data/transfer_auc_long.csv", stringsAsFactors = FALSE)
transfer_roc <- read.csv("data/transfer_roc_long.csv", stringsAsFactors = FALSE)
transfer_feat <- read.csv("data/transfer_features_long.csv", stringsAsFactors = FALSE)
own_roc <- if (file.exists("data/own_roc_long.csv")) read.csv("data/own_roc_long.csv", stringsAsFactors = FALSE) else
  data.frame(fs_method=character(), cohort=character(), varcut=integer(), auc=numeric(),
             n_pred=integer(), fpr=character(), tpr=character(), stringsAsFactors = FALSE)

# transfer summary tables (the all_sources_..._summary format), tagged by fs_method
T14 <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260626_full/mirna_20260714/results"
sum_mad <- read.csv(file.path(T14, "mad/transfer/all_sources_fixed_models_on_other2_summary.csv"), stringsAsFactors = FALSE)
sum_cor <- read.csv(file.path(T14, "correlation/transfer/all_sources_corr_transfer_summary.csv"), stringsAsFactors = FALSE)
sum_mad$fs_method <- "mad"; sum_cor$fs_method <- "cor"
cc <- intersect(names(sum_mad), names(sum_cor))
transfer_summary <- rbind(sum_mad[, cc], sum_cor[, cc])
for (col in grep("AUC|CV", names(transfer_summary), value = TRUE)) transfer_summary[[col]] <- round(transfer_summary[[col]], 4)

# keep ROC payload lean: order columns, drop NA-curve rows
combined_roc <- combined_roc[!is.na(combined_roc$fpr), ]
transfer_roc <- transfer_roc[!is.na(transfer_roc$fpr), ]
own_roc <- own_roc[!is.na(own_roc$fpr), ]

j <- function(x) toJSON(x, dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 6)
meta <- list(generated = as.character(Sys.Date()),
             n_combined = nrow(combined_auc), n_transfer = nrow(transfer_auc))

tpl <- paste(readLines("codes/dashboard_template.html", warn = FALSE), collapse = "\n")
plotly <- paste(readLines("codes/plotly-2.32.0.min.js", warn = FALSE), collapse = "\n")

tpl <- sub("__PLOTLY_JS__", plotly, tpl, fixed = TRUE)
tpl <- sub("__COMBINED_AUC__", j(combined_auc), tpl, fixed = TRUE)
tpl <- sub("__COMBINED_ROC__", j(combined_roc), tpl, fixed = TRUE)
tpl <- sub("__TRANSFER_AUC__", j(transfer_auc), tpl, fixed = TRUE)
tpl <- sub("__TRANSFER_ROC__", j(transfer_roc), tpl, fixed = TRUE)
tpl <- sub("__TRANSFER_SUMMARY__", j(transfer_summary), tpl, fixed = TRUE)
tpl <- sub("__TRANSFER_FEAT__", j(transfer_feat), tpl, fixed = TRUE)
tpl <- sub("__OWN_ROC__", j(own_roc), tpl, fixed = TRUE)
tpl <- sub("__META__", toJSON(meta, auto_unbox = TRUE), tpl, fixed = TRUE)
tpl <- gsub("__DATE__", as.character(Sys.Date()), tpl, fixed = TRUE)

out <- "htmls/mirna_auc_dashboard_20260715.html"
dir.create("htmls", showWarnings = FALSE)
writeLines(tpl, out)
cat("Wrote", out, "(", round(file.info(out)$size/1e6, 1), "MB )\n")
cat("combined:", nrow(combined_auc), "auc /", nrow(combined_roc), "roc | transfer:",
    nrow(transfer_auc), "auc /", nrow(transfer_roc), "roc\n")
