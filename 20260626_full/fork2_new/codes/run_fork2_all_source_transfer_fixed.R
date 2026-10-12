#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(pROC)
  library(jsonlite)
  library(caret)
})

setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2")

varcuts <- c(5, 10, 20, 50, 100, 200, 500, 1000)
sources <- c("EV_RNA", "nyu_166", "nyu_82")
hyperparams <- jsonlite::fromJSON("codes/hyperparams.json", simplifyVector=TRUE)

source_cv_files <- function(source, v) {
  if (source == "EV_RNA") {
    file.path("/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/raw_pred", sprintf("EV_RNA_miRNA_glm_v%d_cvSummary.csv", v))
  } else {
    file.path("/gpfs/scratch/zhouh05/aris/segal/val_20260605/res", sprintf("%s_miRNA_glm_v%d_cvSummary.csv", source, v))
  }
}

old_model_files <- function(source, v) {
  if (source == "EV_RNA") {
    file.path("/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/res", sprintf("EV_RNA_miRNA_glm_v%d.rds", v))
  } else {
    file.path("/gpfs/scratch/zhouh05/aris/segal/val_20260605/res", sprintf("%s_miRNA_glm_v%d.rds", source, v))
  }
}

data_paths <- c(
  EV_RNA = "rdata/EV_RNA_miRNA.RDS",
  nyu_166 = "rdata/nyu_166_miRNA.RDS",
  nyu_82 = "rdata/nyu_82_miRNA.RDS"
)

data_list <- lapply(data_paths, readRDS)
data_list <- lapply(data_list, function(x) {
  x <- x[!is.na(x$recurrence), , drop=FALSE]
  x$recurrence <- as.numeric(as.character(x$recurrence))
  x
})

pretty_name <- function(x) dplyr::recode(x, EV_RNA="EV_RNA", nyu_166="NYU 166", nyu_82="NYU 82")

feature_audit <- bind_rows(feature_audit_rows)
write.csv(feature_audit, "fixed_model_feature_audit.csv", row.names=FALSE, quote=FALSE)
all_summary <- bind_rows(summary_by_source)
write.csv(all_summary, "all_sources_fixed_models_on_other2_summary.csv", row.names=FALSE, quote=FALSE)
print(feature_audit)
print(all_summary)
cat("fork2 all-source fixed transfer complete\n")
