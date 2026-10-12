#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(pROC)
  library(jsonlite)
  library(caret)
})

setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2")
dir.create("models", recursive=TRUE, showWarnings=FALSE)
dir.create("source_transfer_archives", recursive=TRUE, showWarnings=FALSE)

varcuts <- c(5, 10, 20, 50, 100, 200, 500, 1000)
sources <- c("EV_RNA", "nyu_166", "nyu_82")
hyperparams <- jsonlite::fromJSON("codes/hyperparams.json", simplifyVector=TRUE)

source_cv_files <- function(source, v) {
  if (source == "EV_RNA") {
    file.path("/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/res", sprintf("EV_RNA_miRNA_glm_v%d_cvSummary.csv", v))
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

data_paths <- c(EV_RNA="rdata/EV_RNA_miRNA.RDS", nyu_166="rdata/nyu_166_miRNA.RDS", nyu_82="rdata/nyu_82_miRNA.RDS")
data_list <- lapply(data_paths, readRDS)
data_list <- lapply(data_list, function(x) {
  x <- x[!is.na(x$recurrence), , drop=FALSE]
  x$recurrence <- as.numeric(as.character(x$recurrence))
  x
})

select_whole_features <- function(dat, varcut) {
  tissues <- grep("^Tissue_", colnames(dat), value=TRUE)
  genes <- setdiff(colnames(dat), c("recurrence", tissues))
  mad_vals <- apply(dat[, genes, drop=FALSE], 2, mad, na.rm=TRUE)
  c(tissues, genes[order(mad_vals, decreasing=TRUE)[seq_len(min(varcut, length(genes)))]])
}

fit_fixed_glm <- function(dat, features, source, varcut) {
  x <- as.matrix(dat[, features, drop=FALSE])
  y <- as.numeric(as.character(dat$recurrence))
  set.seed(123)  # fresh seed immediately before each fold draw (per-varcut reseed)
  foldid <- caret::createFolds(y, k=10, list=FALSE)  # 10-fold lambda selection (consistent w/ 10-fold pipeline)
  models <- lapply(hyperparams, function(alpha) {
    tryCatch(
      cv.glmnet(x, y, family="binomial", alpha=alpha, foldid=foldid),
      error=function(e) cv.glmnet(x, y, family="binomial", alpha=alpha, foldid=foldid, lambda=exp(seq(log(0.001), log(5), length.out=100)))
    )
  })
  cvms <- sapply(models, function(m) m$cvm[which.min(abs(m$lambda - m$lambda.min))])
  best_i <- which.min(cvms)
  model <- models[[best_i]]
  attr(model, "source") <- source
  attr(model, "varcut") <- varcut
  attr(model, "features") <- features
  attr(model, "alpha") <- hyperparams[[best_i]]
  attr(model, "fixed_final_feature_set") <- "whole_source_dataset_top_mad"
  model
}

extract_coefs <- function(model) {
  cm <- coef(model, s="lambda.min")
  coefs <- as.numeric(cm[-1, 1])
  names(coefs) <- rownames(cm)[-1]
  coefs
}

predict_auc <- function(model, dat, features) {
  missing <- setdiff(features, colnames(dat))
  if (length(missing) > 0) stop("Missing features: ", paste(missing, collapse=", "))
  pred <- as.numeric(predict(model, newx=as.matrix(dat[, features, drop=FALSE]), s="lambda.min", type="response"))
  auc <- as.numeric(pROC::roc(dat$recurrence, pred, quiet=TRUE, direction = "<")$auc)  # direction="<": avoid autoflip inflating AUC on anti-correlated folds
  list(pred=pred, auc=auc)
}

source_cv_auc <- function(source, v) {
  f <- source_cv_files(source, v)
  if (!file.exists(f)) return(NA_real_)
  mean(read.csv(f)$AUC, na.rm=TRUE)
}

feature_audit_rows <- list()
summary_by_source <- list()

for (source in sources) {
  archive_dir <- file.path("source_transfer_archives", paste0(source, "_source"), "Archive")
  summary_dir <- file.path("source_transfer_archives", paste0(source, "_source"), "summary")
  dir.create(archive_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(summary_dir, recursive=TRUE, showWarnings=FALSE)
  summary_rows <- list()
  feature_rows <- list()

  for (v in varcuts) {
    features <- select_whole_features(data_list[[source]], v)
    model <- fit_fixed_glm(data_list[[source]], features, source, v)
    saveRDS(model, file.path("models", sprintf("%s_miRNA_glm_v%d_fixed_whole_features.rds", source, v)))
    coefs <- extract_coefs(model)
    row <- data.frame(Version=v, Source=source, N_Features=length(coefs), N_Nonzero=sum(coefs != 0), Source_CV_AUC=source_cv_auc(source, v), stringsAsFactors=FALSE)
    for (target in setdiff(sources, source)) {
      res <- predict_auc(model, data_list[[target]], names(coefs))
      row[[paste0("AUC_", target)]] <- res$auc
    }
    summary_rows[[as.character(v)]] <- row
    feature_rows[[as.character(v)]] <- data.frame(Version=v, Feature=names(coefs), Coefficient=as.numeric(coefs), AbsCoef=abs(as.numeric(coefs)), Nonzero=as.numeric(coefs) != 0, stringsAsFactors=FALSE)

    old_file <- old_model_files(source, v)
    if (file.exists(old_file)) {
      old_model <- readRDS(old_file)
      old_features <- rownames(coef(old_model, s="lambda.min"))[-1]
      feature_audit_rows[[length(feature_audit_rows)+1]] <- data.frame(Source=source, Version=v, Old_Matches_Whole=setequal(old_features, features), Old_Missing_Whole_Features=length(setdiff(features, old_features)), Old_Extra_vs_Whole=length(setdiff(old_features, features)), Fixed_Matches_Whole=setequal(names(coefs), features), stringsAsFactors=FALSE)
    }
  }
  summary_df <- bind_rows(summary_rows) %>% arrange(Version)
  summary_by_source[[source]] <- summary_df
  write.csv(summary_df, file.path(archive_dir, sprintf("%s_models_on_other2_summary.csv", source)), row.names=FALSE, quote=FALSE)
  write.csv(summary_df, file.path(summary_dir, "models_on_other2_summary.csv"), row.names=FALSE, quote=FALSE)
  write.csv(bind_rows(feature_rows), file.path(archive_dir, sprintf("%s_all_model_features.csv", source)), row.names=FALSE, quote=FALSE)
}

feature_audit <- bind_rows(feature_audit_rows)
write.csv(feature_audit, "fixed_model_feature_audit.csv", row.names=FALSE, quote=FALSE)
all_summary <- bind_rows(summary_by_source)
write.csv(all_summary, "all_sources_fixed_models_on_other2_summary.csv", row.names=FALSE, quote=FALSE)
print(feature_audit)
print(all_summary)
cat("fork2 train/apply fixed complete\n")
