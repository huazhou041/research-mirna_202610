#!/usr/bin/env Rscript
# ============================================================================
# Correlation-FS version of the 3x3 cross-cohort transfer analysis, mirroring
# fork2's run_fork2_train_apply_fixed.R (MAD version) but selecting features by
# |Pearson correlation with recurrence| instead of MAD.
#
#   * transfer model: refit a GLM on the WHOLE source cohort using the top-|corr|
#     features (selection on source labels is fine here -- the model is evaluated
#     on OTHER cohorts, which are inherently held out), then apply to the other
#     two cohorts -> AUC_<target> (pROC direction="<").
#   * Source_CV_AUC: honest within-cohort CV (100x10 bootstrap) with correlation
#     feature selection re-done INSIDE each training fold (leakage-proof).
#
# Output mirrors fork2: all_sources_corr_transfer_summary.csv + per-source files.
# ============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(glmnet); library(pROC); library(caret); library(parallel); library(jsonlite)
})

base   <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260626_full"
rdata  <- file.path(base, "fork2_results/rdata")
outdir <- file.path(base, "mirna_20260714/results/correlation/transfer")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

varcuts <- c(5, 10, 20, 50, 100, 200, 500, 1000)
sources <- c("EV_RNA", "nyu_166", "nyu_82")
alphas  <- fromJSON(file.path(base, "fork2_new/codes/hyperparams.json"), simplifyVector = TRUE)  # 21-alpha grid
NCORES  <- max(1, detectCores() - 1)
NBOOT   <- 100; KFOLD <- 10

data_list <- lapply(setNames(file.path(rdata, paste0(sources, "_miRNA.RDS")), sources), function(p) {
  x <- readRDS(p); x <- x[!is.na(x$recurrence), , drop = FALSE]
  x$recurrence <- as.numeric(as.character(x$recurrence)); x
})

# correlation feature selection: top-|corr(gene, recurrence)|, tissues always kept
select_corr <- function(dat, varcut, y = dat$recurrence) {
  tissues <- grep("^Tissue_", colnames(dat), value = TRUE)
  genes <- setdiff(colnames(dat), c("recurrence", tissues))
  cc <- abs(as.numeric(suppressWarnings(cor(as.matrix(dat[, genes, drop = FALSE]), y))))
  cc[is.na(cc)] <- -Inf
  c(tissues, genes[order(cc, decreasing = TRUE)[seq_len(min(varcut, length(genes)))]])
}

fit_glm <- function(x, y) {
  set.seed(123); foldid <- caret::createFolds(y, k = 10, list = FALSE)
  models <- lapply(alphas, function(a)
    tryCatch(cv.glmnet(x, y, family = "binomial", alpha = a, foldid = foldid),
             error = function(e) cv.glmnet(x, y, family = "binomial", alpha = a, foldid = foldid,
                                           lambda = exp(seq(log(0.001), log(5), length.out = 100)))))
  cvms <- sapply(models, function(m) m$cvm[which.min(abs(m$lambda - m$lambda.min))])
  models[[which.min(cvms)]]
}

auc_dir <- function(y, p) as.numeric(pROC::roc(y, p, quiet = TRUE, direction = "<")$auc)

# honest within-cohort CV with fold-internal correlation FS (leakage-proof)
source_cv_corr <- function(dat, varcut) {
  y <- dat$recurrence
  boot_auc <- unlist(mclapply(seq_len(NBOOT), function(b) {
    set.seed(1000 + b)
    folds <- caret::createFolds(factor(y), k = KFOLD, list = TRUE, returnTrain = FALSE)
    mean(vapply(seq_len(KFOLD), function(i) {
      va <- folds[[i]]; tr <- setdiff(seq_along(y), va)
      feats <- select_corr(dat[tr, , drop = FALSE], varcut, y[tr])   # FS on TRAIN fold only
      m <- fit_glm(as.matrix(dat[tr, feats, drop = FALSE]), y[tr])
      p <- as.numeric(predict(m, as.matrix(dat[va, feats, drop = FALSE]), s = "lambda.min", type = "response"))
      auc_dir(y[va], p)
    }, numeric(1)), na.rm = TRUE)
  }, mc.cores = NCORES))
  mean(boot_auc, na.rm = TRUE)
}

# ---- 1. transfer table (fast): whole-source corr-FS GLM applied to other cohorts ----
message(sprintf("[%s] transfer AUCs ...", format(Sys.time(), "%H:%M:%S")))
rows <- list(); feat_rows <- list()
for (source in sources) {
  for (v in varcuts) {
    feats <- select_corr(data_list[[source]], v)
    x <- as.matrix(data_list[[source]][, feats, drop = FALSE]); y <- data_list[[source]]$recurrence
    m <- fit_glm(x, y)
    cm <- coef(m, s = "lambda.min"); coefs <- setNames(as.numeric(cm[-1, 1]), rownames(cm)[-1])
    row <- data.frame(Version = v, Source = source, N_Features = length(coefs),
                      N_Nonzero = sum(coefs != 0), Source_CV_AUC = NA_real_, stringsAsFactors = FALSE)
    for (target in setdiff(sources, source)) {
      td <- data_list[[target]]
      p <- as.numeric(predict(m, as.matrix(td[, names(coefs), drop = FALSE]), s = "lambda.min", type = "response"))
      row[[paste0("AUC_", target)]] <- auc_dir(td$recurrence, p)
    }
    rows[[paste(source, v)]] <- row
    feat_rows[[paste(source, v)]] <- data.frame(Version = v, Source = source, Feature = names(coefs),
                                                Coefficient = as.numeric(coefs), stringsAsFactors = FALSE)
  }
}
transfer <- bind_rows(rows)
write.csv(transfer, file.path(outdir, "all_sources_corr_transfer_summary_NO_SOURCECV.csv"), row.names = FALSE, quote = FALSE)
write.csv(bind_rows(feat_rows), file.path(outdir, "corr_transfer_model_features.csv"), row.names = FALSE, quote = FALSE)
message(sprintf("[%s] transfer table written (Source_CV_AUC pending)", format(Sys.time(), "%H:%M:%S")))

# ---- 2. Source_CV_AUC (correlation FS, fold-internal) ----
for (source in sources) for (v in varcuts) {
  scv <- source_cv_corr(data_list[[source]], v)
  transfer$Source_CV_AUC[transfer$Source == source & transfer$Version == v] <- scv
  message(sprintf("[%s] Source_CV_AUC %s v%d = %.4f", format(Sys.time(), "%H:%M:%S"), source, v, scv))
}

# reorder columns like fork2 (Source_CV_AUC before the AUC_ targets)
auc_cols <- grep("^AUC_", names(transfer), value = TRUE)
transfer <- transfer[, c("Version", "Source", "N_Features", "N_Nonzero", "Source_CV_AUC", auc_cols)]
transfer <- transfer %>% arrange(match(Source, sources), Version)
write.csv(transfer, file.path(outdir, "all_sources_corr_transfer_summary.csv"), row.names = FALSE, quote = FALSE)
for (source in sources) {
  write.csv(transfer %>% filter(Source == source),
            file.path(outdir, sprintf("%s_corr_models_on_other2_summary.csv", source)), row.names = FALSE, quote = FALSE)
}
message(sprintf("[%s] DONE -> all_sources_corr_transfer_summary.csv", format(Sys.time(), "%H:%M:%S")))
print(transfer, digits = 3)
