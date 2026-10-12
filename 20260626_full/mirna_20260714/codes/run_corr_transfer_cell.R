#!/usr/bin/env Rscript
# One (source, varcut) cell of the correlation-FS 3x3 transfer analysis.
# Self-contained: computes the transfer AUCs (whole-source corr-FS GLM applied to
# the other two cohorts) AND the honest within-cohort Source_CV_AUC (100x10
# bootstrap, correlation feature selection re-done inside each training fold).
# Writes one row to cells/<source>_v<varcut>.csv. Run as a SLURM array; a finalize
# step concatenates the 24 cells. pROC direction="<" throughout.
suppressPackageStartupMessages({
  library(glmnet); library(pROC); library(caret); library(parallel); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
source_i <- args[1]; varcut <- as.integer(args[2])
root <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/corr_transfer"
setwd(root)

sources <- c("EV_RNA", "nyu_166", "nyu_82")
alphas  <- fromJSON("codes/hyperparams.json", simplifyVector = TRUE)
NCORES  <- as.integer(Sys.getenv("SLURM_NTASKS", unset = "4"))
NBOOT   <- 50; KFOLD <- 10   # 50x10 = 500 fold-AUCs -> stable mean; reference column

data_list <- lapply(setNames(file.path("rdata", paste0(sources, "_miRNA.RDS")), sources), function(p) {
  x <- readRDS(p); x <- x[!is.na(x$recurrence), , drop = FALSE]
  x$recurrence <- as.numeric(as.character(x$recurrence)); x
})

select_corr <- function(dat, vc, y = dat$recurrence) {
  tissues <- grep("^Tissue_", colnames(dat), value = TRUE)
  genes <- setdiff(colnames(dat), c("recurrence", tissues))
  cc <- abs(as.numeric(suppressWarnings(cor(as.matrix(dat[, genes, drop = FALSE]), y))))
  cc[is.na(cc)] <- -Inf
  c(tissues, genes[order(cc, decreasing = TRUE)[seq_len(min(vc, length(genes)))]])
}
fit_glm <- function(x, y) {
  set.seed(123); foldid <- caret::createFolds(y, k = 10, list = FALSE)
  ms <- lapply(alphas, function(a)
    tryCatch(cv.glmnet(x, y, family = "binomial", alpha = a, foldid = foldid),
             error = function(e) cv.glmnet(x, y, family = "binomial", alpha = a, foldid = foldid,
                                           lambda = exp(seq(log(0.001), log(5), length.out = 100)))))
  ms[[which.min(sapply(ms, function(m) m$cvm[which.min(abs(m$lambda - m$lambda.min))]))]]
}
auc_dir <- function(y, p) as.numeric(pROC::roc(y, p, quiet = TRUE, direction = "<")$auc)

dat <- data_list[[source_i]]; y <- dat$recurrence
# --- transfer model + apply to other 2 cohorts ---
feats <- select_corr(dat, varcut)
m <- fit_glm(as.matrix(dat[, feats, drop = FALSE]), y)
cm <- coef(m, s = "lambda.min"); coefs <- setNames(as.numeric(cm[-1, 1]), rownames(cm)[-1])
row <- data.frame(Version = varcut, Source = source_i, N_Features = length(coefs),
                  N_Nonzero = sum(coefs != 0), Source_CV_AUC = NA_real_, stringsAsFactors = FALSE)
for (tg in setdiff(sources, source_i)) {
  td <- data_list[[tg]]
  p <- as.numeric(predict(m, as.matrix(td[, names(coefs), drop = FALSE]), s = "lambda.min", type = "response"))
  row[[paste0("AUC_", tg)]] <- auc_dir(td$recurrence, p)
}
# --- Source_CV_AUC: fold-internal correlation FS ---
boot_auc <- unlist(mclapply(seq_len(NBOOT), function(b) {
  set.seed(1000 + b)
  folds <- caret::createFolds(factor(y), k = KFOLD, list = TRUE, returnTrain = FALSE)
  mean(vapply(seq_len(KFOLD), function(i) {
    va <- folds[[i]]; tr <- setdiff(seq_along(y), va)
    ff <- select_corr(dat[tr, , drop = FALSE], varcut, y[tr])
    mm <- fit_glm(as.matrix(dat[tr, ff, drop = FALSE]), y[tr])
    pp <- as.numeric(predict(mm, as.matrix(dat[va, ff, drop = FALSE]), s = "lambda.min", type = "response"))
    auc_dir(y[va], pp)
  }, numeric(1)), na.rm = TRUE)
}, mc.cores = NCORES))
row$Source_CV_AUC <- mean(boot_auc, na.rm = TRUE)

dir.create("cells", showWarnings = FALSE)
write.csv(row[, c("Version","Source","N_Features","N_Nonzero","Source_CV_AUC", grep("^AUC_", names(row), value=TRUE))],
          file.path("cells", sprintf("%s_v%d.csv", source_i, varcut)), row.names = FALSE, quote = FALSE)
cat(sprintf("done %s v%d: Source_CV_AUC=%.4f\n", source_i, varcut, row$Source_CV_AUC))
