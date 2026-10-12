#!/usr/bin/env Rscript
# ============================================================================
# fork3: CV-HONEST version of the combined 3-cohort batch-handling comparison.
#
# The fork2 pipeline built the "combat"/"both" input matrices by running
# limma::removeBatchEffect ONCE on the full pooled dataset (all samples, and
# with design = ~recurrence using EVERY sample's label) BEFORE any CV split.
# Every bootstrap fold then read that pre-corrected matrix, so each held-out
# validation sample's batch correction was informed by (a) other validation
# samples and (b) its own outcome label -> optimistic batch-effect + outcome
# leakage into the reported CV AUC.
#
# fork3 fixes this: for the combat/both variants the batch correction is fit
# INSIDE each training fold, on the training samples only, using only the
# training labels to protect biological signal, and then applied to the
# held-out fold using ONLY that fold's cohort membership (never its labels).
# plain/cov carry no feature correction and are re-run here unchanged for a
# clean, same-seed 4-way comparison.
#
# Usage:
#   Rscript tune_combined_models_cvcombat.R <raw_matrix.RDS> <metadata.csv> \
#           <out_prefix> --variant {plain|cov|combat|both} \
#           --learner {GLM|RF|XGB} --varcut N --nboot 100 --kfold 10
# ============================================================================
suppressPackageStartupMessages({
  library(argparser)
  library(tidyverse)
  library(parallel)
  library(caret)
  library(pROC)
})

p <- arg_parser("CV-honest fold-internal ComBat combined-cohort trainer")
p <- add_argument(p, "matrix", help = "RAW pooled RDS matrix (recurrence + genes, no correction)")
p <- add_argument(p, "metadata", help = "CSV with sample_id, cohort, recurrence (row-aligned to matrix)")
p <- add_argument(p, "prefix", help = "Output prefix (writes res/, raw_pred/)")
p <- add_argument(p, "--variant", default = "combat", help = "plain | cov | combat | both")
p <- add_argument(p, "--learner", default = "GLM", help = "GLM | RF | XGB")
p <- add_argument(p, "--varcut", default = 100, help = "Number of selected gene features")
p <- add_argument(p, "--fs_method", default = "mad", help = "Feature selection: mad (unsupervised variance) | corr (|Pearson corr with recurrence|, supervised)")
p <- add_argument(p, "--nboot", default = 100, help = "Bootstrap repeats")
p <- add_argument(p, "--kfold", default = 10, help = "CV folds")
p <- add_argument(p, "--cores", default = "Slurm", help = "future::availableCores method")
argv <- parse_args(p)

drug <- "recurrence"
variant <- tolower(argv$variant)
learner <- toupper(argv$learner)
varcut <- as.integer(argv$varcut)
fs_method <- tolower(argv$fs_method)
kboot <- as.integer(argv$nboot)
kfold <- as.integer(argv$kfold)
stopifnot(variant %in% c("plain", "cov", "combat", "both"),
          learner %in% c("GLM", "RF", "XGB"),
          fs_method %in% c("mad", "corr"))

do_combat <- variant %in% c("combat", "both")   # fold-internal batch correction?
do_cov    <- variant %in% c("cov", "both")      # cohort indicator covariates?
batch_levels <- c("EV_RNA", "nyu_166", "nyu_82")
Csum <- contr.sum(length(batch_levels))         # limma removeBatchEffect uses sum-to-zero contrasts
rownames(Csum) <- batch_levels

mat <- readRDS(argv$matrix)
meta <- read.csv(argv$metadata, stringsAsFactors = FALSE)
stopifnot(all(meta$sample_id == rownames(mat)))
mat[[drug]] <- as.numeric(as.character(mat[[drug]]))
keep <- !is.na(mat[[drug]])
mat <- mat[keep, , drop = FALSE]
meta <- meta[keep, , drop = FALSE]
cohort_all <- factor(meta$cohort, levels = batch_levels)

genes <- setdiff(colnames(mat), drug)            # raw gene features only (no batch cols in the raw matrix)
G_all <- as.matrix(mat[, genes, drop = FALSE])   # samples x genes (RAW)
y_all <- mat[[drug]]

Ncores <- future::availableCores(method = argv$cores)
options(mc.cores = Ncores)
set.seed(123)

glm_alphas <- tryCatch(jsonlite::fromJSON("codes/hyperparams.json", simplifyVector = TRUE),
                       error = function(e) c(0, 0.25, 0.5, 0.75, 1))

# ---- fold-internal ComBat: fit batch coefficients on TRAIN, apply to both ----
# Replicates limma::removeBatchEffect(x, batch=cohort, design=~recurrence):
#   * sum-to-zero batch contrasts
#   * recurrence protected (in the fit, never subtracted)
#   * only the batch-coefficient contribution is removed
# Fit uses train samples/labels only; the held-out fold is corrected using its
# cohort membership alone (no validation labels enter the correction).
combat_train_apply <- function(Gtr, ytr, coh_tr, Gva, coh_va) {
  Btr <- Csum[as.character(coh_tr), , drop = FALSE]           # ntrain x (L-1)
  Bva <- Csum[as.character(coh_va), , drop = FALSE]           # nvalid x (L-1)
  D <- cbind(`(Intercept)` = 1, recurrence = ytr)             # protected design (train labels)
  M <- cbind(D, Btr)                                          # ntrain x p
  beta <- t(Gtr) %*% M %*% solve(crossprod(M))                # genes x p  (OLS, = lmFit)
  bcols <- (ncol(D) + 1):(ncol(D) + ncol(Btr))
  beta_b <- beta[, bcols, drop = FALSE]                       # genes x (L-1)
  beta_b[is.na(beta_b)] <- 0
  list(train = Gtr - Btr %*% t(beta_b),                       # subtract batch effect only
       valid = Gva - Bva %*% t(beta_b))
}

# Feature selection on TRAIN samples only (leakage-proof). Two criteria:
#   mad  = unsupervised, rank genes by median absolute deviation (variance)
#   corr = supervised, rank genes by |Pearson corr with the TRAIN recurrence
#          label| (point-biserial). Uses training labels only -> must stay
#          inside the fold; never sees validation samples or labels.
FS_top <- function(Gtr_features, top_n, ytr) {
  score <- if (fs_method == "corr") {
    abs(as.numeric(suppressWarnings(cor(Gtr_features, ytr))))   # |point-biserial corr|
  } else {
    apply(Gtr_features, 2, mad, na.rm = TRUE)
  }
  score[is.na(score)] <- -Inf                                   # constant-in-fold genes rank last
  colnames(Gtr_features)[order(score, decreasing = TRUE)[seq_len(min(top_n, ncol(Gtr_features)))]]
}

build_fold_data <- function(tr_idx, va_idx) {
  if (do_combat) {
    cc <- combat_train_apply(G_all[tr_idx, , drop = FALSE], y_all[tr_idx], cohort_all[tr_idx],
                             G_all[va_idx, , drop = FALSE], cohort_all[va_idx])
    Gtr <- cc$train; Gva <- cc$valid
  } else {
    Gtr <- G_all[tr_idx, , drop = FALSE]; Gva <- G_all[va_idx, , drop = FALSE]
  }
  top <- FS_top(Gtr, varcut, y_all[tr_idx])                   # selection on TRAIN samples/labels only
  Xtr <- Gtr[, top, drop = FALSE]; Xva <- Gva[, top, drop = FALSE]
  if (do_cov) {                                               # cohort indicators (leak-free)
    ind <- function(coh) cbind(batch_nyu_166 = as.integer(coh == "nyu_166"),
                               batch_nyu_82  = as.integer(coh == "nyu_82"))
    Xtr <- cbind(ind(cohort_all[tr_idx]), Xtr)
    Xva <- cbind(ind(cohort_all[va_idx]), Xva)
  }
  list(Xtr = Xtr, Xva = Xva,
       ytr = y_all[tr_idx], yva = y_all[va_idx],
       batch_cols = if (do_cov) c("batch_nyu_166", "batch_nyu_82") else character(0))
}

train_predict <- function(dat) {
  features <- colnames(dat$Xtr)
  if (learner == "GLM") {
    library(glmnet)
    pf <- ifelse(features %in% dat$batch_cols, 0, 1)          # batch covariates unpenalized
    models <- lapply(glm_alphas, function(a)
      tryCatch(cv.glmnet(dat$Xtr, dat$ytr, family = "binomial", alpha = a, nfolds = 5, penalty.factor = pf),
               error = function(e) cv.glmnet(dat$Xtr, dat$ytr, family = "binomial", alpha = a, penalty.factor = pf,
                                             lambda = exp(seq(log(0.001), log(5), length.out = 100)))))
    cvms <- sapply(models, function(m) m$cvm[m$lambda == m$lambda.min])
    model <- models[[which.min(cvms)]]
    pred <- as.vector(predict(model, dat$Xva, s = "lambda.min", type = "response"))
  } else if (learner == "RF") {
    library(randomForest)
    model <- randomForest(x = as.data.frame(dat$Xtr), y = as.factor(dat$ytr),
                          ntree = 500, mtry = max(1, floor(sqrt(length(features)))))
    pred <- predict(model, as.data.frame(dat$Xva), type = "prob")[, 2]
  } else {
    library(xgboost)
    dtrain <- xgb.DMatrix(dat$Xtr, label = dat$ytr)
    params <- list(booster = "gbtree", objective = "binary:logistic", eval_metric = "auc",
                   max_depth = 6, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8)
    model <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
    pred <- predict(model, dat$Xva)
  }
  list(model = model, pred = pred,
       auc = as.numeric(roc(dat$yva, pred, quiet = TRUE, direction = "<")$auc))
}

cvtotal <- vector("list", kboot)
cvpred  <- vector("list", kboot)
message(sprintf("[%s] START variant=%s learner=%s varcut=%d fs=%s n=%d genes=%d | combat=%s cov=%s | nboot=%d kfold=%d",
                format(Sys.time(), "%H:%M:%S"), variant, learner, varcut, fs_method,
                nrow(mat), length(genes), do_combat, do_cov, kboot, kfold))
.t0 <- Sys.time()
for (tt in seq_len(kboot)) {
  .tt0 <- Sys.time()
  flds <- createFolds(factor(y_all), k = kfold, list = TRUE, returnTrain = FALSE)
  fold_id <- integer(length(y_all))
  for (i in seq_len(kfold)) fold_id[flds[[i]]] <- i

  fold_rows <- vector("list", kfold)
  pred_rows <- vector("list", kfold)
  for (i in seq_len(kfold)) {
    tr_idx <- which(fold_id != i); va_idx <- which(fold_id == i)
    fit <- train_predict(build_fold_data(tr_idx, va_idx))
    fold_rows[[i]] <- data.frame(AUC = fit$auc, fold_id = i, bootid = tt)
    pred_rows[[i]] <- data.frame(subj = rownames(mat)[va_idx], foldid = i, bootid = tt, prob = fit$pred)
  }
  cvtotal[[tt]] <- bind_rows(fold_rows)
  cvpred[[tt]]  <- bind_rows(pred_rows)
  message(sprintf("[%s] boot %3d/%d | %.1fs/boot | mean fold AUC=%.3f | total %.1f min",
                  format(Sys.time(), "%H:%M:%S"), tt, kboot,
                  as.numeric(difftime(Sys.time(), .tt0, units = "secs")),
                  mean(cvtotal[[tt]]$AUC, na.rm = TRUE),
                  as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
}

cvtotal <- bind_rows(cvtotal)
cvpred  <- bind_rows(cvpred)
dir.create("raw_pred", recursive = TRUE, showWarnings = FALSE)
dir.create("res", recursive = TRUE, showWarnings = FALSE)
base <- basename(argv$prefix)
write.csv(cvtotal, file.path("raw_pred", sprintf("%s_%s_v%d_cvSummary.csv", base, learner, varcut)), row.names = FALSE, quote = FALSE)
write.csv(cvpred,  file.path("raw_pred", sprintf("%s_%s_v%d_pred_prob.csv",  base, learner, varcut)), row.names = FALSE, quote = FALSE)

# Final deployed model: correct on ALL data (no held-out set here), train on all.
final_dat <- if (do_combat) {
  cc <- combat_train_apply(G_all, y_all, cohort_all, G_all, cohort_all)
  Gf <- cc$train
  top <- FS_top(Gf, varcut, y_all); Xf <- Gf[, top, drop = FALSE]
  if (do_cov) Xf <- cbind(batch_nyu_166 = as.integer(cohort_all == "nyu_166"),
                          batch_nyu_82  = as.integer(cohort_all == "nyu_82"), Xf)
  list(Xtr = Xf, Xva = Xf, ytr = y_all, yva = y_all,
       batch_cols = if (do_cov) c("batch_nyu_166", "batch_nyu_82") else character(0))
} else {
  build_fold_data(seq_along(y_all), seq_along(y_all))
}
final_fit <- train_predict(final_dat)
saveRDS(list(model = final_fit$model, learner = learner, variant = variant, varcut = varcut,
             fs_method = fs_method, features = colnames(final_dat$Xtr), training_n = nrow(mat),
             apparent_auc = final_fit$auc, cv_honest_combat = do_combat, created = Sys.time()),
        file.path("res", sprintf("%s_%s_v%d.rds", base, learner, varcut)))

cat(sprintf("Completed variant=%s %s v%d fs=%s: mean CV AUC %.4f (apparent %.4f)\n",
            variant, learner, varcut, fs_method, mean(cvtotal$AUC, na.rm = TRUE), final_fit$auc))
