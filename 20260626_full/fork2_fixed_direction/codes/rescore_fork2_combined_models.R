# Re-score tune_combined_models.R (GLM/RF/XGB fixed-hyperparameter combined runs)
# from saved raw_pred/*_pred_prob.csv using fixed pROC direction="<", instead of the
# original AUC (computed with direction="auto", which can silently flip per fold and
# inflate AUC on folds where predictions are actually anti-correlated with truth).
# No retraining needed. Run from fork2 root (where raw_pred/, res/, rdata/ live).
suppressPackageStartupMessages({ library(pROC); library(data.table) })

variants <- c("plain", "cov", "combat", "both")
learners <- c("GLM", "RF", "XGB")
varcuts  <- c(5, 10, 20, 50, 100, 200, 500, 1000)

fixed_auc <- function(y, prob) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

out_rows <- list()
n_missing <- 0

for (variant in variants) {
  rds_f <- file.path("rdata", sprintf("combined_%s_miRNA.RDS", variant))
  if (!file.exists(rds_f)) { next }
  mat <- readRDS(rds_f)
  y_all <- setNames(as.numeric(as.character(mat$recurrence)), rownames(mat))
  base <- sprintf("combined_%s_miRNA", variant)

  for (learner in learners) for (vc in varcuts) {
    pf <- file.path("raw_pred", sprintf("%s_%s_v%d_pred_prob.csv", base, learner, vc))
    cf <- file.path("raw_pred", sprintf("%s_%s_v%d_cvSummary.csv", base, learner, vc))
    if (!file.exists(pf) || !file.exists(cf)) { n_missing <- n_missing + 1; next }
    pred <- tryCatch(fread(pf, showProgress = FALSE), error = function(e) NULL)
    cvs  <- tryCatch(fread(cf, showProgress = FALSE), error = function(e) NULL)
    if (is.null(pred) || !nrow(pred) || is.null(cvs) || !("AUC" %in% names(cvs))) { n_missing <- n_missing + 1; next }
    pred$y <- y_all[pred$subj]

    grp <- split(pred, list(pred$foldid, pred$bootid), drop = TRUE)
    aucs <- vapply(grp, function(g) fixed_auc(g$y, g$prob), numeric(1))
    aucs <- aucs[is.finite(aucs)]
    orig <- cvs$AUC[is.finite(cvs$AUC)]

    out_rows[[length(out_rows) + 1]] <- data.frame(
      variant = variant, learner = learner, varcut = vc,
      orig_mean = if (length(orig)) mean(orig) else NA_real_,
      fixed_mean = if (length(aucs)) mean(aucs) else NA_real_,
      n = length(aucs),
      stringsAsFactors = FALSE)
  }
  cat("done:", variant, " ", format(Sys.time()), "\n")
}

res <- do.call(rbind, out_rows)
res$diff <- res$fixed_mean - res$orig_mean
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/combined_models_fixed_direction.csv", row.names = FALSE)
cat("wrote results/combined_models_fixed_direction.csv rows=", nrow(res), " missing=", n_missing, "\n")
cat("cells with |diff| > 0.02:", sum(abs(res$diff) > 0.02, na.rm = TRUE), "of", sum(is.finite(res$diff)), "\n")
cat("max |diff|:", max(abs(res$diff), na.rm = TRUE), "\n")
