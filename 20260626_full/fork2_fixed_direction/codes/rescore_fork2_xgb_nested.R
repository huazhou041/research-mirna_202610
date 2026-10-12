# Re-score tune_combined_xgb_nested.R (nested-CV XGB) results from saved
# res/*_pred_prob.csv using fixed pROC direction="<" instead of the original
# AUC (direction="auto", can silently flip per fold). No retraining needed.
# Covers both the XGBcv (original 32-job run) and XGBcvFix (DMatrix spot-check)
# output prefixes. Run from fork2 root.
suppressPackageStartupMessages({ library(pROC); library(data.table) })

fixed_auc <- function(y, prob) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

pred_files <- list.files("res", pattern = "^combined_.*_XGBcv(Fix)?_v[0-9]+_pred_prob\\.csv$", full.names = TRUE)
rx <- "^combined_([a-z]+)_miRNA_(XGBcvFix|XGBcv)_v([0-9]+)_pred_prob\\.csv$"

out_rows <- list()
n_missing <- 0
y_cache <- list()

for (pf in pred_files) {
  fn <- basename(pf)
  m <- regmatches(fn, regexec(rx, fn))[[1]]
  if (length(m) != 4) { n_missing <- n_missing + 1; next }
  variant <- m[2]; tag <- m[3]; vc <- as.integer(m[4])
  cf <- sub("_pred_prob\\.csv$", "_cvSummary.csv", pf)
  if (!file.exists(cf)) { n_missing <- n_missing + 1; next }

  if (is.null(y_cache[[variant]])) {
    rds_f <- file.path("rdata", sprintf("combined_%s_miRNA.RDS", variant))
    mat <- readRDS(rds_f)
    y_cache[[variant]] <- setNames(as.numeric(as.character(mat$recurrence)), rownames(mat))
  }
  y_all <- y_cache[[variant]]

  pred <- tryCatch(fread(pf, showProgress = FALSE), error = function(e) NULL)
  cvs  <- tryCatch(fread(cf, showProgress = FALSE), error = function(e) NULL)
  if (is.null(pred) || !nrow(pred) || is.null(cvs) || !("AUC" %in% names(cvs))) { n_missing <- n_missing + 1; next }
  pred$y <- y_all[pred$subj]

  grp <- split(pred, list(pred$foldid, pred$bootid), drop = TRUE)
  aucs <- vapply(grp, function(g) fixed_auc(g$y, g$prob), numeric(1))
  aucs <- aucs[is.finite(aucs)]
  orig <- cvs$AUC[is.finite(cvs$AUC)]

  out_rows[[length(out_rows) + 1]] <- data.frame(
    variant = variant, tag = tag, varcut = vc,
    orig_mean = if (length(orig)) mean(orig) else NA_real_,
    fixed_mean = if (length(aucs)) mean(aucs) else NA_real_,
    n = length(aucs),
    stringsAsFactors = FALSE)
  cat("done:", fn, "\n")
}

res <- do.call(rbind, out_rows)
res$diff <- res$fixed_mean - res$orig_mean
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/xgb_nested_fixed_direction.csv", row.names = FALSE)
cat("wrote results/xgb_nested_fixed_direction.csv rows=", nrow(res), " missing=", n_missing, "\n")
cat("cells with |diff| > 0.02:", sum(abs(res$diff) > 0.02, na.rm = TRUE), "of", sum(is.finite(res$diff)), "\n")
cat("max |diff|:", max(abs(res$diff), na.rm = TRUE), "\n")
