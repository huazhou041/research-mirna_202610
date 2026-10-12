# Re-score each source's OWN-cohort bootstrap CV AUC (the "Source_CV_AUC" column
# in all_sources_fixed_models_on_other2_summary.csv) using fixed pROC direction="<"
# instead of the original roc(y, y_hat) call (deep-drug-response/tune_drug.R,
# validate_glmnet(), no direction specified -- same vulnerable pattern found
# elsewhere). No retraining needed: predictions already saved per fold x boot.
suppressPackageStartupMessages({ library(pROC); library(data.table) })

fixed_auc <- function(y, prob) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

varcuts <- c(5, 10, 20, 50, 100, 200, 500, 1000)

sources <- list(
  EV_RNA  = list(pred_dir = "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/res",
                 rds = "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2/rdata/EV_RNA_miRNA.RDS"),
  nyu_166 = list(pred_dir = "/gpfs/scratch/zhouh05/aris/segal/val_20260605/res",
                 rds = "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2/rdata/nyu_166_miRNA.RDS"),
  nyu_82  = list(pred_dir = "/gpfs/scratch/zhouh05/aris/segal/val_20260605/res",
                 rds = "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2/rdata/nyu_82_miRNA.RDS")
)

out_rows <- list()
for (src in names(sources)) {
  info <- sources[[src]]
  mat <- readRDS(info$rds)
  mat <- mat[!is.na(mat$recurrence), , drop = FALSE]
  y_all <- setNames(as.numeric(as.character(mat$recurrence)), rownames(mat))

  for (vc in varcuts) {
    pf <- file.path(info$pred_dir, sprintf("%s_miRNA_glm_v%d_pred_prob.csv", src, vc))
    cf <- file.path(info$pred_dir, sprintf("%s_miRNA_glm_v%d_cvSummary.csv", src, vc))
    if (!file.exists(pf) || !file.exists(cf)) { cat("MISSING:", pf, "\n"); next }
    pred <- fread(pf, showProgress = FALSE)
    cvs  <- fread(cf, showProgress = FALSE)
    n_unmatched <- sum(!pred$subj %in% names(y_all))
    pred$y <- y_all[as.character(pred$subj)]

    grp <- split(pred, list(pred$foldid, pred$bootid), drop = TRUE)
    aucs <- vapply(grp, function(g) fixed_auc(g$y, g$prob), numeric(1))
    aucs <- aucs[is.finite(aucs)]
    orig <- cvs$AUC[is.finite(cvs$AUC)]

    out_rows[[length(out_rows) + 1]] <- data.frame(
      source = src, varcut = vc, n_unmatched_subj = n_unmatched,
      orig_mean = if (length(orig)) mean(orig) else NA_real_,
      fixed_mean = if (length(aucs)) mean(aucs) else NA_real_,
      n = length(aucs), stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, out_rows)
res$diff <- res$fixed_mean - res$orig_mean
print(res, digits = 4)
cat("\ncells with |diff|>0.02:", sum(abs(res$diff) > 0.02, na.rm = TRUE), "of", nrow(res), "\n")
cat("max |diff|:", max(abs(res$diff), na.rm = TRUE), "\n")
cat("any unmatched subj (id mismatch risk):", sum(res$n_unmatched_subj), "\n")
write.csv(res, "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2/results/own_cv_fixed_direction.csv", row.names = FALSE)
