# Overwrite the AUC column in each of the 24 own-CV cvSummary.csv copies (already
# staged under fork2/own_cv_fixed/, copied from mirna_202606/res and val_20260605/res)
# using fixed pROC direction="<", recomputed from the matching pred_prob.csv (source
# predictions untouched, still read from their original location). All writes stay
# under /gpfs/data -- no /gpfs/scratch quota involved.
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
out_dir <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606/fork2/own_cv_fixed"

for (src in names(sources)) {
  info <- sources[[src]]
  mat <- readRDS(info$rds)
  mat <- mat[!is.na(mat$recurrence), , drop = FALSE]
  y_all <- setNames(as.numeric(as.character(mat$recurrence)), rownames(mat))

  for (vc in varcuts) {
    pf <- file.path(info$pred_dir, sprintf("%s_miRNA_glm_v%d_pred_prob.csv", src, vc))
    cf <- file.path(out_dir, sprintf("%s_miRNA_glm_v%d_cvSummary.csv", src, vc))
    if (!file.exists(pf) || !file.exists(cf)) { cat("MISSING:", pf, cf, "\n"); next }
    pred <- fread(pf, showProgress = FALSE)
    cvs  <- fread(cf, showProgress = FALSE)
    pred$y <- y_all[as.character(pred$subj)]

    grp <- split(pred, list(pred$foldid, pred$bootid), drop = TRUE)
    fold_auc <- vapply(grp, function(g) fixed_auc(g$y, g$prob), numeric(1))
    fold_key <- names(grp)
    key_lookup <- setNames(fold_auc, fold_key)
    cvs_key <- paste(cvs$fold_id, cvs$bootid, sep = ".")
    stopifnot(all(cvs_key %in% names(key_lookup)))
    old_mean <- mean(cvs$AUC, na.rm = TRUE)
    cvs$AUC <- as.numeric(key_lookup[cvs_key])
    new_mean <- mean(cvs$AUC, na.rm = TRUE)
    fwrite(cvs, cf)
    cat(sprintf("%-8s v%-5d old_mean=%.4f new_mean=%.4f\n", src, vc, old_mean, new_mean))
  }
}
cat("done\n")
