#!/usr/bin/env Rscript
# Gather combined 3-cohort model results (both FS methods) into two compact long
# CSVs for the dashboard, matching the erdn dashboard methodology:
#   * AUC   = per-fold mean AUC over the 100x10 bootstrap (from cvSummary), + 95% CI
#   * ROC   = POOLED curve across all fold x boot predictions (from pred_prob),
#             thresholds Inf..-Inf; fpr/tpr as pipe-delimited strings.
# fs_method: mad (fork3) | cor (fork3_corr_fs). Run on a compute node.
suppressPackageStartupMessages(library(data.table))

root <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/mirna_202606"
forks <- list(mad = file.path(root, "fork3"), cor = file.path(root, "fork3_corr_fs"))
variants <- c("plain","cov","combat","both"); learners <- c("GLM","RF","XGB")
varcuts  <- c(5,10,20,50,100,200,500,1000)
thresholds <- c(Inf, seq(1, 0, by = -0.01), -Inf)
fmt <- function(x) paste(sprintf("%.4f", x), collapse = "|")

# labels (shared across FS methods; same pooled data)
lab <- readRDS(file.path(forks$mad, "rdata/combined_plain_miRNA.RDS"))
y_all <- setNames(as.numeric(as.character(lab$recurrence)), rownames(lab))

pooled_curve <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  pos <- sum(y == 1); neg <- sum(y == 0)
  if (!pos || !neg) return(list(fpr=NA_character_, tpr=NA_character_))
  tpr <- fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) { pr <- p >= thresholds[i]; tpr[i] <- sum(pr & y==1)/pos; fpr[i] <- sum(pr & y==0)/neg }
  list(fpr = fmt(fpr), tpr = fmt(tpr))
}

auc_rows <- list(); roc_rows <- list()
for (fs in names(forks)) {
  rp <- file.path(forks[[fs]], "raw_pred")
  for (va in variants) for (le in learners) for (vc in varcuts) {
    base <- sprintf("combined_%s_miRNA_%s_v%d", va, le, vc)
    cf <- file.path(rp, paste0(base, "_cvSummary.csv")); pf <- file.path(rp, paste0(base, "_pred_prob.csv"))
    if (!file.exists(cf) || !file.exists(pf)) { cat("MISSING", fs, base, "\n"); next }
    a <- read.csv(cf)$AUC; a <- a[is.finite(a)]
    se <- sd(a)/sqrt(length(a))
    auc_rows[[length(auc_rows)+1]] <- data.frame(
      fs_method=fs, variant=va, learner=le, varcut=vc,
      mean=round(mean(a),4), sd=round(sd(a),4), n=length(a),
      lo=round(max(0,mean(a)-1.96*se),4), hi=round(min(1,mean(a)+1.96*se),4), stringsAsFactors=FALSE)
    pred <- fread(pf, showProgress=FALSE)
    cur <- pooled_curve(unname(y_all[as.character(pred$subj)]), as.numeric(pred$prob))
    roc_rows[[length(roc_rows)+1]] <- data.frame(
      fs_method=fs, variant=va, learner=le, varcut=vc,
      auc=round(mean(a),4), n_pred=nrow(pred), fpr=cur$fpr, tpr=cur$tpr, stringsAsFactors=FALSE)
  }
  cat("done fs=", fs, "\n")
}
dir.create("results", showWarnings=FALSE)
write.csv(do.call(rbind, auc_rows), "results/combined_auc_long.csv", row.names=FALSE)
write.csv(do.call(rbind, roc_rows), "results/combined_roc_long.csv", row.names=FALSE)
cat("wrote combined_auc_long.csv (", length(auc_rows), ") and combined_roc_long.csv (", length(roc_rows), ")\n")
