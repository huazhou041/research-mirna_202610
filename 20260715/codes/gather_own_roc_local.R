#!/usr/bin/env Rscript
# Own-cohort ROC curves for the dashboard (replaces transferred ROC in the ROC panel).
# For each (fs_method, cohort, varcut): 10-fold CV with FOLD-INTERNAL feature selection
# (leakage-proof, mirrors the Source_CV_AUC methodology), collect out-of-fold predictions,
# build one pooled ROC curve (direction="<"). AUC label = pooled OOF AUC.
suppressPackageStartupMessages({ library(glmnet); library(pROC); library(caret); library(jsonlite) })
base   <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260626_full"
rdata  <- file.path(base, "fork2_results/rdata")
outdir <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260715/data"
sources <- c("EV_RNA","nyu_166","nyu_82"); varcuts <- c(5,10,20,50,100,200,500,1000)
alphas  <- fromJSON(file.path(base, "fork2_new/codes/hyperparams.json"), simplifyVector = TRUE)
fmt <- function(x) paste(sprintf("%.4f", x), collapse = "|")

dl <- lapply(setNames(file.path(rdata, paste0(sources, "_miRNA.RDS")), sources), function(p) {
  x <- readRDS(p); x <- x[!is.na(x$recurrence), , drop = FALSE]; x$recurrence <- as.numeric(as.character(x$recurrence)); x
})
select_feat <- function(dat, vc, fs, y) {
  genes <- setdiff(colnames(dat), "recurrence")
  score <- if (fs == "cor") abs(as.numeric(suppressWarnings(cor(as.matrix(dat[, genes]), y)))) else apply(dat[, genes], 2, mad, na.rm = TRUE)
  score[is.na(score)] <- -Inf
  genes[order(score, decreasing = TRUE)[seq_len(min(vc, length(genes)))]]
}
fit_glm <- function(x, y) {
  set.seed(123); fid <- createFolds(y, k = 10, list = FALSE)
  ms <- lapply(alphas, function(a) tryCatch(cv.glmnet(x, y, family="binomial", alpha=a, foldid=fid),
              error=function(e) cv.glmnet(x, y, family="binomial", alpha=a, foldid=fid, lambda=exp(seq(log(.001),log(5),length.out=100)))))
  ms[[which.min(sapply(ms, function(m) m$cvm[which.min(abs(m$lambda-m$lambda.min))]))]]
}
# Curve at pROC's EMPIRICAL thresholds (not a fixed 0..1 grid) so the plotted
# curve's area always equals the reported AUC. See gather_transfer_local.R.
curve_of <- function(r) {
  fpr <- 1 - r$specificities; tpr <- r$sensitivities
  o <- order(fpr, tpr)
  list(fpr = fmt(fpr[o]), tpr = fmt(tpr[o]))
}

roc_rows <- list(); pred_rows <- list()
for (fs in c("mad","cor")) for (co in sources) {
  dat <- dl[[co]]; y <- dat$recurrence
  set.seed(123); folds <- createFolds(y, k = 10, list = TRUE)
  for (vc in varcuts) {
    oof <- rep(NA_real_, length(y))
    for (f in folds) {
      tr <- setdiff(seq_along(y), f)
      feats <- select_feat(dat[tr, , drop=FALSE], vc, fs, y[tr])       # fold-internal FS
      m <- fit_glm(as.matrix(dat[tr, feats, drop=FALSE]), y[tr])
      oof[f] <- as.numeric(predict(m, as.matrix(dat[f, feats, drop=FALSE]), s="lambda.min", type="response"))
    }
    r <- pROC::roc(y, oof, quiet=TRUE, direction="<")
    auc <- as.numeric(r$auc)
    cur <- curve_of(r)
    roc_rows[[length(roc_rows)+1]] <- data.frame(fs_method=fs, cohort=co, varcut=vc,
      auc=round(auc,4), n_pred=length(y), fpr=cur$fpr, tpr=cur$tpr, stringsAsFactors=FALSE)
    # keep OOF predictions so curves can be rebuilt without refitting
    pred_rows[[length(pred_rows)+1]] <- data.frame(fs_method=fs, cohort=co, varcut=vc,
      subj=rownames(dat), y=y, prob=oof, stringsAsFactors=FALSE)
    cat("done", fs, co, "v", vc, "auc", round(auc,3), "\n")
  }
}
write.csv(do.call(rbind, roc_rows), file.path(outdir, "own_roc_long.csv"), row.names=FALSE)
write.csv(do.call(rbind, pred_rows), file.path(outdir, "own_pred_long.csv"), row.names=FALSE)
cat("wrote own_roc_long.csv (", length(roc_rows), ")\n")
