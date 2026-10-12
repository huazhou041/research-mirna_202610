#!/usr/bin/env Rscript
# Transfer analysis (3x3 cross-cohort) for the dashboard, both FS methods.
# For each (fs_method, source, varcut): refit a GLM on the WHOLE source cohort
# using the top features (mad = variance, cor = |corr with recurrence|), apply to
# each other cohort -> AUC + ROC curve (single prediction set per target, so the
# ROC is that transfer's actual curve; AUC is its area, direction="<").
# Writes transfer_auc_long.csv + transfer_roc_long.csv.
suppressPackageStartupMessages({ library(glmnet); library(pROC); library(caret); library(jsonlite) })

base   <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260626_full"
rdata  <- file.path(base, "fork2_results/rdata")
outdir <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260715/data"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sources <- c("EV_RNA","nyu_166","nyu_82"); varcuts <- c(5,10,20,50,100,200,500,1000)
alphas  <- fromJSON(file.path(base, "fork2_new/codes/hyperparams.json"), simplifyVector = TRUE)
fmt <- function(x) paste(sprintf("%.4f", x), collapse = "|")

dl <- lapply(setNames(file.path(rdata, paste0(sources, "_miRNA.RDS")), sources), function(p) {
  x <- readRDS(p); x <- x[!is.na(x$recurrence), , drop = FALSE]; x$recurrence <- as.numeric(as.character(x$recurrence)); x
})
select_feat <- function(dat, vc, fs, y = dat$recurrence) {
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
# Build the curve at the EMPIRICAL thresholds from pROC, not a fixed 0..1 probability
# grid: transferred predictions can saturate/bunch, which a fixed grid cannot resolve
# (it collapsed some curves to <12 points and drifted up to 0.16 from the true AUC).
# Using pROC's own operating points guarantees curve area == reported AUC.
curve_of <- function(r) {
  fpr <- 1 - r$specificities; tpr <- r$sensitivities
  o <- order(fpr, tpr)
  list(fpr = fmt(fpr[o]), tpr = fmt(tpr[o]))
}

auc_rows <- list(); roc_rows <- list(); pred_rows <- list()
for (fs in c("mad","cor")) for (src in sources) for (vc in varcuts) {
  feats <- select_feat(dl[[src]], vc, fs)
  m <- fit_glm(as.matrix(dl[[src]][, feats, drop=FALSE]), dl[[src]]$recurrence)
  for (tg in setdiff(sources, src)) {
    td <- dl[[tg]]; y <- td$recurrence
    p <- as.numeric(predict(m, as.matrix(td[, feats, drop=FALSE]), s="lambda.min", type="response"))
    r <- pROC::roc(y, p, quiet=TRUE, direction="<")
    auc <- as.numeric(r$auc)
    auc_rows[[length(auc_rows)+1]] <- data.frame(fs_method=fs, source=src, target=tg, varcut=vc,
                                                 auc=round(auc,4), n_target=length(y), stringsAsFactors=FALSE)
    cur <- curve_of(r)
    roc_rows[[length(roc_rows)+1]] <- data.frame(fs_method=fs, source=src, target=tg, varcut=vc,
                                                 auc=round(auc,4), n_pred=length(y), fpr=cur$fpr, tpr=cur$tpr, stringsAsFactors=FALSE)
    # keep raw predictions so curves can be rebuilt without refitting
    pred_rows[[length(pred_rows)+1]] <- data.frame(fs_method=fs, source=src, target=tg, varcut=vc,
                                                   subj=rownames(td), y=y, prob=p, stringsAsFactors=FALSE)
  }
  cat("done", fs, src, "v", vc, "\n")
}
write.csv(do.call(rbind, auc_rows), file.path(outdir, "transfer_auc_long.csv"), row.names=FALSE)
write.csv(do.call(rbind, roc_rows), file.path(outdir, "transfer_roc_long.csv"), row.names=FALSE)
write.csv(do.call(rbind, pred_rows), file.path(outdir, "transfer_pred_long.csv"), row.names=FALSE)
cat("wrote transfer_auc_long.csv (", length(auc_rows), ") transfer_roc_long.csv (", length(roc_rows), ")\n")
