#!/usr/bin/env Rscript
# Build a unified per-model feature+coefficient table for the transfer models,
# both FS methods, used for the Venn (feature-set overlap across the 3 cohorts)
# and feature-importance (consensus-feature |coef|) panels.
#   mad -> extract coef() from the saved fixed-whole-feature GLM models
#   cor -> already computed (corr_transfer_model_features.csv)
# Output: data/transfer_features_long.csv  (fs_method, cohort, varcut, feature, coef)
suppressPackageStartupMessages(library(glmnet))
D <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260626_full"
out <- "/Volumes/hua_mac/research/aris/harvey/mirna_202605/20260715/data"
sources <- c("EV_RNA","nyu_166","nyu_82"); varcuts <- c(5,10,20,50,100,200,500,1000)

rows <- list()
# --- MAD: coefficients from saved models ---
for (co in sources) for (vc in varcuts) {
  f <- file.path(D, "fork2_results/models", sprintf("%s_miRNA_glm_v%d_fixed_whole_features.rds", co, vc))
  if (!file.exists(f)) { cat("MISSING mad model", co, vc, "\n"); next }
  m <- readRDS(f); cm <- coef(m, s = "lambda.min")
  feats <- rownames(cm)[-1]; coefs <- as.numeric(cm[-1, 1])
  rows[[length(rows)+1]] <- data.frame(fs_method="mad", cohort=co, varcut=vc,
                                       feature=feats, coef=coefs, stringsAsFactors=FALSE)
}
# --- COR: from the already-saved feature/coef CSV ---
cr <- read.csv(file.path(D, "mirna_20260714/results/correlation/transfer/corr_transfer_model_features.csv"), stringsAsFactors=FALSE)
names(cr) <- tolower(names(cr))  # Version,Source,Feature,Coefficient
rows[[length(rows)+1]] <- data.frame(fs_method="cor", cohort=cr$source, varcut=cr$version,
                                     feature=cr$feature, coef=cr$coefficient, stringsAsFactors=FALSE)
d <- do.call(rbind, rows)
d$coef <- round(d$coef, 6)
write.csv(d, file.path(out, "transfer_features_long.csv"), row.names=FALSE)
cat("wrote transfer_features_long.csv rows", nrow(d),
    "| mad", sum(d$fs_method=="mad"), "cor", sum(d$fs_method=="cor"), "\n")
# quick check: selected-set sizes per (fs, cohort) at v20
for (fs in c("mad","cor")) { s <- d[d$fs_method==fs & d$varcut==20, ]
  cat(fs, "v20 set sizes:", paste(sapply(sources, function(c) sum(s$cohort==c)), collapse="/"), "\n") }
