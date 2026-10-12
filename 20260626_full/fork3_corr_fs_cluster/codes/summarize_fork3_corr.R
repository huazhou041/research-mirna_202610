#!/usr/bin/env Rscript
# Summarize the fork3_corr_fs runs (CV-honest ComBat + supervised correlation
# feature selection) and compare, cell for cell, against the MAD-selection fork3
# results. Both use fold-internal selection + direction="<", so the difference
# isolates the effect of supervised (corr) vs unsupervised (MAD) feature ranking.
# Run from fork3_corr_fs root.
suppressPackageStartupMessages(library(tidyverse))

variants <- c("plain", "cov", "combat", "both")
learners <- c("GLM", "RF", "XGB")
varcuts  <- c(5, 10, 20, 50, 100, 200, 500, 1000)

rows <- list(); missing <- character(0)
for (v in variants) for (l in learners) for (vc in varcuts) {
  f <- file.path("raw_pred", sprintf("combined_%s_miRNA_%s_v%d_cvSummary.csv", v, l, vc))
  if (!file.exists(f)) { missing <- c(missing, basename(f)); next }
  a <- read.csv(f)$AUC
  rows[[length(rows) + 1]] <- data.frame(variant = v, learner = l, varcut = vc,
                                         corr_fs = mean(a, na.rm = TRUE), n = sum(is.finite(a)))
}
corr <- bind_rows(rows)
if (length(missing)) cat("MISSING corr_fs cells:", length(missing), "->", paste(head(missing, 20), collapse=", "), "\n")

wide <- corr %>% select(learner, variant, varcut, corr_fs) %>%
  pivot_wider(names_from = varcut, values_from = corr_fs, names_prefix = "v") %>%
  arrange(match(learner, learners), match(variant, variants))
write.csv(wide, "results/fork3_corr_fs_4way.csv", row.names = FALSE)

# compare vs MAD-selection fork3 (leakage comparison file has the CV-honest MAD means)
mad_path <- "../fork3/results/fork3_vs_fork2_combat_leakage_comparison.csv"
if (file.exists(mad_path)) {
  mad <- read.csv(mad_path) %>% transmute(variant, learner, varcut, mad_fs = cvhonest)
  cmp <- corr %>% inner_join(mad, by = c("variant","learner","varcut")) %>%
    mutate(delta_corr_minus_mad = corr_fs - mad_fs) %>%
    arrange(match(variant, variants), match(learner, learners), varcut)
  write.csv(cmp, "results/corr_fs_vs_mad_fs_comparison.csv", row.names = FALSE)
  cat("\n=== corr-FS minus MAD-FS (mean CV AUC), by learner x varcut (avg over variants) ===\n")
  cmp %>% group_by(learner, varcut) %>% summarise(delta = mean(delta_corr_minus_mad), .groups="drop") %>%
    pivot_wider(names_from = varcut, values_from = delta, names_prefix="v") %>% as.data.frame() %>% print(digits=2)
} else {
  cat("NOTE: MAD comparison file not found at", mad_path, "\n")
}
cat("\nwrote results/fork3_corr_fs_4way.csv and results/corr_fs_vs_mad_fs_comparison.csv\n")
