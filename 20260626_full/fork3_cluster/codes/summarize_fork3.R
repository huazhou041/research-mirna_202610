#!/usr/bin/env Rscript
# Summarize the fork3 CV-honest ComBat runs and compare, cell for cell, against
# the fork2 leaky-ComBat numbers (both already use pROC direction="<", so the
# ONLY difference is fold-internal vs whole-dataset batch correction = the
# leakage effect). Run from fork3 root after the array completes.
suppressPackageStartupMessages({ library(tidyverse) })

variants <- c("plain", "cov", "combat", "both")
learners <- c("GLM", "RF", "XGB")
varcuts  <- c(5, 10, 20, 50, 100, 200, 500, 1000)

# ---- fork3 CV-honest means from raw_pred cvSummaries ----
rows <- list()
missing <- character(0)
for (v in variants) for (l in learners) for (vc in varcuts) {
  f <- file.path("raw_pred", sprintf("combined_%s_miRNA_%s_v%d_cvSummary.csv", v, l, vc))
  if (!file.exists(f)) { missing <- c(missing, basename(f)); next }
  a <- read.csv(f)$AUC
  rows[[length(rows) + 1]] <- data.frame(variant = v, learner = l, varcut = vc,
                                         cvhonest = mean(a, na.rm = TRUE), n = sum(is.finite(a)))
}
fork3 <- bind_rows(rows)
if (length(missing)) cat("MISSING fork3 cells:", length(missing), "->", paste(head(missing, 20), collapse=", "), "\n")

# ---- wide 4-way table (matches combined_4way_comparison.csv layout) ----
wide <- fork3 %>%
  select(learner, variant, varcut, cvhonest) %>%
  pivot_wider(names_from = varcut, values_from = cvhonest, names_prefix = "v") %>%
  arrange(match(learner, learners), match(variant, variants))
write.csv(wide, "results/fork3_combined_4way_cvhonest.csv", row.names = FALSE)

# ---- side-by-side vs fork2 leaky-ComBat (direction-fixed) ----
f2_path <- "../fork2/results/combined_models_fixed_direction.csv"
if (file.exists(f2_path)) {
  f2 <- read.csv(f2_path) %>%
    transmute(variant, learner, varcut, leaky = fixed_mean)
  cmp <- fork3 %>%
    inner_join(f2, by = c("variant", "learner", "varcut")) %>%
    mutate(delta = cvhonest - leaky) %>%
    arrange(match(variant, variants), match(learner, learners), varcut)
  write.csv(cmp, "results/fork3_vs_fork2_combat_leakage_comparison.csv", row.names = FALSE)

  cat("\n=== leakage effect (fork3 CV-honest minus fork2 leaky-ComBat), corrected variants only ===\n")
  cmp %>% filter(variant %in% c("combat", "both")) %>%
    group_by(variant, learner) %>%
    summarise(mean_delta = mean(delta), min_delta = min(delta), max_delta = max(delta), .groups = "drop") %>%
    as.data.frame() %>% print(digits = 3)
  cat("\nplain/cov (should be ~0, sanity check that no-correction path is unchanged):\n")
  cmp %>% filter(variant %in% c("plain", "cov")) %>%
    summarise(max_abs_delta = max(abs(delta))) %>% as.data.frame() %>% print(digits = 3)
} else {
  cat("NOTE: fork2 comparison file not found at", f2_path, "- wrote CV-honest table only.\n")
}
cat("\nwrote results/fork3_combined_4way_cvhonest.csv and results/fork3_vs_fork2_combat_leakage_comparison.csv\n")
