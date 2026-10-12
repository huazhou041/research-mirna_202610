library(tidyverse)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

out_dir = '20260626_full/results'
summary_dir = file.path(out_dir, 'summary')
plot_dir = file.path(out_dir, 'plots')
archive_dir = file.path(out_dir, 'Archive')

dir.create(archive_dir, recursive=TRUE, showWarnings=FALSE)

summary_df = read.csv(file.path(summary_dir, 'edrn_all_models_summary.csv')) %>%
  select(Version, N_Features, N_Nonzero, EDRN_AUC, AUC_166, AUC_82)

write.csv(summary_df,
          file.path(archive_dir, 'edrn_all_models_summary.csv'),
          row.names=FALSE, quote=FALSE)

plot_files = c(
  'venn_top20_3models.pdf',
  'venn_top50_3models.pdf',
  'feature_importance_top20_3models.pdf',
  'feature_importance_top50_3models.pdf'
)

file.copy(file.path(plot_dir, plot_files),
          file.path(archive_dir, plot_files),
          overwrite=TRUE)

cat(sprintf('Archive outputs copied to %s\n', archive_dir))
