library(tidyverse)
library(ggplot2)
library(glmnet)
library(pROC)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

cat('=== FIXING Feature Importance: Consensus Features with 3 Bars per Feature ===\n\n')

out_dir = 'results/codes/results'
plot_dir = file.path(out_dir, 'plots')
summary_dir = file.path(out_dir, 'summary')
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(summary_dir, recursive=TRUE, showWarnings=FALSE)

make_foldid = function(y, nfolds=5) {
  set.seed(20260618)
  sample(rep(seq_len(nfolds), length.out=length(y)))
}

# Load data
d_166 = readRDS('data/model_data/nyu_166_miRNA.RDS')
d_82 = readRDS('data/model_data/nyu_82_miRNA.RDS')

d_166_clean = d_166[!is.na(d_166$recurrence), ]
d_82_clean = d_82[!is.na(d_82$recurrence), ]

# ============================================================
# TOP 20 FEATURES - CONSENSUS
# ============================================================
cat('Processing TOP 20 - Consensus Features...\n\n')

# Get top 20 for each cohort
gene_cols_166 = setdiff(colnames(d_166_clean), 'recurrence')
mad_vals_166 = apply(d_166_clean[, gene_cols_166], 2, mad, na.rm=TRUE)
top_20_genes_166 = gene_cols_166[order(mad_vals_166, decreasing=TRUE)[1:20]]

gene_cols_82 = setdiff(colnames(d_82_clean), 'recurrence')
mad_vals_82 = apply(d_82_clean[, gene_cols_82], 2, mad, na.rm=TRUE)
top_20_genes_82 = gene_cols_82[order(mad_vals_82, decreasing=TRUE)[1:20]]

# EDRN v20
edrn_v20 = readRDS('prev/models/EV_RNA_miRNA_glm_v20.rds')
coef_edrn_20 = coef(edrn_v20, s='lambda.min')
edrn_20_feat = rownames(coef_edrn_20)[-1]
coef_edrn_20_vec = as.vector(coef_edrn_20[-1])
names(coef_edrn_20_vec) = edrn_20_feat

# Train GLM models to get coefficients
y_166 = d_166_clean$recurrence
X_166_20 = as.matrix(d_166_clean[, top_20_genes_166])
m_166_20 = cv.glmnet(X_166_20, y_166, family='binomial', nfolds=5, alpha=0.5,
                     foldid=make_foldid(y_166, nfolds=5))
coef_166_20 = as.vector(coef(m_166_20, s='lambda.min')[-1])
names(coef_166_20) = top_20_genes_166

y_82 = d_82_clean$recurrence
X_82_20 = as.matrix(d_82_clean[, top_20_genes_82])
m_82_20 = cv.glmnet(X_82_20, y_82, family='binomial', nfolds=5, alpha=0.5,
                    foldid=make_foldid(y_82, nfolds=5))
coef_82_20 = as.vector(coef(m_82_20, s='lambda.min')[-1])
names(coef_82_20) = top_20_genes_82

# Find consensus features (in all 3 models)
consensus_20 = intersect(intersect(top_20_genes_166, top_20_genes_82), edrn_20_feat)

cat('Consensus features (top 20, in all 3 models):', length(consensus_20), '\n\n')

# Create table with all 3 coefficients for consensus features
plot_data_20_consensus = data.frame(
  Feature = consensus_20,
  coef_166 = abs(coef_166_20[consensus_20]),
  coef_82 = abs(coef_82_20[consensus_20]),
  coef_EDRN = abs(coef_edrn_20_vec[consensus_20]),
  stringsAsFactors = FALSE
) %>%
  mutate(Mean_AbsCoef = (coef_166 + coef_82 + coef_EDRN) / 3) %>%
  arrange(desc(Mean_AbsCoef)) %>%
  head(15)  # Top 15 consensus features

write.csv(plot_data_20_consensus,
          file.path(summary_dir, 'feature_importance_top20_3models.csv'),
          row.names=FALSE, quote=FALSE)

# Reshape for plotting (3 bars per feature)
plot_20_long = plot_data_20_consensus %>%
  select(-Mean_AbsCoef) %>%
  pivot_longer(cols=starts_with('coef_'),
               names_to='Model', values_to='AbsCoef') %>%
  mutate(Model = gsub('coef_', '', Model)) %>%
  mutate(Model = factor(Model, levels=c('166', '82', 'EDRN')))

p_imp_20_consensus = ggplot(plot_20_long,
                             aes(x=reorder(Feature, -AbsCoef, FUN=mean), y=AbsCoef, fill=Model)) +
  geom_col(position='dodge', width=0.7, alpha=0.85) +
  scale_fill_manual(values=c('166'='#E41A1C', '82'='#377EB8', 'EDRN'='#4DAF4A')) +
  coord_flip() +
  labs(
    title='Top 15 Consensus Features: Coefficient Comparison (Top 20 Analysis)',
    subtitle='Features appearing in all 3 models - 3 bars per feature show model differences',
    x='Feature',
    y='|Coefficient|',
    fill='Model'
  ) +
  theme_bw() +
  theme(axis.text.y=element_text(size=9),
        plot.title=element_text(face='bold', size=12),
        plot.subtitle=element_text(size=10))

ggsave(file.path(plot_dir, 'feature_importance_top20_3models.pdf'), p_imp_20_consensus, width=13, height=9)
cat('✓ Saved: feature_importance_top20_3models.pdf\n')
cat(sprintf('  Showing %d consensus features with 3 bars each\n\n', nrow(plot_data_20_consensus)))

# ============================================================
# TOP 50 FEATURES - CONSENSUS
# ============================================================
cat('Processing TOP 50 - Consensus Features...\n\n')

# Get top 50 for each cohort
top_50_genes_166 = gene_cols_166[order(mad_vals_166, decreasing=TRUE)[1:50]]
top_50_genes_82 = gene_cols_82[order(mad_vals_82, decreasing=TRUE)[1:50]]

# EDRN v50
edrn_v50 = readRDS('prev/models/EV_RNA_miRNA_glm_v50.rds')
coef_edrn_50 = coef(edrn_v50, s='lambda.min')
edrn_50_feat = rownames(coef_edrn_50)[-1]
coef_edrn_50_vec = as.vector(coef_edrn_50[-1])
names(coef_edrn_50_vec) = edrn_50_feat

# Train GLM models to get coefficients
X_166_50 = as.matrix(d_166_clean[, top_50_genes_166])
m_166_50 = cv.glmnet(X_166_50, y_166, family='binomial', nfolds=5, alpha=0.5,
                     foldid=make_foldid(y_166, nfolds=5))
coef_166_50 = as.vector(coef(m_166_50, s='lambda.min')[-1])
names(coef_166_50) = top_50_genes_166

X_82_50 = as.matrix(d_82_clean[, top_50_genes_82])
m_82_50 = cv.glmnet(X_82_50, y_82, family='binomial', nfolds=5, alpha=0.5,
                    foldid=make_foldid(y_82, nfolds=5))
coef_82_50 = as.vector(coef(m_82_50, s='lambda.min')[-1])
names(coef_82_50) = top_50_genes_82

# Find consensus features (in all 3 models)
consensus_50 = intersect(intersect(top_50_genes_166, top_50_genes_82), edrn_50_feat)

cat('Consensus features (top 50, in all 3 models):', length(consensus_50), '\n\n')

# Create table with all 3 coefficients for consensus features
plot_data_50_consensus = data.frame(
  Feature = consensus_50,
  coef_166 = abs(coef_166_50[consensus_50]),
  coef_82 = abs(coef_82_50[consensus_50]),
  coef_EDRN = abs(coef_edrn_50_vec[consensus_50]),
  stringsAsFactors = FALSE
) %>%
  mutate(Mean_AbsCoef = (coef_166 + coef_82 + coef_EDRN) / 3) %>%
  arrange(desc(Mean_AbsCoef)) %>%
  head(18)  # Top 18 consensus features

write.csv(plot_data_50_consensus,
          file.path(summary_dir, 'feature_importance_top50_3models.csv'),
          row.names=FALSE, quote=FALSE)

# Reshape for plotting (3 bars per feature)
plot_50_long = plot_data_50_consensus %>%
  select(-Mean_AbsCoef) %>%
  pivot_longer(cols=starts_with('coef_'),
               names_to='Model', values_to='AbsCoef') %>%
  mutate(Model = gsub('coef_', '', Model)) %>%
  mutate(Model = factor(Model, levels=c('166', '82', 'EDRN')))

p_imp_50_consensus = ggplot(plot_50_long,
                             aes(x=reorder(Feature, -AbsCoef, FUN=mean), y=AbsCoef, fill=Model)) +
  geom_col(position='dodge', width=0.7, alpha=0.85) +
  scale_fill_manual(values=c('166'='#E41A1C', '82'='#377EB8', 'EDRN'='#4DAF4A')) +
  coord_flip() +
  labs(
    title='Top 18 Consensus Features: Coefficient Comparison (Top 50 Analysis)',
    subtitle='Features appearing in all 3 models - 3 bars per feature show model differences',
    x='Feature',
    y='|Coefficient|',
    fill='Model'
  ) +
  theme_bw() +
  theme(axis.text.y=element_text(size=9),
        plot.title=element_text(face='bold', size=12),
        plot.subtitle=element_text(size=10))

ggsave(file.path(plot_dir, 'feature_importance_top50_3models.pdf'), p_imp_50_consensus, width=13, height=10)
cat('✓ Saved: feature_importance_top50_3models.pdf\n')
cat(sprintf('  Showing %d consensus features with 3 bars each\n\n', nrow(plot_data_50_consensus)))

cat('=== FIXED PLOTS (NOW WITH CONSENSUS FEATURES) ===\n')
cat('✓ feature_importance_top20_3models.pdf\n')
cat('  - Shows consensus features (in all 3 models)\n')
cat('  - 3 bars per feature (Red=nyu_166, Blue=nyu_82, Green=EDRN)\n')
cat('  - Directly comparable coefficients across models\n\n')

cat('✓ feature_importance_top50_3models.pdf\n')
cat('  - Shows consensus features (in all 3 models)\n')
cat('  - 3 bars per feature (Red=nyu_166, Blue=nyu_82, Green=EDRN)\n')
cat('  - Directly comparable coefficients across models\n')
