library(tidyverse)
library(ggplot2)
library(ggvenn)
library(glmnet)
library(pROC)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

cat('=== REGENERATING: Venn Diagrams & Feature Importance (Top 20 & 50) ===\n\n')

out_dir = 'results/codes/results'
plot_dir = file.path(out_dir, 'plots')
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

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
# TOP 20 FEATURES ANALYSIS
# ============================================================
cat('Processing TOP 20 FEATURES...\n\n')

# Get top 20 for each new cohort
gene_cols_166 = setdiff(colnames(d_166_clean), 'recurrence')
mad_vals_166 = apply(d_166_clean[, gene_cols_166], 2, mad, na.rm=TRUE)
top_20_genes_166 = gene_cols_166[order(mad_vals_166, decreasing=TRUE)[1:20]]

gene_cols_82 = setdiff(colnames(d_82_clean), 'recurrence')
mad_vals_82 = apply(d_82_clean[, gene_cols_82], 2, mad, na.rm=TRUE)
top_20_genes_82 = gene_cols_82[order(mad_vals_82, decreasing=TRUE)[1:20]]

# Get EDRN v20 features
edrn_v20_model = readRDS('prev/models/EV_RNA_miRNA_glm_v20.rds')
coef_edrn_v20 = coef(edrn_v20_model, s='lambda.min')
edrn_v20_feat = rownames(coef_edrn_v20)[-1]

cat('Top 20 feature sets:\n')
cat('  nyu_166:', length(top_20_genes_166), '\n')
cat('  nyu_82: ', length(top_20_genes_82), '\n')
cat('  EDRN v20:', length(edrn_v20_feat), '\n\n')

# Venn diagram for TOP 20
venn_20_list = list(
  'nyu_166\n(top 20)' = top_20_genes_166,
  'nyu_82\n(top 20)' = top_20_genes_82,
  'EDRN v20\n(20 feat)' = edrn_v20_feat
)

p_venn_20 = ggvenn(venn_20_list,
                    fill_color = c('#E41A1C', '#377EB8', '#4DAF4A'),
                    stroke_size = 2.5,
                    set_name_size = 6,
                    text_size = 5) +
  labs(title = 'Feature Overlap: Top 20 Features (3 Models)',
       subtitle = 'nyu_166 vs nyu_82 vs EDRN v20') +
  theme_bw() +
  theme(plot.title = element_text(face='bold', size=14),
        plot.subtitle = element_text(size=11))

ggsave(file.path(plot_dir, 'venn_top20_3models.pdf'), p_venn_20, width=11, height=10)
cat('Saved: venn_top20_3models.pdf\n')

# Feature importance for TOP 20
# Train GLM on each using top 20
train_glm_top20 = function(data, features) {
  X = as.matrix(data[, features, drop=FALSE])
  y = data$recurrence
  m = cv.glmnet(X, y, family='binomial', nfolds=5, alpha=0.5,
                foldid=make_foldid(y, nfolds=5))
  coef_m = coef(m, s='lambda.min')
  coefs = as.vector(coef_m[-1])
  names(coefs) = features
  list(coefs = coefs, features = features)
}

glm_166_top20 = train_glm_top20(d_166_clean, top_20_genes_166)
glm_82_top20 = train_glm_top20(d_82_clean, top_20_genes_82)
edrn_v20_coefs = as.vector(coef(edrn_v20_model, s='lambda.min')[-1])

# Create comparison table for true consensus features
consensus_feat_20 = Reduce(intersect, list(glm_166_top20$features, glm_82_top20$features, edrn_v20_feat))

feat_20_table = data.frame(Feature = consensus_feat_20)
for (f in consensus_feat_20) {
  feat_20_table$coef_166[feat_20_table$Feature == f] =
    if (f %in% names(glm_166_top20$coefs)) abs(glm_166_top20$coefs[f]) else NA
  feat_20_table$coef_82[feat_20_table$Feature == f] =
    if (f %in% names(glm_82_top20$coefs)) abs(glm_82_top20$coefs[f]) else NA
  feat_20_table$coef_EDRN[feat_20_table$Feature == f] =
    if (f %in% edrn_v20_feat) abs(edrn_v20_coefs[which(edrn_v20_feat == f)]) else NA
}

feat_20_table = feat_20_table %>%
  mutate(Mean_AbsCoef = rowMeans(.[2:4], na.rm=TRUE)) %>%
  arrange(desc(Mean_AbsCoef)) %>%
  head(15)

plot_20 = feat_20_table %>%
  pivot_longer(cols=starts_with('coef_'),
               names_to='Model', values_to='AbsCoef') %>%
  mutate(Model = gsub('coef_', '', Model)) %>%
  mutate(Model = factor(Model, levels=c('166', '82', 'EDRN')))

p_imp_20 = ggplot(plot_20, aes(x=reorder(Feature, -AbsCoef, FUN=mean), y=AbsCoef, fill=Model)) +
  geom_col(position='dodge', width=0.7, alpha=0.8) +
  scale_fill_manual(values=c('166'='#E41A1C', '82'='#377EB8', 'EDRN'='#4DAF4A')) +
  coord_flip() +
  labs(
    title='Consensus Features: Coefficients (Top 20 Analysis)',
    subtitle='Features appearing in all 3 models',
    x='Feature',
    y='|Coefficient|',
    fill='Model'
  ) +
  theme_bw() +
  theme(axis.text.y=element_text(size=9), plot.title=element_text(face='bold', size=12))

ggsave(file.path(plot_dir, 'feature_importance_top20_3models.pdf'), p_imp_20, width=13, height=8)
cat('Saved: feature_importance_top20_3models.pdf\n')

# ============================================================
# TOP 50 FEATURES ANALYSIS
# ============================================================
cat('\n\nProcessing TOP 50 FEATURES...\n\n')

# Get top 50 (already have from earlier)
top_50_genes_166 = gene_cols_166[order(mad_vals_166, decreasing=TRUE)[1:50]]
top_50_genes_82 = gene_cols_82[order(mad_vals_82, decreasing=TRUE)[1:50]]

# Get EDRN v50 features
edrn_v50_model = readRDS('prev/models/EV_RNA_miRNA_glm_v50.rds')
coef_edrn_v50 = coef(edrn_v50_model, s='lambda.min')
edrn_v50_feat = rownames(coef_edrn_v50)[-1]

cat('Top 50 feature sets:\n')
cat('  nyu_166:', length(top_50_genes_166), '\n')
cat('  nyu_82: ', length(top_50_genes_82), '\n')
cat('  EDRN v50:', length(edrn_v50_feat), '\n\n')

# Venn diagram for TOP 50
venn_50_list = list(
  'nyu_166\n(top 50)' = top_50_genes_166,
  'nyu_82\n(top 50)' = top_50_genes_82,
  'EDRN v50\n(50 feat)' = edrn_v50_feat
)

p_venn_50 = ggvenn(venn_50_list,
                    fill_color = c('#E41A1C', '#377EB8', '#4DAF4A'),
                    stroke_size = 2.5,
                    set_name_size = 6,
                    text_size = 5) +
  labs(title = 'Feature Overlap: Top 50 Features (3 Models)',
       subtitle = 'nyu_166 vs nyu_82 vs EDRN v50') +
  theme_bw() +
  theme(plot.title = element_text(face='bold', size=14),
        plot.subtitle = element_text(size=11))

ggsave(file.path(plot_dir, 'venn_top50_3models.pdf'), p_venn_50, width=11, height=10)
cat('Saved: venn_top50_3models.pdf\n')

# Feature importance for TOP 50
train_glm_top50 = function(data, features) {
  X = as.matrix(data[, features, drop=FALSE])
  y = data$recurrence
  m = cv.glmnet(X, y, family='binomial', nfolds=5, alpha=0.5,
                foldid=make_foldid(y, nfolds=5))
  coef_m = coef(m, s='lambda.min')
  coefs = as.vector(coef_m[-1])
  names(coefs) = features
  list(coefs = coefs, features = features)
}

glm_166_top50 = train_glm_top50(d_166_clean, top_50_genes_166)
glm_82_top50 = train_glm_top50(d_82_clean, top_50_genes_82)
edrn_v50_coefs = as.vector(coef(edrn_v50_model, s='lambda.min')[-1])

# Create comparison table for true consensus features
consensus_feat_50 = Reduce(intersect, list(glm_166_top50$features, glm_82_top50$features, edrn_v50_feat))

feat_50_table = data.frame(Feature = consensus_feat_50)
for (f in consensus_feat_50) {
  feat_50_table$coef_166[feat_50_table$Feature == f] =
    if (f %in% names(glm_166_top50$coefs)) abs(glm_166_top50$coefs[f]) else NA
  feat_50_table$coef_82[feat_50_table$Feature == f] =
    if (f %in% names(glm_82_top50$coefs)) abs(glm_82_top50$coefs[f]) else NA
  feat_50_table$coef_EDRN[feat_50_table$Feature == f] =
    if (f %in% edrn_v50_feat) abs(edrn_v50_coefs[which(edrn_v50_feat == f)]) else NA
}

feat_50_table = feat_50_table %>%
  mutate(Mean_AbsCoef = rowMeans(.[2:4], na.rm=TRUE)) %>%
  arrange(desc(Mean_AbsCoef)) %>%
  head(20)

plot_50 = feat_50_table %>%
  pivot_longer(cols=starts_with('coef_'),
               names_to='Model', values_to='AbsCoef') %>%
  mutate(Model = gsub('coef_', '', Model)) %>%
  mutate(Model = factor(Model, levels=c('166', '82', 'EDRN')))

p_imp_50 = ggplot(plot_50, aes(x=reorder(Feature, -AbsCoef, FUN=mean), y=AbsCoef, fill=Model)) +
  geom_col(position='dodge', width=0.7, alpha=0.8) +
  scale_fill_manual(values=c('166'='#E41A1C', '82'='#377EB8', 'EDRN'='#4DAF4A')) +
  coord_flip() +
  labs(
    title='Consensus Features: Coefficients (Top 50 Analysis)',
    subtitle='Features appearing in all 3 models',
    x='Feature',
    y='|Coefficient|',
    fill='Model'
  ) +
  theme_bw() +
  theme(axis.text.y=element_text(size=9), plot.title=element_text(face='bold', size=12))

ggsave(file.path(plot_dir, 'feature_importance_top50_3models.pdf'), p_imp_50, width=13, height=9)
cat('Saved: feature_importance_top50_3models.pdf\n')

cat('\n\n=== ALL PLOTS REGENERATED ===\n')
cat('Top 20:\n')
cat('  ✓ venn_top20_3models.pdf\n')
cat('  ✓ feature_importance_top20_3models.pdf\n')
cat('\nTop 50:\n')
cat('  ✓ venn_top50_3models.pdf\n')
cat('  ✓ feature_importance_top50_3models.pdf\n')
