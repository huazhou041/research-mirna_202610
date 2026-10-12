library(tidyverse)
library(glmnet)
library(pROC)
library(ggplot2)
library(ggvenn)

setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

cat('=== UNIFIED GLM ANALYSIS: 3 COHORTS ===\n\n')

# Load new data
d_166 = readRDS('data/model_data/nyu_166_miRNA.RDS')
d_82 = readRDS('data/model_data/nyu_82_miRNA.RDS')

# Load previous model (best: v5 with 4 features)
prev_model = readRDS('prev/models/EV_RNA_miRNA_glm_v5.rds')
prev_features = rownames(prev_model$glmnet.fit$beta)[-1]  # Exclude intercept

set.seed(123)
d_166_clean = d_166[!is.na(d_166$recurrence), ]
d_82_clean = d_82[!is.na(d_82$recurrence), ]

cat('=== 3 COHORT OVERVIEW ===\n\n')
cat('1. EDRN Model (EV_RNA_miRNA_glm_v5 - BEST):\n')
cat('   - 4 features trained on EDRN data\n')
cat('   - Applied to nyu_166: AUC 0.614\n')
cat('   - Applied to nyu_82: AUC 0.827\n\n')

cat('2. nyu_166_miRNA (Best: RF v20):\n')
cat('   - 20 features (MAD-selected)\n')
cat('   - RF AUC: 0.789\n')
cat('   - GLM AUC: 0.757\n\n')

cat('3. nyu_82_miRNA (Best: RF v50):\n')
cat('   - 50 features (MAD-selected)\n')
cat('   - RF AUC: 0.974\n')
cat('   - GLM AUC: 0.946\n\n')

# ===== Extract features from each dataset =====
cat('Extracting features...\n\n')

# 166 - top 20 MAD
gene_cols_166 = setdiff(colnames(d_166_clean), 'recurrence')
mad_vals_166 = apply(d_166_clean[, gene_cols_166], 2, mad, na.rm=TRUE)
top_20_genes_166 = gene_cols_166[order(mad_vals_166, decreasing=TRUE)[1:20]]

# 82 - top 50 MAD
gene_cols_82 = setdiff(colnames(d_82_clean), 'recurrence')
mad_vals_82 = apply(d_82_clean[, gene_cols_82], 2, mad, na.rm=TRUE)
top_50_genes_82 = gene_cols_82[order(mad_vals_82, decreasing=TRUE)[1:50]]

# Previous model - 10 features
prev_feat = prev_features

# Find overlaps
overlap_all_3 = intersect(intersect(top_20_genes_166, top_50_genes_82), prev_feat)
overlap_166_82 = intersect(top_20_genes_166, top_50_genes_82)
overlap_166_prev = intersect(top_20_genes_166, prev_feat)
overlap_82_prev = intersect(top_50_genes_82, prev_feat)

cat('Feature overlaps:\n')
cat('  166 ∩ 82:', length(overlap_166_82), '\n')
cat('  166 ∩ Prev:', length(overlap_166_prev), '\n')
cat('  82 ∩ Prev:', length(overlap_82_prev), '\n')
cat('  All 3:', length(overlap_all_3), '\n\n')

if (length(overlap_all_3) > 0) {
  cat('Features in all 3:\n')
  for (f in overlap_all_3) {
    cat('  ', f, '\n')
  }
  cat('\n')
}

# ===== Extract GLM coefficients =====
cat('Computing GLM coefficients for each cohort...\n\n')

extract_glm_coefs = function(data, feature_set, name) {
  X = as.matrix(data[, feature_set])
  y = data$recurrence

  m = cv.glmnet(X, y, family='binomial', nfolds=5, alpha=0.5)
  coef_m = coef(m, s='lambda.min')
  coef_vec = as.vector(coef_m[-1])
  names(coef_vec) = feature_set

  # AUC on full data
  pred = predict(m, newx=X, s='lambda.min', type='response')
  auc = roc(y, as.vector(pred), quiet=TRUE)$auc

  list(coefs = coef_vec, auc = auc, features = names(coef_vec))
}

glm_166 = extract_glm_coefs(d_166_clean, top_20_genes_166, '166')
glm_82 = extract_glm_coefs(d_82_clean, top_50_genes_82, '82')

# For previous, extract coefficients
coef_prev = coef(prev_model, s='lambda.min')
glm_prev_coefs = as.vector(coef_prev[-1])
names(glm_prev_coefs) = prev_features

cat('GLM AUC on full datasets:\n')
cat('  166 miRNA:', round(glm_166$auc, 3), '\n')
cat('  82 miRNA: ', round(glm_82$auc, 3), '\n')
cat('  Previous: (0.614 on 166, 0.827 on 82 from old data)\n\n')

# ===== Find consensus features =====
all_features = unique(c(names(glm_166$coefs), names(glm_82$coefs), names(glm_prev_coefs)))

feature_consensus = data.frame(
  Feature = all_features,
  AbsCoef_166 = NA,
  AbsCoef_82 = NA,
  AbsCoef_EDRN = NA,
  Count_Models = NA
)

for (i in 1:nrow(feature_consensus)) {
  f = feature_consensus$Feature[i]
  feature_consensus$AbsCoef_166[i] = abs(glm_166$coefs[f])
  feature_consensus$AbsCoef_82[i] = abs(glm_82$coefs[f])
  feature_consensus$AbsCoef_EDRN[i] = abs(glm_prev_coefs[f])
  feature_consensus$Count_Models[i] = sum(!is.na(c(
    if (f %in% names(glm_166$coefs)) 1 else NA,
    if (f %in% names(glm_82$coefs)) 1 else NA,
    if (f %in% names(glm_prev_coefs)) 1 else NA
  )))
}

feature_consensus = feature_consensus %>%
  mutate(Mean_AbsCoef = rowMeans(.[2:4], na.rm=TRUE)) %>%
  arrange(desc(Count_Models), desc(Mean_AbsCoef))

cat('Top 20 Consensus Features (appear in 2+ models):\n\n')
top_consensus = feature_consensus %>% filter(Count_Models >= 2) %>% head(20)
print(top_consensus)

write.csv(feature_consensus, 'results/glm_consensus_3cohorts.csv', row.names=FALSE)
cat('\n\nSaved: glm_consensus_3cohorts.csv\n')

# ===== Improved Venn Diagram =====
cat('\nCreating 3-way Venn diagram...\n')

venn_list = list(
  '166 miRNA\n(top 20)' = top_20_genes_166,
  '82 miRNA\n(top 50)' = top_50_genes_82,
  'EDRN\n(v5)' = prev_feat
)

p_venn = ggvenn(venn_list,
                 fill_color = c('#E41A1C', '#377EB8', '#4DAF4A'),
                 stroke_size = 2.5,
                 set_name_size = 6,
                 text_size = 5) +
  labs(title = 'Feature Overlap: All 3 Cohorts (nyu_166, nyu_82, EDRN)',
       subtitle = 'Top MAD-selected features from each cohort') +
  theme_bw() +
  theme(plot.title = element_text(face='bold', size=14),
        plot.subtitle = element_text(size=11))

ggsave('results/plots/venn_3cohorts_improved.pdf', p_venn, width=11, height=10)
cat('Saved: venn_3cohorts_improved.pdf\n')

# ===== Feature Importance Plot =====
cat('Creating feature importance plot...\n\n')

plot_data = feature_consensus %>%
  head(25) %>%
  pivot_longer(cols=c(AbsCoef_166, AbsCoef_82, AbsCoef_EDRN),
               names_to='Cohort', values_to='AbsCoef') %>%
  mutate(Cohort = factor(Cohort, levels=c('AbsCoef_166', 'AbsCoef_82', 'AbsCoef_EDRN'),
                         labels=c('nyu_166', 'nyu_82', 'EDRN')))

p_imp = ggplot(plot_data, aes(x=reorder(Feature, -AbsCoef, FUN=mean), y=AbsCoef, fill=Cohort)) +
  geom_col(position='dodge', width=0.7, alpha=0.8) +
  scale_fill_manual(values=c('nyu_166'='#E41A1C', 'nyu_82'='#377EB8', 'EDRN'='#4DAF4A')) +
  coord_flip() +
  labs(
    title='GLM Feature Coefficients: 3 Cohorts (nyu_166, nyu_82, EDRN)',
    x='Feature',
    y='|Coefficient|',
    fill='Cohort'
  ) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=9),
    axis.text.x=element_text(size=10),
    plot.title=element_text(face='bold', size=12),
    legend.position='right'
  )

ggsave('results/plots/glm_3cohorts_importance_improved.pdf', p_imp, width=13, height=10)
cat('Saved: glm_3cohorts_importance_improved.pdf\n')

cat('\n\n=== RECOMMENDATION ===\n')
consensus_count = nrow(feature_consensus %>% filter(Count_Models >= 2))
cat('Consensus features (2+ cohorts):', consensus_count, '\n')
cat('Consensus features (all 3): ', nrow(feature_consensus %>% filter(Count_Models == 3)), '\n')
cat('\nUse top 20 consensus features for unified GLM model\n')
