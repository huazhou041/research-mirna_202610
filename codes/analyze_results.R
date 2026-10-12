rm(list=ls())
library(ggplot2)
library(pROC)
setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

resdir  = 'results'
outdir  = 'results/plots'
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

varlist   = c(5, 10, 20, 50, 100, 200, 500, 1000)
datasets  = list(
  nyu_166 = c('Clinical','miRNA','Both'),
  nyu_82  = c('miRNA')
)

# ============================================================
# 1. AUC SUMMARY TABLE
# ============================================================
rows = list()
for (dname in names(datasets)) {
  for (subset in datasets[[dname]]) {
    row = data.frame(dataset=dname, subset=subset,
                     stringsAsFactors=FALSE)
    for (v in varlist) {
      f = file.path(resdir, 'raw_pred', sprintf('%s_%s_glm_v%d_cvSummary.csv', dname, subset, v))
      if (file.exists(f)) {
        tmp  = read.csv(f, stringsAsFactors=FALSE)
        tmp  = tmp[!is.na(tmp$AUC), ]
        tauc = mean(tmp$AUC)
        tsd  = sd(tmp$AUC)
        tl   = max(0,   tauc - 1.96*tsd/sqrt(nrow(tmp)))
        th   = min(1.0, tauc + 1.96*tsd/sqrt(nrow(tmp)))
        row[[paste0('v', v)]] = sprintf('%.2f (%.2f~%.2f)', tauc, tl, th)
      } else {
        row[[paste0('v', v)]] = NA
      }
    }
    rows[[length(rows)+1]] = row
  }
}
auc_table = do.call(rbind, rows)
write.csv(auc_table, file.path(resdir, 'auc_summary.csv'), row.names=FALSE, quote=FALSE)
cat('=== AUC Summary ===\n')
print(auc_table)

# ============================================================
# 2. ROC CURVES — best varcut per dataset/subset
# ============================================================
roc_theme = theme(
  panel.border      = element_rect(fill=NA, color='black', linewidth=1),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  panel.background  = element_blank(),
  axis.text         = element_text(size=14),
  axis.title        = element_text(size=14),
  legend.position   = c(.6, .1),
  legend.title      = element_blank(),
  legend.text       = element_text(size=12)
)

plot_roc = function(dname, subset, v) {
  prob_f = file.path(resdir, 'raw_pred', sprintf('%s_%s_glm_v%d_pred_prob.csv', dname, subset, v))
  rds_f  = file.path(resdir, 'raw_pred', sprintf('%s_%s_glm_v%d_cvSummary.csv', dname, subset, v))
  if (!file.exists(prob_f)) return(NULL)

  dat_f  = file.path('data/model_data', sprintf('%s_%s.RDS', dname, subset))
  d      = readRDS(dat_f)
  probs  = read.csv(prob_f, stringsAsFactors=FALSE)

  # use bootid==1 as representative fold
  p1     = probs[probs$bootid == 1, ]
  truth  = d$recurrence[match(p1$subj, rownames(d))]
  roc_obj = roc(truth, p1$prob, smooth=FALSE)

  # mean AUC with CI across all bootstraps
  tmp    = read.csv(rds_f, stringsAsFactors=FALSE)
  tmp    = tmp[!is.na(tmp$AUC), ]
  tauc   = mean(tmp$AUC); tsd = sd(tmp$AUC)
  tl     = max(0,   tauc - 1.96*tsd/sqrt(nrow(tmp)))
  th     = min(1.0, tauc + 1.96*tsd/sqrt(nrow(tmp)))
  lbl    = sprintf('AUC = %.2f (%.2f~%.2f)', tauc, tl, th)

  ciobj  = tryCatch(ci.se(roc_obj, specificities=seq(0,1,l=25)), error=function(e) NULL)
  p = ggroc(setNames(list(roc_obj), lbl), color='blue', linewidth=1.1) +
    geom_abline(slope=1, intercept=1, linetype='dashed', alpha=0.7, color='red') +
    labs(x='Specificity', y='Sensitivity',
         title=sprintf('%s | %s | top-%d miRNAs', dname, subset, v)) +
    roc_theme + coord_equal()
  if (!is.null(ciobj)) {
    dat_ci = data.frame(x=as.numeric(rownames(ciobj)),
                        lower=ciobj[,1], upper=ciobj[,3])
    p = p + geom_ribbon(data=dat_ci, aes(x=x, ymin=lower, ymax=upper),
                        fill='grey', alpha=0.2, inherit.aes=FALSE)
  }
  p
}

# find best varcut per dataset/subset (highest mean AUC)
best_v = list()
for (dname in names(datasets)) {
  for (subset in datasets[[dname]]) {
    aucs = sapply(varlist, function(v) {
      f = file.path(resdir, 'raw_pred', sprintf('%s_%s_glm_v%d_cvSummary.csv', dname, subset, v))
      if (!file.exists(f)) return(NA)
      tmp = read.csv(f, stringsAsFactors=FALSE)
      mean(tmp$AUC[!is.na(tmp$AUC)])
    })
    best_v[[paste(dname, subset)]] = varlist[which.max(aucs)]
  }
}

cat('\n=== Best varcut per model ===\n')
print(best_v)

# plot ROC for every varcut, save best separately
for (dname in names(datasets)) {
  for (subset in datasets[[dname]]) {
    plots = list()
    for (v in varlist) {
      p = plot_roc(dname, subset, v)
      if (!is.null(p)) plots[[as.character(v)]] = p
    }
    # save individual ROCs as one PDF (8 panels)
    pdf(file.path(outdir, sprintf('roc_%s_%s_allv.pdf', dname, subset)),
        width=16, height=8)
    n = length(plots)
    par_rows = ceiling(n/4)
    gridExtra::grid.arrange(grobs=plots, nrow=par_rows)
    dev.off()

    # save best varcut ROC separately
    bv = best_v[[paste(dname, subset)]]
    p_best = plot_roc(dname, subset, bv)
    if (!is.null(p_best)) {
      ggsave(file.path(outdir, sprintf('roc_%s_%s_best_v%d.pdf', dname, subset, bv)),
             p_best, width=6, height=6)
    }
  }
}

# ============================================================
# 3. BETA COEFFICIENTS — best model per dataset/subset
# ============================================================
library(glmnet)
plot_beta = function(dname, subset, v) {
  f = file.path(resdir, sprintf('%s_%s_glm_v%d.rds', dname, subset, v))
  if (!file.exists(f)) return(NULL)
  m    = readRDS(f)
  beta = m$glmnet.fit$beta[, which.min(abs(m$lambda - m$lambda.min))]
  beta = beta[abs(beta) > 0]
  if (length(beta) == 0) return(NULL)
  df = data.frame(gene=names(beta), beta=as.numeric(beta), stringsAsFactors=FALSE)
  df$gene = gsub('Tissue_', '', df$gene)
  df = df[order(df$beta), ]
  df$gene = factor(df$gene, levels=df$gene)
  df$dir  = ifelse(df$beta > 0, 'pos', 'neg')
  ggplot(df, aes(y=gene, x=beta, fill=dir)) +
    geom_bar(stat='identity') +
    scale_fill_manual(values=c(pos='red', neg='blue')) +
    labs(title=sprintf('%s | %s | top-%d', dname, subset, v),
         x='Coefficient', y=NULL) +
    theme_bw() + theme(legend.position='none',
                       axis.text.y=element_text(size=max(4, 10-length(beta)/20)))
}

for (dname in names(datasets)) {
  for (subset in datasets[[dname]]) {
    bv = best_v[[paste(dname, subset)]]
    p  = plot_beta(dname, subset, bv)
    if (!is.null(p)) {
      ht = max(4, min(12, 0.2 * length(coef(readRDS(
        file.path(resdir, sprintf('%s_%s_glm_v%d.rds', dname, subset, bv))
      ), s='lambda.min')[,1])))
      ggsave(file.path(outdir, sprintf('beta_%s_%s_best_v%d.pdf', dname, subset, bv)),
             p, width=6, height=ht)
    }
  }
}

cat('\nDone. Output in results/plots/\n')
