rm(list=ls())
library(ggplot2)
library(readxl)
library(pROC)
setwd('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/')

outpath='result/'
if (!dir.exists(outpath))
  dir.create(outpath,recursive = T)
##doing coefficient of best model both
  curname='/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/ML_res/EV_RNA_Both_glm_v50.rds'
  temp=readRDS(curname)
  ind=which.min(abs(temp$lambda-temp$lambda.min))
  beta=temp$glmnet.fit$beta[,ind]
  if (sum(abs(beta)>0)==0){
    ind=round(length(temp$lambda)/2)
    beta=temp$glmnet.fit$beta[,ind]
  }
 
  beta=beta[abs(beta)>0]
  beta=data.frame(gene=names(beta),beta=beta,stringsAsFactors = F)
  beta$direction='neg'
  beta$direction[beta[,2]>0]='pos'
  beta=beta[order(beta[,2]),]
  beta$gene=gsub('Tissue_','',beta$gene)
  beta$gene=factor(beta$gene, levels = beta$gene)
 p1= ggplot(beta,aes(y=gene,x=beta,fill=direction))+geom_bar(stat="identity")+scale_fill_manual(values=c(pos='red',neg='blue'))

  ggsave(paste0(outpath,'beta_v50_both.pdf'),p1,width=6,height=7)
  
  
  
  ##doing coefficient of best model miRNA only
  curname='/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/ML_res/EV_RNA_miRNA_glm_v50.rds'
  temp=readRDS(curname)
  ind=which.min(abs(temp$lambda-temp$lambda.min))
  beta=temp$glmnet.fit$beta[,ind]
  if (sum(abs(beta)>0)==0){
    ind=round(length(temp$lambda)/2)
    beta=temp$glmnet.fit$beta[,ind]
  }
  
  beta=beta[abs(beta)>0]
  beta=data.frame(gene=names(beta),beta=beta,stringsAsFactors = F)
  beta$direction='neg'
  beta$direction[beta[,2]>0]='pos'
  beta=beta[order(beta[,2]),]
  beta$gene=gsub('Tissue_','',beta$gene)
  beta$gene=factor(beta$gene, levels = beta$gene)
  p1= ggplot(beta,aes(y=gene,x=beta,fill=direction))+geom_bar(stat="identity")+scale_fill_manual(values=c(pos='red',neg='blue'))
  write.csv(beta,paste0(outpath,'beta_v50_miRNA.csv'),row.names=F,quote=F)
  ggsave(paste0(outpath,'beta_v50_miRNA.pdf'),p1,width=6,height=7)
  
  
  ###doing ttest
  data=readRDS('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/data/model_data/EV_RNA_miRNA.RDS')
  ##
  res=data.frame(matrix(NaN,ncol(data)-1,4))
  names(res)=c('gene','lfc','pval','fdr')

  for (i in 2:ncol(data)) {
    res[i-1,1]=colnames(data)[i]
    d1=data[data[,1]==1,i]
    d2=data[data[,1]==0,i]
    t1=t.test(d1,d2)
    res[i-1,2]=log2(mean(d1)/mean(d2))
    res[i-1,3]=t1$p.value
  }
  res[,4]=p.adjust(res[,3],method='fdr')
  
  write.csv(res,paste0(outpath,'ttest_evrna_lfc.csv'),row.names=F,quote=F)
  
  ###doing roc both
  d1=readRDS('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/data/model_data/EV_RNA_Both.RDS')
  temp=read.csv('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/ML_res/EV_RNA_Both_glm_v50_pred_prob.csv',header=T,stringsAsFactors = F,check.names=F)
  temp=temp[temp$bootid==1,]
  roc.list=roc(d1$recurrence[match(temp$subj,rownames(d1))],temp$prob,smooth=F,thresholds=30)
  ciobj <- ci.se(roc.list, specificities=seq(0, 1, l=25))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                      upper = ciobj[, 3])
  auc_label='AUC = 0.95 (0.92~0.97)'
  roclist=list()
  roclist[[auc_label]]=roc.list
  p1=ggroc(roclist,color='blue',size=1.1,aes=c('linetype')) +
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "red") + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "grey", alpha= 0.2, inherit.aes = F) +
    labs(x='Specificity',y='Sensitivity',linetype='auc')+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          panel.background = element_blank(),
          legend.position =c(.6,.1),
          legend.title = element_blank(),
          legend.text = element_text(size=20))+
    coord_equal()
  filename=paste(outpath,'/roc_50features_and_clinical.pdf',collapse='',sep='')
  ggsave(filename = filename,p1,width=8,height=6)
  
  
  
  
  ###doing roc only gene
  d1=readRDS('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/data/model_data/EV_RNA_miRNA.RDS')
  temp=read.csv('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/ML_res/EV_RNA_miRNA_glm_v50_pred_prob.csv',header=T,stringsAsFactors = F,check.names=F)
  temp=temp[temp$bootid==1,]
  roc.list=roc(d1$recurrence[match(temp$subj,rownames(d1))],temp$prob,smooth=F,thresholds=30)
  ciobj <- ci.se(roc.list, specificities=seq(0, 1, l=25))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  auc_label='AUC = 0.91 (0.88~0.94)'
  roclist=list()
  roclist[[auc_label]]=roc.list
  p1=ggroc(roclist,color='blue',size=1.1,aes=c('linetype')) +
    geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "red") + 
    geom_ribbon(data = dat.ci, aes(x = x, ymin = lower, ymax = upper), fill = "grey", alpha= 0.2, inherit.aes = F) +
    labs(x='Specificity',y='Sensitivity',linetype='auc')+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.0, linetype="solid"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          panel.background = element_blank(),
          legend.position =c(.6,.1),
          legend.title = element_blank(),
          legend.text = element_text(size=20))+
    coord_equal()
  filename=paste(outpath,'/roc_50features.pdf',collapse='',sep='')
  ggsave(filename = filename,p1,width=8,height=6)