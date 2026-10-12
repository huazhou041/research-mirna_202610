rm(list=ls())
setwd('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/')


datanames=c( 'Clinical', 'miRNA', 'Both')
varlist=c(5,10,20,50,100,200,500,1000)



res=data.frame(matrix(NaN,length(datanames),length(varlist)+2))
names(res)=c('Model','dataname',paste0('var_',varlist))
cnt=0
res_file=paste('auc_summary_evRNA.csv',sep='')

        for (i in 1:length(datanames)) {
            cnt=cnt+1
            res[cnt,1]='GLM'
            res[cnt,2]=datanames[i]
            for (k in 1:length(varlist)) {
              curname=paste('./ML_res/EV_RNA_',datanames[i],'_glm_v',varlist[k],'_cvSummary.csv',
                            sep='')
              if (file.exists(curname)) {
                temp=read.csv(curname,header=T,stringsAsFactors = F)
                temp=temp[order(temp$AUC,decreasing = T),]
                temp=temp[!is.na(temp$AUC),]
                auc_sd=sd(temp$AUC)
                ###doing cis
                tauc=mean(temp$AUC)
                tauc_l=max(0,tauc-1.96*auc_sd/sqrt(10))
                tauc_h=min(1.0,tauc+1.96*auc_sd/sqrt(10))
                res[cnt,2+k]=sprintf('%.2f (%.2f~%.2f)',c(round(tauc,2)), round(tauc_l,2),round(tauc_h,2))
              }
            }
          }


  write.csv(res,res_file,row.names=F)






