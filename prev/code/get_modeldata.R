rm(list=ls())
library(readxl)
setwd('/Volumes/hua_mac/research/aris/harvey/erdn/erdn_202503/')

outdir='data/model_data/'
if (!dir.exists(outdir))
  dir.create(outdir,recursive = T)

expr=data.frame(read_excel('data/Norm.miRNA.xlsx',1,skip=1),stringsAsFactors = F,check.names=F)
genes=expr[,1]
genes=gsub('[+]','_',genes)
#genes=gsub('hsa-','',genes)
#genes=gsub('[-|+]','_',genes)
rownames(expr)=genes
expr=data.frame(t(expr[,-c(1:5)]),check.names=F)
meta=data.frame(read_excel('data/EDRN_Plasma_Mapping_checked.xlsx',1),stringsAsFactors = F,check.names=F)
rownames(meta)=meta[,1]
meta$height=-9
for (i in 1:nrow(meta)) {
  if (meta$height_feet[i]!='NA') {
    if (meta$height_inch[i]=='NA')
      inch=0
    else
      inch=as.numeric(meta$height_inch[i])
    meta$height[i]=as.numeric(meta$height_feet[i])*12+inch
  }}

subjs=meta[match(rownames(expr),meta[,1]),c(5,15,16,17,18,20,23,24,25,26,27)]
#1st 
meta_names=c('Sex','Race','Smoking','Pack_years','Age','Weight','Smoke_mariju','Smoke_mariju_frq','Tumor_size','Height')
names(subjs)=c('recurrence',paste0('Tissue_',meta_names))
subjs[,1]=ifelse(subjs[,1]=='Lungcancer', 1, 0)

subjs[subjs[,2]=='NA',2]=-9
subjs[subjs[,2]=='Female',2]=0
subjs[subjs[,2]=='Male',2]=1

#subjs[subjs[,3]=='white']=1
subjs[is.na(subjs[,3]) | subjs[,3]=='NA',3]='Other'
subjs[subjs[,3] %in% c('Unknown'),3]='Other'
subjs[subjs[,3]=='Hispanic or Latino',3]='Hispanic'
subjs[subjs[,3]=='Black or African American',3]='Black'
#subjs[,3]=factor(subjs[,3],levels=c('Other','Hispanic','Black','Asian','White'))
subjs[subjs[,3]=='Other',3]=1
subjs[subjs[,3]=='Hispanic',3]=2
subjs[subjs[,3]=='Black',3]=3
subjs[subjs[,3]=='Asian',3]=4
subjs[subjs[,3]=='White',3]=5
subjs[,3]=as.numeric(subjs[,3])

subjs[subjs[,4]=='NA',4]=-9
subjs[subjs[,4]=='Former',4]=0
subjs[subjs[,4]=='Current',4]=1

subjs[is.na(subjs[,5]) | subjs[,5]=='NA',5]=-9
subjs[subjs[,5]=='less than 1ppd x 30 yrs',5]=0
subjs[subjs[,5]=='30-15yrs, 5-5yrs, 1.5-3yrs',5]=1.5
subjs[subjs[,5]=='>30',5]=30
subjs[,5]=as.numeric(subjs[,5])


subjs[subjs[,6]=='NA',6]=-9
subjs[,6]=as.numeric(subjs[,6])

subjs[subjs[,7]=='NA',7]=-9
subjs[,7]=as.numeric(subjs[,7])

subjs[subjs[,8]=='NA',8]=-9
subjs[subjs[,8]=='Former',8]=0
subjs[subjs[,8]=='Current',8]=1


subjs[subjs[,9]=='NA',9]='Unknown'
#subjs[,9]=factor(subjs[,9],levels=c('Unknown','Occasionally','Monthly','Weekly','Daily'))
subjs[subjs[,9]=='Unknown',9]=1
subjs[subjs[,9]=='Occasionally',9]=2
subjs[subjs[,9]=='Monthly',9]=3
subjs[subjs[,9]=='Weekly',9]=4
subjs[subjs[,9]=='Daily',9]=5
subjs[,9]=as.numeric(subjs[,9])

subjs[,10]=gsub('cm','',subjs[,10])
subjs[,10]=as.numeric(subjs[,10])
subjs[is.na(subjs[,10]),10]=-9

for (i in 2:ncol(subjs)) 
  subjs[,i]=as.numeric(subjs[,i])
data=cbind.data.frame(subjs,expr)




data=data[,-c(7:11)]
d1=data[,c(1:6)]

saveRDS(d1,paste0(outdir,'EV_RNA_Clinical.RDS'))
d1=data[,c(1,7:ncol(data))]
saveRDS(d1,paste0(outdir,'EV_RNA_miRNA.RDS'))
d1=data
saveRDS(d1,paste0(outdir,'EV_RNA_Both.RDS'))



for (i in 2:11) {
  t1= t.test(data[data[,1]==1,i],data[data[,1]==0,i])
  cat(gsub('Tissue_','',names(data)[i]),' : ',t1$p.value,'\n')
}

# summ=summ[apply(summ=='NaN',1,sum)<ncol(summ),]
# write.csv(summ,file=paste0('data/model_data/counts_summary.csv'),row.names=F)

