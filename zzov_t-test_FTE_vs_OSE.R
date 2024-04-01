library(openxlsx)
readin<-read.xlsx("cellreports_FTEvsOSE_proteinMatrix.xlsx",sheet = 1)
readin$prot<-paste(readin$Uniprot,readin$Gene.names,sep = "_")
df<-data.frame(t(readin[,10:18]))
names(df)<-readin$prot

readP<-df[grepl("FTE",row.names(df)),]
readR<-df[grepl("OSE",row.names(df)),]

FC<-pvalue<-vector()
for(i in 1:ncol(readP)){
  
  if(sum(!is.na(readP[,i]))>1&sum(!is.na(readR[,i]))>1){
    FC[i]<-mean(2^readP[,i],na.rm=T)/mean(2^readR[,i],na.rm = T)
    tryCatch(pvalue[i]<-t.test(readP[,i],readR[,i], paired = F,  var.equal = F)$p.value,error=function(e) {skip_to_next <<- TRUE})
  }else{
    FC[i]<-pvalue[i]<-NA
  }
  
}
pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd,pvalue)
row.names(FC_PAd)<-colnames(readP)[1:ncol(readP)]
write.csv(FC_PAd,"20240130_FTEvsOSE_t_test_allprot_FC_PAd.csv")