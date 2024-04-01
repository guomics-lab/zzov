library(xlsx)
readin<-read.csv("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.csv",header = T,row.names = 1)
# min(readin,na.rm = T)
na_to<-2389.27
readin$missing<-apply(is.na(readin),1,sum)
readin$na<-readin$missing/(ncol(readin)-1)
readin1<-readin[!readin$na>0.7,1:923]
readin1[is.na(readin1)]<-na_to
log_na_to<-log2(na_to)
readin1<-data.frame(t(readin1))
readin1<-log2(readin1)
inf1<-read.xlsx("20221228genoeme_ZZOVA_todo.xlsx",sheetIndex = 1)

sample_inf<-read.xlsx("20210728ZZOV_patient_sample_info.xlsx",sheetIndex = 4)
readin1$patiend_id<-sample_inf$bcr_patient_barcode[match(row.names(readin1),sample_inf$Sample_name)]
inf5<-read.xlsx("20221228genoeme_ZZOVA_todo.xlsx",sheetIndex = 5)

df5<-readin1[readin1$patiend_id%in%inf5$group5A,]
df5$type<-inf5$NA.[match(df5$patiend_id,inf5$group5A)]
# write.csv(df4,"zzov_group4A_4B_allprot.csv",row.names = T)
readP<-df5[grepl("group5A",df5$type),1:(ncol(df5)-2)]
readR<-df5[grepl("group5B",df5$type),1:(ncol(df5)-2)]

FC<-pvalue<-vector()
for(i in 1:ncol(readP)){
  
  if(sum(!readP[,i]==log_na_to)>3&sum(!readR[,i]==log_na_to)>3){
    FC[i]<-mean(2^readP[,i],na.rm=T)/mean(2^readR[,i],na.rm = T)
    pvalue[i]<-t.test(readP[,i],readR[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
  
}
pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-cbind(FC,pvalueAd,pvalue)
row.names(FC_PAd)<-colnames(readP)[1:ncol(readP)]
write.csv(FC_PAd,"20230424_zzov_923sample_t_test_allprot_FC_PAd_group5A_vs_group5B.csv")
test_list<-read.xlsx("20230424fig7Binput.xlsx",sheetIndex = 2)
FC_PAd<-data.frame(FC_PAd)
FC_PAd$prot<-row.names(FC_PAd)
FC_PAd$prot<-sapply(strsplit(FC_PAd$prot,"_"),function(e){e[1]})
test_list1<-test_list[test_list$label==1,]
test_lab<-FC_PAd[FC_PAd$prot%in%test_list1$X33prots,]
color_label<-FC_PAd[FC_PAd$prot%in%test_list$X33prots,]
hist(pvalue)
hist(log2(FC),breaks = 1000)
pdf("20230424_zzovc_42tiss_sen_vs_res_t_test_volcan_0.8min_70miss_1_5fc_1.pdf",width =6,height = 6)
plot(log2(FC),-log10(pvalue),
     pch=19, col="gray",
     xlim=c(-5,5),
     main="sen vs res"
     ,xlab="log2 FC",ylab="-log10 pvalue",
     cex.axis = 1.5,cex.lab=1.5)
up <- which(pvalue< 0.05 & FC >1.5)
Uppoint1<-which(-log10(pvalue)>1.30103&FC > 1.5)
down <- which(pvalue < 0.05 & FC <  (1/1.5))
Downpoint1<-which(-log10(pvalue)>1.30103&FC < (1/1.5))
points(log2(FC[up]), -log10(pvalue[up]), col=1, bg = "red", pch=21, cex=1)
# text(log2(FC[Uppoint1]), -log10(pvalue[Uppoint1]),labels = colnames(readP)[Uppoint1], col=1, bg = "grey", cex=0.8,font = 3)
points(log2(FC[down]), -log10(pvalue[down]), col = "blue", pch=19,cex=1)
points(log2(color_label$FC), -log10(color_label$pvalue), col = "orange", pch=19,cex=1)
text(log2(test_lab$FC), -log10(test_lab$pvalue),labels =row.names(test_lab), col=1, bg = "grey", cex=0.5,font = 3)
abline(h=-log10(0.05),v=c(-log2(1.5),log2(1.5)),lty=2,lwd=1)
dev.off()

protUP<-FC_PAd[up,]
protdown<-FC_PAd[down,]
ipainput<-rbind(protUP,protdown)
write.csv(ipainput,"zzovc_42tiss_sen_vs_res_t_test_diffprot_FC_PAd_p2p_50miss_1_5fc.csv")