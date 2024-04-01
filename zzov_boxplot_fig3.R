library(xlsx)
input1<-read.xlsx("Supplementary Table 4_histotype specific DEPs20230417(1).xlsx",sheetIndex = 2)
input1<-input1[-1,]
library(xlsx)
samp_inf<-read.xlsx("20221105ZZOV_patient_sample_info(1)(1).xlsx",sheetIndex = 2)
# prot_list_3b<-read.xlsx("20230501ref_exampleProt.xlsx",sheetIndex = 2)
prot_mat<-read.table("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.txt",sep = "\t",header = T)
samp<-read.xlsx("input3_group.xlsx",sheetIndex = 2)
samp$samp<-substr(samp$Sample.name,1,nchar(samp$Sample.name)-1)
prot_mat$prot<-sapply(strsplit(row.names(prot_mat),"_"),function(e){e[1]})
df_prot<-prot_mat[prot_mat$prot%in%input1$UniprotID,-ncol(prot_mat)]

df_prot<-data.frame(t(df_prot))

df_prot$pat_id<-samp$Bcr.patient.barcode[match(row.names(df_prot),samp$samp)]
df_prot$Histology8<-samp_inf$Histology8[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$Group<-samp_inf$Group[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$RFS<-samp_inf$Recurrence.free.survival.RFS..month[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$dis<-samp_inf$Recurrent.disease[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$age<-samp_inf$Age_at_diagnosis[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$type<-NA
df_prot$type[df_prot$Histology8=="NM"]<-"Normal"
df_prot$type[df_prot$Histology8=="HS"&df_prot$Group=="Primary carcinoma"]<-"HGSOC"
df_prot$type[df_prot$Histology8=="LS"&df_prot$Group=="Primary carcinoma"]<-"LGSOC"
df_prot$type[df_prot$Histology8=="CC"&df_prot$Group=="Primary carcinoma"]<-"CCOC"
df_prot$type[df_prot$Histology8=="MC"&df_prot$Group=="Primary carcinoma"]<-"MCOC"
df_prot$type[df_prot$Histology8=="EC"&df_prot$Group=="Primary carcinoma"]<-"EMOC"

df_prot1<-df_prot[!is.na(df_prot$type),]
df_prot1$type<-factor(df_prot1$type,levels = c("Normal","HGSOC","LGSOC","CCOC","MCOC","EMOC"))
df_prot1<-df_prot1[order(df_prot1$type),]
mean_temp<-data.frame(matrix(NA,nrow =length(unique(df_prot1$type)),ncol = (ncol(df_prot1)-4)))
row.names(mean_temp)<-unique(df_prot1$type)
colnames(mean_temp)<-colnames(df_prot1)[1:(ncol(df_prot1)-4)]
for (i in 1:nrow(mean_temp)) {
  indsel<-df_prot1[df_prot1$type==row.names(mean_temp)[i],]
  mean_temp[i,]<-apply(indsel[,1:(ncol(indsel)-4)],2,mean,na.rm=T)
}
write.csv(mean_temp,"20230529_zzov_table4_prot_intensity_mean_group.csv",row.names = T)


ccoc<-df_prot[grepl("CCOC",df_prot$type),]
ccoc<-ccoc[!grepl("NOT available",ccoc$RFS),]
write.csv(ccoc,"test.csv",row.names = T)
ccoc<-read.csv("test.csv",header = T,row.names = 1)
int <- function(x){
  #inverse normal transformation
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

## rank-based INT
# Age residuals + standardization
res_stan <- function(column){
  lm.fit <- lm(formula = column~age,data=ccoc)
  res <- lm.fit$residuals
  res <- int(res)
}

pro_res <- data.frame(sapply(ccoc[,1:4327], res_stan))
# pro_res <-data.frame(pro_res)
# row.names(pro_res) = readin$patient_barcode
library("survival")
library("survminer")

pro_res<-cbind(pro_res,ccoc[,4331:4332])
# yes=2 no=1
pro_res$dis<-factor(pro_res$dis, labels = c(0,1)) %>% as.numeric()
pro_res$RFS<-as.numeric(pro_res$RFS)


pvalue_cox<-exp<-ecp_low<-exp_up<-vector()
for (i in 1:(ncol(pro_res)-2)) {
  pvalue_cox[i]<-summary(coxph(Surv(RFS, dis) ~ pro_res[,i], data = pro_res))$logtest[3]
  exp[i]<-summary(coxph(Surv(RFS, dis) ~ pro_res[,i], data = pro_res))$conf.int[1]
  ecp_low[i]<-summary(coxph(Surv(RFS, dis) ~ pro_res[,i], data = pro_res))$conf.int[3]
  exp_up[i]<-summary(coxph(Surv(RFS, dis) ~ pro_res[,i], data = pro_res))$conf.int[4]
  
}
pvalueAd<-p.adjust(pvalue_cox,method = "BH")
P_cox<-data.frame(cbind(pvalue_cox,pvalueAd,exp,ecp_low,exp_up))
row.names(P_cox)<-colnames(pro_res)[1:(ncol(pro_res)-2)]
# write.csv(P_cox,"20230530_zzovb_ccoc_cov_p.csv",row.names = T)

######################################################
#boxplot
HGSOC<-df_prot[grepl("HGSOC",df_prot$type),]
Normal<-df_prot[grepl("Normal",df_prot$type),]

HGSOC_Nor<-rbind(Normal,HGSOC)
# write.csv(HGSOC_Nor,"sdjfp.csv",row.names = T)
prot_list<-read.xlsx("IHC_antibody230531(2)_1.xlsx",sheetIndex = 1)
df_plot<-HGSOC_Nor[,colnames(HGSOC_Nor)%in%prot_list$prot]
df_plot<-log2(df_plot)
# df_plot1<-HGSOC_Nor[,c(3022,479,3015,449,9,2176)]

df_plot<-cbind(df_plot,HGSOC_Nor[,4328:4334])
df_plot$type<-factor(df_plot$type,levels = c("Normal","HGSOC"))
for (i in 1:6) {
  indsel<-df_plot[,c(i,13)]
  a<-ggplot(indsel,aes(x=type,y=indsel[,1],color=type))+geom_boxplot(size=2)+scale_color_brewer(palette = "Dark2")+
    theme_classic()+geom_jitter(size=2.5)+ylab(colnames(df_plot)[i])
  a
  ggsave(paste0("20230531_zzov_HGSOC_",colnames(df_plot)[i],"_log2.pdf"),a,width = 8,height = 8)
}















