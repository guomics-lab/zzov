
#########################
#boxplot

library(xlsx)
samp_inf<-read.xlsx("20221105ZZOV_patient_sample_info(1)(1).xlsx",sheetIndex = 2)
prot_list_3b<-read.xlsx("20230501ref_exampleProt.xlsx",sheetIndex = 2)
prot_mat<-read.table("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.txt",sep = "\t",header = T)
samp<-read.xlsx("input3_group.xlsx",sheetIndex = 2)
samp$samp<-substr(samp$Sample.name,1,nchar(samp$Sample.name)-1)
pvalue<-read.csv("pvalue.csv",header = T,row.names = 1)

df_prot<-prot_mat[row.names(prot_mat)%in%paste0(prot_list_3b$uniprotID,"_",prot_list_3b$geneName),]
df_prot<-data.frame(t(df_prot))

df_prot$pat_id<-samp$Bcr.patient.barcode[match(row.names(df_prot),samp$samp)]
df_prot$Histology8<-samp_inf$Histology8[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$Group<-samp_inf$Group[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
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
pvalue<-read.csv("pvalue.csv",header = T,row.names = 1)


#3c


prot_list_3c<-read.xlsx("20230501ref_exampleProt.xlsx",sheetIndex = 3)
df_prot<-prot_mat[row.names(prot_mat)%in%paste0(prot_list_3c$uniprotID,"_",prot_list_3c$geneName),]
df_prot<-data.frame(t(df_prot))

df_prot$pat_id<-samp$Bcr.patient.barcode[match(row.names(df_prot),samp$samp)]
df_prot$Histology8<-samp_inf$Histology8[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
df_prot$Group<-samp_inf$Group[match(df_prot$pat_id,samp_inf$Bcr_patient_barcode)]
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


library(ggplot2)
library( RColorBrewer)
for (i in 1:(ncol(df_prot1)-4)) {
  indsel<-df_prot1[,c(i,ncol(df_prot1))]
  prot_p<-data.frame(pvalue[,names(pvalue)==unlist(strsplit(colnames(df_prot1)[i],"_"))[1]])
  row.names(prot_p)<-row.names(pvalue)
  indsel$p<-match(indsel$type,row.names(prot_p))
  indsel$p[indsel$p==1]<-prot_p$pvalue...names.pvalue.....unlist.strsplit.colnames.df_prot1..i...[1]
  indsel$p[indsel$p==2]<-prot_p$pvalue...names.pvalue.....unlist.strsplit.colnames.df_prot1..i...[2]
  indsel$p[indsel$p==3]<-prot_p$pvalue...names.pvalue.....unlist.strsplit.colnames.df_prot1..i...[3]
  indsel$p[indsel$p==4]<-prot_p$pvalue...names.pvalue.....unlist.strsplit.colnames.df_prot1..i...[4]
  indsel$p[indsel$p==5]<-prot_p$pvalue...names.pvalue.....unlist.strsplit.colnames.df_prot1..i...[5]
  indsel$type1<-paste0(indsel$type,"_",indsel$p)
  indsel[,1]<-log2(indsel[,1])
  a<-ggplot(indsel,aes(x=type,y=indsel[,1],color=type))+geom_boxplot(size=2)+scale_color_brewer(palette = "Dark2")+
    theme_classic()+geom_jitter(size=2.5)+ylab(colnames(df_prot1)[i])
  a
  ggsave(paste0("20230506_3C_",colnames(df_prot1)[i],".pdf"),a,width = 8,height = 8)
  b<-ggplot(indsel,aes(x=type1,y=indsel[,1],color=type1))+geom_boxplot(size=2)+scale_color_brewer(palette = "Dark2")+
    theme_classic()+geom_jitter(size=2.5)+ylab(colnames(df_prot1)[i])
  b
  ggsave(paste0("20230506_3C_",colnames(df_prot1)[i],"_1.pdf"),b,width = 8,height = 8)
}




