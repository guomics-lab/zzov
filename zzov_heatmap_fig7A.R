
#########
#20240130
library(openxlsx)
df1<-read.xlsx("E:/Qsync/work/zzov/20230409_figure/table/2021.9.2  300基因测序-germline mutation list20230418update.xlsx",sheet = 2)
df3<-read.xlsx("E:/Qsync/work/zzov/20230409_figure/table/Supplementary Table 2_global proteome20230416.xlsx",sheet = 2,rowNames = T)
df4<-read.xlsx("E:/Qsync/work/zzov/20230409_figure/table/Supplementary Table 1_patient sample_info20230416.xlsx",sheet = 3)
df3_name<-data.frame(names(df3))
df3_name$type<-df4$Bcr.patient.barcode[match(df3_name$names.df3.,df4$Sample.name)]
df3_name$label<-df1$label[match(df3_name$type,df1$bcr_patient_barcode)]
df1$proteome<-df3_name$names.df3.[match(df1$bcr_patient_barcode,df3_name$type)]


###############
#20230426

library(reshape2)

all_gene<-read.xlsx("20230422ZZOV_genome_summary_fig7A.xlsx",sheet = 6)
# all_gene<-all_gene[-1,]
# top_gene<-all_gene[all_gene$count0_WT>=5,]
#del MUC16
all_gene<-all_gene[!grepl("MUC16",all_gene$X1),]
all_gene$X1<-gsub("_summary","",all_gene$X1)
germiline<-read.xlsx("20230422ZZOV_genome_summary_fig7A.xlsx",sheet = 3)

ger<-germiline[,c(1,3,4,7)]

somatic<-read.xlsx("20230422ZZOV_genome_summary_fig7A.xlsx",sheet = 4)
som<-somatic[,c(2,4,5,12)]
names(som)<-names(ger)

ger_som<-rbind(ger,som)
ger_som1<-dcast(ger_som,`#Sample_ID`+Variant_Classification~Gene_Name)
ger_som2<-ger_som1[,names(ger_som1)%in%all_gene$X1]
ger_som2<-cbind(ger_som1[,1:2],ger_som2)
ger_som2$group<-ger_som$group[match(ger_som2$`#Sample_ID`,ger_som$`#Sample_ID`)]
write.csv(ger_som2,"20230426_zzoc_fig7A_all.csv",row.names = F)

#Change 2 to 1 in ger_som2 ("2" represents a gene that is present in both germline and somatic cells (of the same type).)
ger_som2<-read.csv("20230426_zzoc_fig7A_all_1.csv",header = T)
for (i in 1:nrow(ger_som2)) {
  ger_som2[i,3:152]<-gsub(1,ger_som2[i,2],ger_som2[i,3:152])
}

ger_som_uni<-data.frame(matrix(NA,nrow = length(unique(ger_som2$X.Sample_ID)),ncol = 150))
row.names(ger_som_uni)<-unique(ger_som2$X.Sample_ID)
names(ger_som_uni)<-names(ger_som2)[3:152]
for (i in 1:nrow(ger_som_uni)) {
  indsel<-ger_som2[ger_som2$X.Sample_ID==row.names(ger_som_uni)[i],3:152]
  df<-c()
  for (j in 1:nrow(indsel)) {
   df<-paste0(df,indsel[j,],sep=";") 
  }
  ger_som_uni[i,]<-df
}
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("0;","",e)}))
row.names(ger_som_uni)<-unique(ger_som2$X.Sample_ID)
unique(ger_som2$Variant_Classification)
ger_som_uni$group<-ger_som$group[match(row.names(ger_som_uni),ger_som$`#Sample_ID`)]
# write.csv(ger_som_uni,"20230425_zzov_fig7A_gene_summary.csv",row.names = T)

ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("frameshift_variant;",1,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("missense_variant;",2,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("splice_donor_variant;",3,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("inframe_deletion;",4,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("splice_acceptor_variant;",5,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("stop_gained;",6,e)}))
ger_som_uni<-data.frame(sapply(ger_som_uni,function(e){gsub("stop_lost;",7,e)}))

ger_som_uni$TP53<-gsub(23,8,ger_som_uni$TP53)
ger_som_uni$TP53<-gsub(26,8,ger_som_uni$TP53)
ger_som_uni$FANCI<-gsub(25,8,ger_som_uni$TP53)
ger_som_uni$ATM<-gsub(12,8,ger_som_uni$ATM)
row.names(ger_som_uni)<-unique(ger_som2$X.Sample_ID)
group1<-ger_som_uni[ger_som_uni$group=="primary_sensitive",]
group2<-ger_som_uni[ger_som_uni$group=="primary_resistant",]
group3<-ger_som_uni[ger_som_uni$group=="relapsing_sensitive",]
group4<-ger_som_uni[ger_som_uni$group=="relapsing_resistant",]
ger_som_uni<-rbind(group1,group2,group3,group4)


# write.csv(ger_som_uni,"test.csv",row.names = T)
# ger_som_uni<-read.csv("test.csv",header = T,row.names = 1)
HR_path<-read.xlsx("20230422ZZOV_genome_summary_fig7A.xlsx",sheet= 7)
ger_som_uni_hr_num<-ger_som_uni[,names(ger_som_uni)%in%HR_path$Gene]
ger_som_uni_hr_num$hr_num<-apply(!is.na(ger_som_uni_hr_num),1,sum)
ger_som_uni_hr_num$hr_num[ger_som_uni_hr_num$hr_num>1]<-8
#The two mutations of OV029 are of the same type, so they should be changed to its own type.
ger_som_uni_hr_num$hr_num[row.names(ger_som_uni_hr_num)=="OV029"]<-6
ger_som_uni_hr_num$max<-apply(ger_som_uni_hr_num[,1:11],1,max,na.rm=T)

ger_som_uni_hr_num$hr_num[ger_som_uni_hr_num$hr_num==1]<-ger_som_uni_hr_num$max[ger_som_uni_hr_num$hr_num==1]


ger_som_uni$hr_num<-ger_som_uni_hr_num$hr_num
ger_som_uni1<-data.frame(t(ger_som_uni[,c(1:150,152)]))
ger_som_uni1[is.na(ger_som_uni1)]<-0
top_gene<-all_gene[all_gene$count0_WT>=5,]
ger_som_uni2<-ger_som_uni1[row.names(ger_som_uni1)%in%top_gene$X1,]
ger_som_uni2["HR_pathway",]<-ger_som_uni$hr_num
ger_som_uni2$HR_path<-0
ger_som_uni2$HR_path[match(HR_path$Gene,row.names(ger_som_uni2))]<-1
ger_som_uni2<-ger_som_uni2[c(17,3,1,4,14,7,8,9,11,6,10,12,13,15,16,2,5,18),]
# write.csv(ger_som_uni2,"20240206_ger_som_uni2.csv",row.names = T)
# ger_som_uni2<-read.csv("20240206_ger_som_uni2.csv",header = T,row.names = 1)

#1）Put all HRR genes together after TP53, arranged according to their mutation percentages.
#2）Sort the samples based on whether the HRR gene is mutated, with the mutated samples grouped together and the non-mutated samples placed afterwards. 
#For the samples with mutations, prioritize the gene with the highest mutation rate and place it at the front.

ger_som_uni3<-data.frame(sapply(ger_som_uni2,function(e){gsub("[1-9]",1,e)}))
# write.csv(ger_som_uni3,"test2.csv",row.names = F)
# 
# ger_som_uni3<-read.csv("test2.csv",header = T)
ger_som_uni3$pesent<-round(apply(ger_som_uni3[c(1:53,55:97)],1,sum)/96,3)*100
row.names(ger_som_uni3)<-row.names(ger_som_uni2)


row.names(ger_som_uni2)[c(1:18)]<-paste0(row.names(ger_som_uni3),"_",ger_som_uni3$pesent,"%")   


annotation_row<-data.frame(ger_som_uni3$HR_path)

row.names(annotation_row)<-row.names(ger_som_uni2)
annotation_col<-data.frame(ger_som_uni$group)

row.names(annotation_col)<-row.names(ger_som_uni)
annotation_col$port<-df1$proteome[match(row.names(annotation_col),df1$label)]
annotation_col$port[!is.na(annotation_col$port)]<-"yes"
annotation_col$port[is.na(annotation_col$port)]<-"no"
# write.csv(annotation_col,"20240206_annotation_col.csv",row.names = T)
pdf("20240324_zzov_fig7A.pdf",width = 10,height = 8)
pheatmap::pheatmap(ger_som_uni2[,1:97],cluster_rows = F,cluster_cols = F,
                   color =c("white","red","darkseagreen3","orange","skyblue","darkblue","coral2","green","chocolate4","cornflowerblue"),
                   legend_breaks=seq(0,9,1),border = T,
                   border_color = "black",annotation_row=annotation_row,annotation_col=annotation_col,annotation_legend = T)
dev.off()

#
table(annotation_col$ger_som_uni.group)





























