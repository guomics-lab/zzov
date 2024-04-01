library(xlsx)
# readin<-read.csv("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.csv",header = T,row.names = 1)
readin<-read.table("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.txt",header = T,row.names = 1,sep = "\t",check.names = F)
# min_mat<-2389.27
prot_list<-read.xlsx("20230423DNAdamage_repairFig7E_input.xlsx",sheetIndex = 1)

heat_prot<-data.frame(t(readin))
sample_inf<-read.xlsx("20210728ZZOV_patient_sample_info.xlsx",sheetIndex = 4)
heat_prot$patiend_id<-sample_inf$bcr_patient_barcode[match(row.names(heat_prot),sample_inf$Sample_name)]
heat_prot1<-heat_prot[heat_prot$patiend_id%in%prot_list$samp,]
heat_prot1$type<-prot_list$type[match(heat_prot1$patiend_id,prot_list$samp)]


inf<-read.xlsx("20230423DNAdamage_repairFig7E_input.xlsx",sheetIndex = 2)
inf$prot<-paste(inf$uniprot,inf$gene,sep = "_")
heat_prot2<-heat_prot1[,names(heat_prot1)%in%inf$prot]
heat_prot2<-cbind(heat_prot2,heat_prot1[,10528:10529])
heat_prot2<-heat_prot2[order(heat_prot2$type),]
heat_prot3<-data.frame(t(heat_prot2[,1:18]))
heat_prot3$pathway1<-inf$pathwayLabel1[match(row.names(heat_prot3),inf$prot)]
heat_prot3$pathway2<-inf$pathwayLabel2[match(row.names(heat_prot3),inf$prot)]


heat_prot3[,1:42]<-log2(heat_prot3[,1:42])
plot(density(unlist(heat_prot3[,1:42]),na.rm = T))
annotation_row<-data.frame(heat_prot3[,43:44])
row.names(annotation_row)<-row.names(heat_prot3)

annotation_col<-data.frame(heat_prot2$type)
row.names(annotation_col)<-row.names(heat_prot2)
bk <- c(seq(-2.5,-0.1,by=0.01),seq(0,2.5,by=0.01))
pdf("20230423_zzov_fig7E_1.pdf",width = 10,height = 8)
pheatmap::pheatmap(heat_prot3[,1:42],cluster_rows = T,cluster_cols = F,scale = "row",
                   color =c(colorRampPalette(colors = c("#004C99","white"))(length(bk)/2),colorRampPalette(colors = c("white","#990000"))(length(bk)/2)),legend_breaks=seq(-2.5,2.5,1),
                   breaks=bk,
                   border_color = "NA",annotation_row = annotation_row,annotation_col=annotation_col,annotation_legend = F)
dev.off()

table(annotation_col$heat_prot2.type)


########################################
#fig7E2
library(xlsx)
# readin<-read.csv("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.csv",header = T,row.names = 1)
readin<-read.table("20210615ZZOV_923sample_quanNor_NA=0.8min_diann_prot.txt",header = T,row.names = 1,sep = "\t",check.names = F)
# min_mat<-2389.27
prot_list<-read.xlsx("20230423DNAdamage_repairFig7E_input.xlsx",sheetIndex = 3)

heat_prot<-data.frame(t(readin))
sample_inf<-read.xlsx("20210728ZZOV_patient_sample_info.xlsx",sheetIndex = 4)
heat_prot$patiend_id<-sample_inf$bcr_patient_barcode[match(row.names(heat_prot),sample_inf$Sample_name)]
heat_prot1<-heat_prot[heat_prot$patiend_id%in%prot_list$samp,]
heat_prot1$type<-prot_list$type[match(heat_prot1$patiend_id,prot_list$samp)]


inf<-read.xlsx("20230423DNAdamage_repairFig7E_input.xlsx",sheetIndex = 4)
inf$prot<-paste(inf$uniprot,inf$gene,sep = "_")
heat_prot2<-heat_prot1[,names(heat_prot1)%in%inf$prot]
heat_prot2<-cbind(heat_prot2,heat_prot1[,10528:10529])
heat_prot2<-heat_prot2[order(heat_prot2$type,decreasing = T),]
heat_prot3<-data.frame(t(heat_prot2[,1:18]))
heat_prot3$pathway<-inf$pathway[match(row.names(heat_prot3),inf$prot)]


heat_prot3[,1:46]<-log2(heat_prot3[,1:46])
plot(density(unlist(heat_prot3[,1:46]),na.rm = T))
annotation_row<-data.frame(heat_prot3$pathway)
row.names(annotation_row)<-row.names(heat_prot3)

annotation_col<-data.frame(heat_prot2$type)
row.names(annotation_col)<-row.names(heat_prot2)
bk <- c(seq(-2.5,-0.1,by=0.01),seq(0,2.5,by=0.01))
pdf("20230423_zzov_fig7E2.pdf",width = 10,height = 8)
pheatmap::pheatmap(heat_prot3[,1:46],cluster_rows = T,cluster_cols = F,scale = "row",
                   color =c(colorRampPalette(colors = c("#004C99","white"))(length(bk)/2),colorRampPalette(colors = c("white","#990000"))(length(bk)/2)),
                   legend_breaks=seq(-2.5,2.5,1),
                   breaks=bk,
                   border_color = "NA",annotation_row = annotation_row,annotation_col=annotation_col,annotation_legend = F)
dev.off()

table(annotation_col$heat_prot2.type)




























