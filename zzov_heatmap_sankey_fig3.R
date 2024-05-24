

############################
library(openxlsx)
library(plyr)


############

file_list<-c("20240205_zzov_t-test_ Clear cell carcinoma _vs_Benign_Borderline.xlsx",
             "20240205_zzov_t-test_ Mucinous adenocarcinoma _vs_Benign_Borderline.xlsx",
             "20240205_zzov_t-test_ Endometrioid carcinoma _vs_Benign_Borderline.xlsx",
             "20240205_zzov_t-test_ High grade serous carcinoma _vs_Benign_Borderline.xlsx",
             "20240205_zzov_t-test_ Low grade serous carcinoma _vs_Benign_Borderline.xlsx")

#!is.na(readin$CCOC)&is.na(readin$EMOC)&is.na(readin$HGSOC)&is.na(readin$LGSOC)&is.na(readin$MCOC)



#i=1
i=1
indsel<-read.xlsx(file_list[i],sheet = 1)
indsel<-indsel[!is.na(indsel$p_value),]
up<-indsel[indsel$p_value_adjusted<0.05&indsel$fc_B>1.5,]
down<-indsel[indsel$p_value_adjusted<0.05&indsel$fc_B<2/3,]
temp<-rbind(up,down)
names(temp)<-c("prot",paste0(strsplit(file_list[i],"_")[[1]][4],"p_value"),paste0(strsplit(file_list[i],"_")[[1]][4],"p_value_adjusted"),
               paste0(strsplit(file_list[i],"_")[[1]][4],"fc_B"),paste0(strsplit(file_list[i],"_")[[1]][4],"log2fc"))


for (i in 2:length(file_list)) {
  indsel<-read.xlsx(file_list[i],sheet = 1)
  indsel<-indsel[!is.na(indsel$p_value),]
  up<-indsel[indsel$p_value_adjusted<0.05&indsel$fc_B>1.5,]
  down<-indsel[indsel$p_value_adjusted<0.05&indsel$fc_B<2/3,]
  diff<-rbind(up,down)
  names(diff)<-c("prot",paste0(strsplit(file_list[i],"_")[[1]][4],"p_value"),paste0(strsplit(file_list[i],"_")[[1]][4],"p_value_adjusted"),
                 paste0(strsplit(file_list[i],"_")[[1]][4],"fc_B"),paste0(strsplit(file_list[i],"_")[[1]][4],"log2fc"))
  temp<-merge(temp,diff,by="prot",all=T)
  
}
# unique(temp$prot)
write.csv(temp,"20240207_zzov_x_vs_Benign_Borderline_difprot_1_5fc.csv",row.names = F)

df1<-read.xlsx("Supplementary Table 4_histotype specific DEPs20240206.xlsx",sheet = 4)

temp$prot_name<-sapply(strsplit(temp$prot,"_"),function(e){e[1]})
upset_prot<-temp[temp$prot_name%in%df1$UniprotID,]


library(UpSetR)

listinput <- list(Clear_cell_carcinoma=upset_prot$prot_name[!is.na(upset_prot$` Clear cell carcinoma p_value`)],
                  Mucinous_adenocarcinoma=upset_prot$prot_name[!is.na(upset_prot$` Mucinous adenocarcinoma p_value`)],
                  Endometrioid_carcinoma=upset_prot$prot_name[!is.na(upset_prot$` Endometrioid carcinoma p_value`)],
                  High_grade_serous_carcinoma=upset_prot$prot_name[!is.na(upset_prot$` High grade serous carcinoma p_value`)],
                  Low_grade_serous_carcinoma=upset_prot$prot_name[!is.na(upset_prot$` Low grade serous carcinoma p_value`)])

# pdf('20240207_zzov_x_vs_Benign_Borderline_difprot_overlap.pdf',height = 8,width = 8)
# upset(fromList(listinput),nsets =5, order.by = "freq",text.scale =1.5)
# dev.off()


# CCOC=Clear cell carcinoma
# EMOC=Endometrioid carcinoma
# HGSOC=High grade serous carcinoma
# LGSOC=Low grade serous carcinoma
# MCOC=Mucinous adenocarcinoma

df2<-read.xlsx("20240205_zzov_t-test_ carcinoma _vs_Benign_Borderline_sum.xlsx",sheet =8 )
df2$label<-gsub("MC","MCOC",df2$label)
df2$label<-gsub("EC","EMOC",df2$label)

upset_uni<-data.frame(upset_prot$prot_name[!is.na(upset_prot$` Clear cell carcinoma p_value`)&is.na(upset_prot$` Mucinous adenocarcinoma p_value`)&
                                             is.na(upset_prot$` Endometrioid carcinoma p_value`)&is.na(upset_prot$` High grade serous carcinoma p_value`)&
                                             is.na(upset_prot$` Low grade serous carcinoma p_value`)])

upset_uni$type<-"CCOC"
names(upset_uni)[1]<-"prot"


upset_uni1<-data.frame(upset_prot$prot_name[is.na(upset_prot$` Clear cell carcinoma p_value`)&!is.na(upset_prot$` Mucinous adenocarcinoma p_value`)&
                                             is.na(upset_prot$` Endometrioid carcinoma p_value`)&is.na(upset_prot$` High grade serous carcinoma p_value`)&
                                             is.na(upset_prot$` Low grade serous carcinoma p_value`)])

upset_uni1$type<-"MCOC"
names(upset_uni1)[1]<-"prot"


upset_uni2<-data.frame(upset_prot$prot_name[is.na(upset_prot$` Clear cell carcinoma p_value`)&is.na(upset_prot$` Mucinous adenocarcinoma p_value`)&
                                              !is.na(upset_prot$` Endometrioid carcinoma p_value`)&is.na(upset_prot$` High grade serous carcinoma p_value`)&
                                              is.na(upset_prot$` Low grade serous carcinoma p_value`)])

upset_uni2$type<-"EMOC"
names(upset_uni2)[1]<-"prot"


upset_uni3<-data.frame(upset_prot$prot_name[is.na(upset_prot$` Clear cell carcinoma p_value`)&is.na(upset_prot$` Mucinous adenocarcinoma p_value`)&
                                              is.na(upset_prot$` Endometrioid carcinoma p_value`)&!is.na(upset_prot$` High grade serous carcinoma p_value`)&
                                              is.na(upset_prot$` Low grade serous carcinoma p_value`)])

upset_uni3$type<-"HGSOC"
names(upset_uni3)[1]<-"prot"


upset_uni4<-data.frame(upset_prot$prot_name[is.na(upset_prot$` Clear cell carcinoma p_value`)&is.na(upset_prot$` Mucinous adenocarcinoma p_value`)&
                                              is.na(upset_prot$` Endometrioid carcinoma p_value`)&is.na(upset_prot$` High grade serous carcinoma p_value`)&
                                              !is.na(upset_prot$` Low grade serous carcinoma p_value`)])

upset_uni4$type<-"LGSOC"
names(upset_uni4)[1]<-"prot"

upset_uni<-rbind(upset_uni,upset_uni1,upset_uni2,upset_uni3,upset_uni4)

setdiff(upset_uni$prot,df2$uniprotID)
setdiff(df2$uniprotID,upset_uni$prot)

setdiff(upset_uni$prot[upset_uni$type=="CCOC"],df2$uniprotID[df2$label=="CCOC"])
setdiff(upset_uni$prot[upset_uni$type=="MCOC"],df2$uniprotID[df2$label=="MCOC"])
setdiff(upset_uni$prot[upset_uni$type=="EMOC"],df2$uniprotID[df2$label=="EMOC"])
setdiff(upset_uni$prot[upset_uni$type=="HGSOC"],df2$uniprotID[df2$label=="HGSOC"])
setdiff(upset_uni$prot[upset_uni$type=="LGSOC"],df2$uniprotID[df2$label=="LGSOC"])



####################
#heatmap

prot_mat<-read.xlsx("Supplementary Table 2_global proteome20230416.xlsx",sheet = 2,rowNames = T)
prot_mat$prot<-sapply(strsplit(row.names(prot_mat),"_"),function(e){e[1]})
hp_mat<-prot_mat[prot_mat$prot%in%upset_uni$prot,]

patinf<-read.xlsx("20210421ZZOV_patient_923sample_info.xlsx",sheet = 1)
df4<-read.xlsx("Supplementary Table 1_patient sample_info20230416.xlsx",sheet = 3)
name_prot<-df4$Bcr.patient.barcode[match(names(hp_mat)[1:1041],df4$Sample.name)]

hp_dat<-data.frame(matrix(NA,nrow = nrow(hp_mat),ncol = 6))
row.names(hp_dat)<-row.names(hp_mat)
names(hp_dat)<-c("Benign&borderline","CCOC","EMOC","HGSOC","LGSOC","MCOC")
hp_mat1<-hp_mat[,-ncol(hp_mat)]
hp_dat$`Benign&borderline`<-rowMeans(hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="NM"]],na.rm = T)
hp_dat$CCOC<-rowMeans(hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="CC"]],na.rm = T)
hp_dat$EMOC<-rowMeans(hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="EC"]],na.rm = T)
hp_dat$HGSOC<-rowMeans(hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="HS"]],na.rm = T)
hp_dat$LGSOC<-rowMeans(hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="LS"]],na.rm = T)
hp_dat$MCOC<-rowMeans(aaa<-hp_mat1[,name_prot%in%patinf$Bcr_patient_barcode[patinf$Histology=="MC"]],na.rm = T)
hp_dat[is.na(hp_dat)]<-NA
# unique(name_prot)
# setdiff(unique(name_prot),patinf$Bcr_patient_barcode)
library(pheatmap)

### Z Score normalization
centermatrix <- scale(t(hp_dat))
testdata <- data.frame(centermatrix)
# testdata = testdata[match(anncol1$SampleName, row.names(testdata)),]
# names(testdata) <- names(df_filter)

o1 <- (fivenum(unlist(testdata))[4]-fivenum(unlist(testdata))[2])*2+fivenum(unlist(testdata))[4]
o2 <- fivenum(unlist(testdata))[2]-(fivenum(unlist(testdata))[4]-fivenum(unlist(testdata))[2])*2

library(magrittr)
library(RColorBrewer)
testdata[testdata>o1]=o1 %>%as.numeric()
testdata[testdata<o2]=o2 %>%as.numeric() 

row.names(testdata)<-c("Benign&borderline","CC","EC","HS","LS","MC")
# testdata2 <- t(testdata)

paletteLength = 11
myColor <- c(brewer.pal(11,"RdYlBu")[11:7],"azure1",brewer.pal(11,"RdYlBu")[5:1])
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
# myBreaks <- c(seq(min(testdata), 0, length.out=ceiling(paletteLength/2) + 1), http://127.0.0.1:32293/graphics/plot_zoom_png?width=2544&height=1338
#               seq(max(testdata)/paletteLength, max(testdata), length.out=floor(paletteLength/2)))
myBreaks = seq(-3.5,3.5,length.out = 12)

anncol2<-data.frame(row.names(testdata) )
names(anncol2)<-"Histology"
row.names(anncol2)<-anncol2$Histology
testdata2 <- as.matrix(testdata)
#Change the site of EC and HS
testdata2<-testdata2[c(1,2,4,3,5,6),]
ann_colors <- list(Histology = c(`Benign&borderline` = "#63B797" ,HS="#C75E1B",LS="#7E7EA8",MC="#C4518A",EC="#5F9624",CC="#E5AC39"))
pdf("20240207_zzov_heatmap_5_mean_4.pdf",
    height = 8,width = 15)
p<-pheatmap(testdata2,
         cellwidth = 1, cellheight = 6,
         # fontsize_col = 3,
         color = c(brewer.pal(11,"RdYlBu")[11:7],"azure1",brewer.pal(11,"RdYlBu")[5:1]),
         annotation_row  = anncol2,
         annotation_colors = ann_colors,
         na_col = "#CCCCCC",
         #scale = "row", #nomlization
         scale = "none",
         cluster_rows = T, 
         cluster_cols = T,
         treeheight_col=5,
         treeheight_row=5,
         clustering_method = "ward.D2",
         show_rownames = T, 
         show_colnames = T,
         fontsize_row = 4,
         fontsize_col = 1,
         # cutree_rows = 8,
         border_color = NA,
         #filename = "allprotein_heatmap.pdf",       
         main = "")
dev.off()
# prot_names_new= names(testdata2[,p$tree_col[["order"]]])
row.order = p$tree_row$order 
col.order = testdata2[,p$tree_col$order]
write.csv(col.order,"20240207_zzov_prot_pheatmap_list_new.csv",row.names = T)


###############
#20240218 sankey
library(openxlsx)
sank_list<-read.xlsx("20240201_zzov_prot_pheatmap_list(1).xlsx",sheet = 2)
sank_list$pathway_clus<-sapply(strsplit(sank_list$pathway2,"\\."),function(e){e[1]})
sank_list$pathway_clus1<-sapply(strsplit(sank_list$pathway2,"\\."),function(e){e[2]})
# sank_list$pathway_clus2<-gsub("[a-z,A-Z]","",sank_list$pathway_clus1)
# sank_list$pathway_clus1<-gsub("[0-9]","",sank_list$pathway_clus1)
sank_list$num<-1
sank_list$UniprotID<-factor(sank_list$UniprotID,levels = sank_list$UniprotID)
sank_list$pathway_clus2<-sank_list$pathway_clus
# sank_list$pathway_clus1<-gsub("[0-9]NA",NA,sank_list$pathway_clus1)
sank_list$pathway_clus2[is.na(sank_list$pathway)]<-NA
# sank_list$pathway2<-factor(sank_list$pathway2,levels = sank_list$pathway2)
# sank_plot<-sank_list[,c(1,13,14,15)]
# write.csv(sank_plot,"sank_plot.csv",row.names = F)
# 
# 
# sank_plot<-read.csv("sank_plot.csv",header = T)
# sank_plot$UniprotID<-as.factor(sank_plot$UniprotID)
# sank_plot$num<-as.numeric(sank_plot$num)
# sank_plot$type<-as.factor(sank_plot$type)
# sank_plot$des<-as.numeric(sank_plot$des)
# sank_plot$des1<-as.numeric(sank_plot$des1)
# class(sank_plot$des1)
#sanky
library(reshape2)

library(ggplot2)
library(ggalluvial)
library(dplyr)

a<-ggplot(sank_list,
          aes(y = num, axis1 = UniprotID, axis2 = pathway2)) +
  geom_alluvium(aes(fill = pathway_clus2), width = 0.2) +
  geom_stratum(width = 0.2, fill = "white", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=2,color="black")+
  #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Uniprot", "Pathway"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") 
a
ggsave("20240220_zzov_sankey.pdf",a,width = 8,height = 8)


