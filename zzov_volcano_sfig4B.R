library(xlsx)
listfile<-list.files("E:/Qsync/work/zzov/volcano_plot/volcano_input_zz",pattern = ".xlsx")
for (i in listfile) {
  indsel<-read.xlsx(i,sheetIndex = 1,row.names=1)
  pdf(paste0(i,"zzovc_volcan_1.pdf"),width =6,height = 6)
  plot(indsel$log2fc,-log10(indsel$p_value_adjusted),
       pch=19, col="gray",
       xlim=c(-10,10),
       main="volcano plot"
       ,xlab="log2 FC",ylab="-log10 pvalueAD",
       cex.axis = 1.5,cex.lab=1.5)
  up <- which(indsel$p_value_adjusted< 0.05 & indsel$log2fc >1)
  Uppoint1<-which(-log10(indsel$p_value_adjusted)>1.30103&indsel$log2fc > 1)
  down <- which(indsel$p_value_adjusted < 0.05 & indsel$log2fc <  -1)
  Downpoint1<-which(-log10(indsel$p_value_adjusted)>1.30103& indsel$log2fc <  -1)
  points(indsel$log2fc[up], -log10(indsel$p_value_adjusted[up]), col=1, bg = "red", pch=21, cex=1)
   text(indsel$log2fc[Uppoint1], -log10(indsel$p_value_adjusted[Uppoint1]),labels = row.names(indsel)[Uppoint1], col=1, bg = "grey", cex=0.8,font = 3)
  points(indsel$log2fc[down], -log10(indsel$p_value_adjusted[down]), col = "blue", pch=19,cex=1)
   text(indsel$log2fc[Downpoint1]+0.1, -log10(indsel$p_value_adjusted[Downpoint1]),labels =row.names(indsel)[Downpoint1], col=1, bg = "grey", cex=0.8,font = 3)
  abline(h=-log10(0.05),v=c(-1,1),lty=2,lwd=1)
  dev.off()
  
  # protUP<-FC_PAd[up,]
  # protdown<-FC_PAd[down,]
  # ipainput<-rbind(protUP,protdown)
  # write.csv(ipainput,"zzovc_42tiss_sen_vs_res_t_test_diffprot_FC_PAd_p2p_50miss_1_5fc.csv")
}




