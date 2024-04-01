options(digits = 4)
library(xlsx)
readin<-read.xlsx("20211228ZZOVplasma_70%NA_batchfree_ComRep_168samp_1660Prot.xlsx",sheetIndex = 1)
readinf<-read.xlsx("20211228ZZOVplasma_70%NA_batchfree_ComRep_168samp_1660Prot.xlsx",sheetIndex = 3)
readin$label<-readinf$label[match(readinf$patientID,readin$SwissprotID)]
protlist<-read.xlsx("8feature.xlsx",sheetIndex = 1)

#HS=1 BN&NM=0
readprot<-readin[,colnames(readin)%in%protlist$prot]
readprot<-cbind(readprot,readin[,c(1,ncol(readin))])
readprot[is.na(readprot)]<-0
readprot$label<-as.factor(readprot$label)
row.names(readprot)<-readprot$SwissprotID
readprot<-readprot[,-9]
#####################################################################################
#Machine learning

#Random forest

library(readr)
library(readxl)
library(stringr)
library(magrittr)
library("randomForest")
library(pROC)
library(caret)

accu<-pred_list<-c()

for (i in 1:8) {
  indsel<-data.frame(cbind(readprot[,i],readprot$label))
  row.names(indsel)<-row.names(readprot)
  names(indsel)<-c(names(readprot)[i],"label")
  indsel$label<-as.factor(indsel$label)
  set.seed(2022)
  
  folds <- createFolds(indsel$label,5)
  n=0
  for(fold in folds){
    n=n+1
    #fold=folds[[8]]
    valids <- indsel[fold,]
    # valids$label <- train_set$label[fold]
    trains <- indsel[setdiff(1:dim(indsel)[1],fold),]
    # trains$label <- indsel$label[setdiff(1:dim(train_set)[1],fold)]
    # trains$label <- as.factor(trains$label)
    # for (ntree in seq(600,1000,200)) {
    set.seed(2022.12)
    tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
    
        predicted <- predict(tmpRF,valids,type='prob')
        predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
        colnames(predicts) <- colnames(predicted)
        predicts <- data.frame(predicts,check.names=F)
        predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
        predicts$observed <- valids$label
        predicts$samp_id<-row.names(valids)
        pred_list<-rbind(pred_list,predicts)
        ROC <- roc(predicts$observed, as.numeric(predicts$`1`))
        auc <- as.numeric(ROC[["auc"]])
        acc <- sum(predicts$predicted==predicts$observed)
        accu <- rbind(accu,c(names(readprot)[i],n,acc/length(fold),auc))
      }
}
accu<-data.frame(accu)
names(accu)<-c("prot","fold","acc","auc")

set.seed(2022)

folds <- createFolds(indsel$label,5)
n=0
fea_total<-pred_list_total<-c()
for(fold in folds){
  n=n+1
  #fold=folds[[8]]
  # fold=c(1,3,4 ,11,20,29,34,47,48,54,59,61,68,82,84,85,92,93,97,106,107,113,115,116,124,131,132,138,139,150,152,153,154,168)
  valids <- readprot[fold,]
  # valids$label <- train_set$label[fold]
  trains <- readprot[setdiff(1:dim(readprot)[1],fold),]
  length(trains$label[trains$label==1])
  # trains$label <- indsel$label[setdiff(1:dim(train_set)[1],fold)]
  # trains$label <- as.factor(trains$label)
  # for (ntree in seq(600,1000,200)) {
  set.seed(2022.12)
  tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
  fea <- data.frame(importance(tmpRF,type=1))
  fea$feature<-row.names(fea)
  fea_total<-rbind(fea_total,fea)
  predicted <- predict(tmpRF,valids,type='prob')
  predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
  colnames(predicts) <- colnames(predicted)
  predicts <- data.frame(predicts,check.names=F)
  predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
  predicts$observed <- valids$label
  pred_list_total<-rbind(pred_list_total,predicts)
  ROC <- roc(predicts$observed, as.numeric(predicts$`1`))
  auc <- as.numeric(ROC[["auc"]])
  acc <- sum(predicts$predicted==predicts$observed)
  accu <- rbind(accu,c("allprot",n,acc/length(fold),auc))
  
}
fea_total<-data.frame(fea_total)
write.csv(accu,"20221206_9model_5fold.csv",row.names = F)

##########################################################
#ROC plot
ROC_all <- roc(pred_list_total$observed, as.numeric(pred_list_total$`1`))
auc_all <- as.numeric(ROC_all[["auc"]])
acc_all <- sum(pred_list_total$predicted==pred_list_total$observed)/nrow(pred_list_total)

list_lab<-list()
for (i in 1:8) {
  list_lab[[i]]<-pred_list[((i-1)*168+1):(i*168),]
  
}

#1

ROC_1 <- roc(list_lab[[1]]$observed, as.numeric(list_lab[[1]]$`1`))
auc_1 <- as.numeric(ROC_1[["auc"]])
acc_1 <- sum(list_lab[[1]]$predicted==list_lab[[1]]$observed)/nrow(list_lab[[1]])

#2

ROC_2 <- roc(list_lab[[2]]$observed, as.numeric(list_lab[[2]]$`1`))
auc_2 <- as.numeric(ROC_2[["auc"]])
acc_2 <- sum(list_lab[[2]]$predicted==list_lab[[2]]$observed)/nrow(list_lab[[2]])

#3

ROC_3 <- roc(list_lab[[3]]$observed, as.numeric(list_lab[[3]]$`1`))
auc_3 <- as.numeric(ROC_3[["auc"]])
acc_3 <- sum(list_lab[[3]]$predicted==list_lab[[3]]$observed)/nrow(list_lab[[3]])

#4

ROC_4 <-roc(list_lab[[4]]$observed, as.numeric(list_lab[[4]]$`1`))
auc_4 <-as.numeric(ROC_4[["auc"]])
acc_4 <-sum(list_lab[[4]]$predicted==list_lab[[4]]$observed)/nrow(list_lab[[4]])

#5

ROC_5 <- roc(list_lab[[5]]$observed, as.numeric(list_lab[[5]]$`1`))
auc_5 <- as.numeric(ROC_5[["auc"]])
acc_5 <- sum(list_lab[[5]]$predicted==list_lab[[5]]$observed)/nrow(list_lab[[5]])

#6

ROC_6 <- roc(list_lab[[6]]$observed, as.numeric(list_lab[[6]]$`1`))
auc_6 <- as.numeric(ROC_6[["auc"]])
acc_6 <- sum(list_lab[[6]]$predicted==list_lab[[6]]$observed)/nrow(list_lab[[6]])

#7

ROC_7 <- roc(list_lab[[7]]$observed, as.numeric(list_lab[[7]]$`1`))
auc_7 <- as.numeric(ROC_7[["auc"]])
acc_7 <- sum(list_lab[[7]]$predicted==list_lab[[7]]$observed)/nrow(list_lab[[7]])

#8
ROC_8 <- roc(list_lab[[8]]$observed, as.numeric(list_lab[[8]]$`1`))
auc_8 <- as.numeric(ROC_8[["auc"]])
acc_8 <- sum(list_lab[[8]]$predicted==list_lab[[8]]$observed)/nrow(list_lab[[8]])

###########################################
#roc
pdf("20230220_zzovc_ROC_zzov_plasma_9model.pdf", width = 8, height = 6)
plot(ROC_all, print.auc = TRUE, print.thres = F, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve of 9 model", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5)

plot(ROC_1, add = T, col = 'blue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_2, add = T, col = 'black', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_3, add = T, col = 'pink', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_4, add = T, col = 'green', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_5, add = T, col = 'orange', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_6, add = T, col = 'purple', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_7, add = T, col = 'darkblue', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_8, add = T, col = 'gray', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)

legend("bottomright", legend = c(paste0("auc_all:",auc_all),paste0("SAA2:",auc_1),paste0("CHI3L1:",auc_2),
                                 paste0("MUC16:",auc_3),paste0("SPINT1:",auc_4),paste0("HPSE:",auc_5),
                                 paste0("MUC1:",auc_6),paste0("SMPDL3B:",auc_7),paste0("MXRA5:",auc_8)),
       col=c("red", "blue",'black','pink','green','orange','purple','darkblue','gray'), lwd=2)
dev.off()

######################
#roc compair
roc_temp<-data.frame(matrix(NA,nrow = 14,ncol =2 ))
# row.names(roc_temp)<-c("All_SAA2","All_CHI3L1","All_MUC16","All_SPINT1","All_HPSE","All_MUC1","All_SMPDL3B","All_MXRA5",
#                        "All_CHI3L1_7","All_MUC1_6","All_SMPDL3B_5","All_SAA2_4","All_MXRA5_3","All_HPSE_2")
names(roc_temp)<-c("model","p_value")
roc_temp$model<-c("All_SAA2","All_CHI3L1","All_MUC16","All_SPINT1","All_HPSE","All_MUC1","All_SMPDL3B","All_MXRA5",
                   "All_CHI3L1_7","All_MUC1_6","All_SMPDL3B_5","All_SAA2_4","All_MXRA5_3","All_HPSE_2")
set.seed(123)
roc_temp[1,2]<-roc.test(ROC_all , ROC_1 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[2,2]<-roc.test(ROC_all , ROC_2 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[3,2]<-roc.test(ROC_all , ROC_3 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[4,2]<-roc.test(ROC_all , ROC_4 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[5,2]<-roc.test(ROC_all , ROC_5 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[6,2]<-roc.test(ROC_all , ROC_6 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[7,2]<-roc.test(ROC_all , ROC_7 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[8,2]<-roc.test(ROC_all , ROC_8 , reuse.auc=FALSE)$p.value



#############################
fea_score<-data.frame(matrix(NA,nrow = length(unique(fea_total$feature)),ncol = 1))
row.names(fea_score)<-unique(fea_total$feature)
names(fea_score)<-"total_score"
for (i in 1:nrow(fea_score)) {
  indsel<-fea_total[fea_total$feature==row.names(fea_score)[i],]
  fea_score$total_score[i]<-mean(indsel$MeanDecreaseAccuracy)
  
  
}

fea_score1<-fea_score
fea_score1$feature<-row.names(fea_score1)
fea_score1<-fea_score1[order(fea_score1$total_score),]
fea_score1$order<-c(1:8)
write.csv(fea_score1,"20230221_total_meandecacc_RF_important.csv",row.names = F)
pdf("20230221_total_meandecacc_RF_important.pdf",width = 4,height = 8)
plot(x=fea_score1$total_score,y=fea_score1$order,ylim = c(1,8))
dev.off()

prot_order<-row.names(fea_score)[order(fea_score$total_score)]
accu1<-pred_list1<-c()
for (i in 1:(length(prot_order)-2)) {
  indsel<-data.frame(cbind(readprot[,prot_order[(i+1):8]],readprot$label))
  row.names(indsel)<-row.names(readprot)
  names(indsel)<-c(prot_order[(i+1):8],"label")
  feat_selt<-names(indsel)[1:(ncol(indsel)-1)]
  
  indsel$label<-as.factor(indsel$label)
  set.seed(2022)
  
  folds <- createFolds(indsel$label,5)
  n=0
  for(fold in folds){
    n=n+1
    #fold=folds[[8]]
    valids <- indsel[fold,]
    # valids$label <- train_set$label[fold]
    trains <- indsel[setdiff(1:dim(indsel)[1],fold),]
    # trains$label <- indsel$label[setdiff(1:dim(train_set)[1],fold)]
    # trains$label <- as.factor(trains$label)
    # for (ntree in seq(600,1000,200)) {
    set.seed(2022.12)
    tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
    
    predicted <- predict(tmpRF,valids,type='prob')
    predicts <- t(apply(predicted,1,function(v){v/sum(v)}))
    colnames(predicts) <- colnames(predicted)
    predicts <- data.frame(predicts,check.names=F)
    predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
    predicts$observed <- valids$label
    predicts$samp_id<-row.names(valids)
    pred_list1<-rbind(pred_list1,predicts)
    ROC <- roc(predicts$observed, as.numeric(predicts$`1`))
    auc <- as.numeric(ROC[["auc"]])
    acc <- sum(predicts$predicted==predicts$observed)
    accu1 <- rbind(accu1,c(prot_order[(i+1)],n,acc/length(fold),auc))
  }
}
accu1<-data.frame(accu1)
names(accu1)<-c("start_prot","fold","acc","auc")

write.csv(accu1,"20230220_8_prot_model_5_flod.csv",row.names = F)


###############################################
list_lab<-list()
for (i in 1:6) {
  list_lab[[i]]<-pred_list[((i-1)*168+1):(i*168),]
  
}

#1

ROC_1 <- roc(list_lab[[1]]$observed, as.numeric(list_lab[[1]]$`1`))
auc_1 <- as.numeric(ROC_1[["auc"]])
acc_1 <- sum(list_lab[[1]]$predicted==list_lab[[1]]$observed)/nrow(list_lab[[1]])

#2

ROC_2 <- roc(list_lab[[2]]$observed, as.numeric(list_lab[[2]]$`1`))
auc_2 <- as.numeric(ROC_2[["auc"]])
acc_2 <- sum(list_lab[[2]]$predicted==list_lab[[2]]$observed)/nrow(list_lab[[2]])

#3

ROC_3 <- roc(list_lab[[3]]$observed, as.numeric(list_lab[[3]]$`1`))
auc_3 <- as.numeric(ROC_3[["auc"]])
acc_3 <- sum(list_lab[[3]]$predicted==list_lab[[3]]$observed)/nrow(list_lab[[3]])

#4

ROC_4 <-roc(list_lab[[4]]$observed, as.numeric(list_lab[[4]]$`1`))
auc_4 <-as.numeric(ROC_4[["auc"]])
acc_4 <-sum(list_lab[[4]]$predicted==list_lab[[4]]$observed)/nrow(list_lab[[4]])

#5

ROC_5 <- roc(list_lab[[5]]$observed, as.numeric(list_lab[[5]]$`1`))
auc_5 <- as.numeric(ROC_5[["auc"]])
acc_5 <- sum(list_lab[[5]]$predicted==list_lab[[5]]$observed)/nrow(list_lab[[5]])

#6

ROC_6 <- roc(list_lab[[6]]$observed, as.numeric(list_lab[[6]]$`1`))
auc_6 <- as.numeric(ROC_6[["auc"]])
acc_6 <- sum(list_lab[[6]]$predicted==list_lab[[6]]$observed)/nrow(list_lab[[6]])

set.seed(123)
roc_temp[9,2]<-roc.test(ROC_all , ROC_1 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[10,2]<-roc.test(ROC_all , ROC_2 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[11,2]<-roc.test(ROC_all , ROC_3 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[12,2]<-roc.test(ROC_all , ROC_4 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[13,2]<-roc.test(ROC_all , ROC_5 , reuse.auc=FALSE)$p.value
set.seed(123)
roc_temp[14,2]<-roc.test(ROC_all , ROC_6 , reuse.auc=FALSE)$p.value
write.csv(roc_temp,"20230222_zzovc_14rf_model_auc_test.csv",row.names = F)


###########################################
#roc
pdf("20230221_zzovc_ROC_zzov_plasma_9model_total_auc.pdf", width = 8, height = 6)
plot(ROC_all, print.auc = TRUE, print.thres = F, xlab = "Specificity", ylab = "Sensitivity",
     main = "ROC Curve of 9 model", col = "red", 
     print.auc.x = 0.5, print.auc.y = 0.5)

plot(ROC_1, add = T, col = 'black', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_2, add = T, col = 'black', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_3, add = T, col = 'pink', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_4, add = T, col = 'green', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_5, add = T, col = 'orange', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)
plot(ROC_6, add = T, col = 'purple', print.auc = T, print.auc.x = 0.5, print.auc.y = 0.45)


legend("bottomright", legend = c(paste0("8 prot:",auc_all),paste0("CHI3L1:",auc_1),paste0("MUC1:",auc_2),paste0("SMPDL3B:",auc_3),
                                 paste0("SAA2:",auc_4),paste0("MXRA5:",auc_5),paste0("HPSE:",auc_6)),
       col=c("red", "blue",'black','green','orange','purple','darkblue',"pink"), lwd=2)
dev.off()










