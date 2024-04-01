setwd("E:/Qsync/work/zzov/ZZOV/ZZOV/experiment20220120/matrix_20221013/xgboxt/plasma_MRM")
#
library("xgboost")
library("Matrix")


protmat<-read.csv("20221031zzova_plasma_MRM_51prot_matrix_delete3CiRT.csv",header = T,row.names = 1)


protmat<-protmat[,!grepl("*pool",colnames(protmat))]
# protmat<-log2(protmat)
protmat1<-data.frame(t(protmat))
library(xlsx)
pat_inf<-read.xlsx("20221116sample_info_final_20230201.xlsx",sheetIndex = 2)
samp_inf<-read.csv("sample_inf_plasma.csv",header = T)

protmat1$batch<-row.names(protmat1)
protmat1$batch<-gsub("train_","",protmat1$batch)
protmat1$batch<-gsub("test_","",protmat1$batch)
protmat1$pat_id<-samp_inf$patient_barcode[match(protmat1$batch,samp_inf$batchID)]
protmat1$in_test<-row.names(protmat1)
##################################
#bio rep corrlation
techRep<-unique(protmat1$pat_id[duplicated(protmat1$pat_id)])

cortechRep<-vector()
m<-1
for(i in 1:length(techRep)){
  indSel<-which(protmat1$pat_id==techRep[i])
  if(length(indSel)>0){
    matCor<-as.matrix(protmat1[indSel,1:(ncol(protmat1)-3)])
    comAll<-combn(nrow(matCor),2)
    for(j in 1:ncol(comAll)){
      #做相关性分???
      cortechRep[m]<-cor(matCor[comAll[1,j],],matCor[comAll[2,j],], use = "na.or.complete")
      label_name<-row.names(matCor)[comAll[,j]]
      names(cortechRep)[m]<-paste0(label_name,collapse = ";")
      m<-m+1
      
    }
  }else{
    print(i)
    
  }
}


library(vioplot)

allsample<-protmat1

# write.csv(allsample,"zzoc_unique_matrix_plasma_MRM.csv",row.names = T)
# uniq_pat<-pat_inf[!duplicated(pat_inf$patientID),]
# allsample$patientID<-row.names(allsample)

unisamp<-merge(allsample,pat_inf,by.x ="pat_id" ,by.y ="patientID" ,all.x = T)

#############################################################################
MRM_mat<-unisamp

vali_mat<-MRM_mat[grepl("vali",MRM_mat$train_test_vali),]

train_mat<-MRM_mat[!grepl("vali",MRM_mat$train_test_vali),]
#################################################################################################
#ttest
MRM_mat1<-log2(train_mat[,2:52])
readP<-MRM_mat1[grepl("1",train_mat$label),]

readR<-MRM_mat1[grepl("0",train_mat$label),]

FC<-pvalue<-vector()
for(i in 1:ncol(readP)){
  
  if(sum(!is.na(readP[,i]))>1&sum(!is.na(readR[,i]))>1){
    FC[i]<-mean(2^readR[,i],na.rm=T)/mean(2^readP[,i],na.rm = T)
    pvalue[i]<-t.test(readR[,i],readP[,i], paired = F,  var.equal = F)$p.value
  }else{
    FC[i]<-pvalue[i]<-NA
  }
  
}
pvalueAd<-p.adjust(pvalue,method = "BH",n=sum(!is.na(pvalue)))
FC_PAd<-data.frame(cbind(FC,pvalueAd,pvalue))
row.names(FC_PAd)<-colnames(readP)[1:ncol(readP)]

# write.csv(FC_PAd,"20230201_zzov_plasma_mrm_0_vs_1_t_test_FC_PAd.csv")


MRM_mat1$age<-train_mat$Age
# MRM_mat2<-MRM_mat1[,4:ncol(MRM_mat1)]
# INT function
int <- function(x){
  #inverse normal transformation
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

# Age residuals + standardization

res_stan <- function(column){
  res <- residuals(lm(formula = column~age,data=MRM_mat1, na.action=na.exclude)) %>% int()
  # res <- int(res)
  return(res)
}

library("survival")
library("survminer")

pro_res <- data.frame(sapply(MRM_mat1[,1:51], res_stan))
pro_res$time<-train_mat$cox_RFS
pro_res$disease<-train_mat$cox_recurrence
pro_res$disease<-gsub("Yes",1,pro_res$disease)
pro_res$disease<-gsub("No",0,pro_res$disease)
pro_res$disease<-gsub("NO",0,pro_res$disease)
pro_res$disease<-as.numeric(pro_res$disease)
pvalue_cox<-exp<-ecp_low<-exp_up<-vector()
for (i in 1:(ncol(pro_res)-2)) {
  pvalue_cox[i]<-summary(coxph(Surv(time, disease) ~ pro_res[,i], data = pro_res))$logtest[3]
  exp[i]<-summary(coxph(Surv(time, disease) ~ pro_res[,i], data = pro_res))$conf.int[1]
  ecp_low[i]<-summary(coxph(Surv(time, disease) ~ pro_res[,i], data = pro_res))$conf.int[3]
  exp_up[i]<-summary(coxph(Surv(time, disease) ~ pro_res[,i], data = pro_res))$conf.int[4]
  
}
pvalueAd<-p.adjust(pvalue_cox,method = "BH")
P_cox<-data.frame(cbind(pvalue_cox,pvalueAd,exp,ecp_low,exp_up))
row.names(P_cox)<-colnames(pro_res)[1:(ncol(pro_res)-2)]
# write.csv(P_cox,"20230201_zzov_plasma_mrm_cov_p.csv",row.names = T)

################################################
#
ttestp<-FC_PAd[FC_PAd$pvalue<0.05,]
coxp<-P_cox[P_cox$pvalue_cox<0.05,]
#合集
prot_list<-unique(c(row.names(ttestp),row.names(coxp)))
# write.csv(prot_list,"20230201_zzov_plasma_mrm_input_feature.csv",row.names = F)


######################################################################################################
train_prot<-train_mat[,colnames(train_mat)%in%prot_list]
train_prot<-cbind(train_prot,train_mat[,c(1,53:62)])
row.names(train_prot)<-train_prot$batch
# test_mat<-train_prot[grepl("test",train_prot$train_test_vali),-ncol(train_prot)]
# train_mat1<-train_prot[!grepl("test",train_prot$train_test_vali),-ncol(train_prot)]

temp_mean<-apply(train_prot[,1:34],2, mean, na.rm = TRUE)
temp_sd<-apply(train_prot[,1:34],2, sd, na.rm = TRUE)
train_prot[,1:34]<-scale(train_prot[,1:34],scale = T,center = T)

################
vali_mat1<-vali_mat[,colnames(vali_mat)%in%prot_list]
vali_mat2<-data.frame(matrix(NA,nrow = nrow(vali_mat1),ncol = ncol(vali_mat1)))
names(vali_mat2)<-names(vali_mat1)
# row.names(Temp1)<-row.names(Temp)
for (i in 1:length(temp_mean)) {
  aa<-vali_mat1[,i]
  vali_mat2[,i]<-sapply(aa,function(e){(e-temp_mean[i])/temp_sd[i]} )
  
}
vali_mat2<-cbind(vali_mat2,vali_mat[,c(1,53:62)])
row.names(vali_mat2)<-vali_mat2$batch


################

train_prot<-train_prot[,-c(36,37)]
train_prot<-train_prot[,c(1:34,36:43,35)]
# write.csv(train_prot,"20230201_zzov_palsma_mrm_ML_train_input.csv",row.names = T)

vali_mat2<-vali_mat2[,-c(36,37)]
vali_mat2<-vali_mat2[,c(1:34,36:43,35)]
# write.csv(vali_mat2,"20230201_zzov_plasma_mrm_ML_vali_input.csv",row.names = T)
# write.csv(vali_mat2,"20230306_zzov_plasma_mrm_ML_vali_input_tims1789.csv",row.names = T)
##################

library("xgboost")
library("Matrix")
library(caret)



# library(caret)
mat1<-train_prot[train_prot$label==0,]
mat2<-train_prot[train_prot$label==1,]
impt<-data.frame(matrix(1,nrow = 1,ncol = 8))
names(impt)<-c("Feature" , "Gain" , "Cover","Frequency","eta_s","ss" , "tims","seed"  )
set.seed(123)
seed <- runif(100,1,1000)
tims<-runif(10,1,1000)
# tim<-609.1262

for (tim in tims) {
  set.seed(tim)
  mat1_pat<-unique(mat1$pat_id)
  mat2_pat<-unique(mat2$pat_id)
  train_0<-sample(mat1_pat,48)
  train_1<-sample(mat2_pat,51)
  train_mat<-rbind(train_prot[train_prot$pat_id%in%train_0,-ncol(train_prot)],train_prot[train_prot$pat_id%in%train_1,-ncol(train_prot)])
  test_mat<-train_prot[!row.names(train_prot)%in%row.names(train_mat),-ncol(train_prot)]
  

for (eta_s in seq(0.1,0.3,0.04)) {
  for (ss in seq(0.5,1,0.05)) {
  # eta_s<-0.3
  # ss=0.5
    
    for (seed1 in seed) {
      # seed1=123
         set.seed(seed1)
      samp<-sample(nrow(train_mat),ceiling(nrow(train_mat)*0.6))
          
              trains <- train_mat[samp,]
              trains$label <- train_mat$label[samp]

              train_matrix <- sparse.model.matrix(label ~ ., data = trains)
              
              train_label <- as.numeric(trains$label)
              
              train_fin <- list(data=train_matrix,label=train_label) 
              
              dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label) 
              
              
              xgb <- xgboost(data = dtrain, eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
                           set.seed(1234),eval_metric='logloss')
              
              
              #
              importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb) 
              importance$eta_s<-eta_s
              importance$ss<-ss
              importance$tims<-tim
              importance$seed<-seed1
              impt<-rbind(impt,importance)
    }       
  }

 }
}

impt<-impt[-1,]

table_temp<-data.frame(table(impt$Feature,impt$eta_s,impt$ss,impt$tims))
names(table_temp)<-c("Feature","eta_","ss","tims","Freq")
write.csv(table_temp,"20240130_zzov_plasma_MRM_xdboost_step1_table.csv",row.names = F)
library(plyr)
library(magrittr)
library(rlist)
impt<-impt %>% mutate(prot = unclass(factor(Feature)))
impt1<-impt[,c(1,9)]
impt2<-impt[,-1]
mean_temp <- aggregate(impt2,by=list(impt2$tims,impt2$eta_s,impt2$ss,impt2$prot),mean)


impt1<-impt1[!duplicated(impt1$Feature),]
mean_temp$feature<-impt1$Feature[match(mean_temp$prot,impt1$prot)]
mean_temp$type<-paste(mean_temp$eta_s,mean_temp$ss,mean_temp$tims,sep = "_")
write.csv(mean_temp,"20240130_zzov_plasma_MRM_xdboost_step1_mean.csv",row.names = F)


tmp3<-unique(mean_temp$type)

##################################################################################################
accu3<-c()
for (i in 1:length(tmp3)) {
  # i=1
  indsel<-mean_temp[mean_temp$type==tmp3[i],]
  set.seed(indsel$tims[1])
  mat1_pat<-unique(mat1$pat_id)
  mat2_pat<-unique(mat2$pat_id)
  train_0<-sample(mat1_pat,48)
  train_1<-sample(mat2_pat,51)
  train_mat<-rbind(train_prot[train_prot$pat_id%in%train_0,-ncol(train_prot)],train_prot[train_prot$pat_id%in%train_1,-ncol(train_prot)])
  test_mat<-train_prot[!row.names(train_prot)%in%row.names(train_mat),-ncol(train_prot)]
  
  eta_s<-indsel$eta_s[1]
  ss<-indsel$ss[1]
  for(j in 5:15){
    # j=10
  feature_choose<-indsel$feature[order(indsel$Frequency,decreasing=T)][1:j]
  
  # feature_choose<-importance$Feature[importance$Gain>fea_selt]
  train_mat2<-train_mat[,colnames(train_mat)%in%feature_choose]
  train_mat2$label<-train_mat$label
  test2<-test_mat[,colnames(test_mat)%in%feature_choose]
  test2$label<-test_mat$label
  vali2<-vali_mat2[,colnames(vali_mat2)%in%feature_choose]
  vali2$label<-vali_mat2$label
  
  
  train_matrix2 <- sparse.model.matrix(label ~ ., data = train_mat2)
  test_matrix2 <- sparse.model.matrix(label ~ ., data =test2)
  vali_matrix2 <- sparse.model.matrix(label ~ ., data =vali2)
  
  train_label2 <- as.numeric(train_mat2$label)
  test_label2 <-  as.numeric(test2$label)
  vali_label2 <-  as.numeric(vali2$label)
  
  train_fin2 <- list(data=train_matrix2,label=train_label2) 
  test_fin2 <- list(data=test_matrix2,label=test_label2) 
  vali_fin2 <- list(data=vali_matrix2,label=vali_label2) 
  
  dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
  dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)
  dvali2 <- xgb.DMatrix(data = vali_fin2$data, label = vali_fin2$label)
  
  
  xgb1 <- xgboost(data = dtrain2, eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
              
                  scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
  
  
  
    
    
    ########################################
    #train
   
    
    pre_xgb = round(predict(xgb1,newdata = dtrain2))
    pre_xgb_train<-predict(xgb1,newdata = dtrain2)
    aaa<-data.frame(table(train_label2,pre_xgb,dnn=c("true","pre")))
    acc_train<-sum(pre_xgb==train_label2)/length(train_label2)
    
    #ROC
    library(pROC)
    xgboost_roc_train <- roc(train_label2,as.numeric(pre_xgb_train))
    
    
    auc_train<-xgboost_roc_train[["auc"]]
    ##################################################
    
    ########################################
    #test
    
    
    pre_xgb = round(predict(xgb1,newdata = dtest2))
    pre_xgb_test = predict(xgb1,newdata = dtest2)
    aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
    acc_test<-sum(pre_xgb==test_label2)/length(test_label2)
    
    #ROC曲线
    library(pROC)
    xgboost_roc_test <- roc(test_label2,as.numeric(pre_xgb_test))
    
    
    auc_test<-xgboost_roc_test[["auc"]]
    ##################################################
    #vali
   
    
    pre_xgb = round(predict(xgb1,newdata =dvali2))
    pre_xgb_vali = predict(xgb1,newdata =dvali2)
    aaa_vali<-data.frame(table(vali_label2,pre_xgb_vali,dnn=c("true","pre")))
    acc_vali<-sum(pre_xgb==vali_label2)/length(vali_label2)
    
    #ROC曲线
    library(pROC)
    xgboost_roc_vali <- roc(vali_label2,as.numeric(pre_xgb_vali))
    
    
    auc_vali<-xgboost_roc_vali[["auc"]]
    
    accu3 <- rbind(accu3,c(indsel$tims[1],eta_s,ss,j,acc_train,auc_train,acc_test,auc_test,acc_vali,auc_vali))
    
  }
}

accu3<-data.frame(accu3)
colnames(accu3)<-c("tims","eta_s","ss","fea_selt","acc_train","auc_train","acc_test","auc_test","acc_vali","auc_vali")

write.csv(accu3,"20230201_zzova_xgboost_plasma_MRM_step2_mean.csv",row.names = F)

#############################################################

mean_temp<-read.csv("20230201_zzov_plasma_MRM_xdboost_step1_mean.csv",header = T)
accu3<-read.csv("20230201_zzova_xgboost_plasma_MRM_step2_mean.csv",header = T)

accu4<-accu3[accu3$acc_test>0.7,]
a<-unique(accu4$tims)
accu4_1<-accu4[(accu4$tims==a[1]|accu4$tims==a[2]),]
accu4_2<-accu4[(accu4$tims==a[4]),]
accu4_3<-accu4[(accu4$tims==a[3]|accu4$tims==a[5]|accu4$tims==a[6])|(accu4$tims==a[7]|accu4$tims==a[8]),]

accu5<-c()


for (i in 1:nrow(accu4_2)) {
    for (gam in seq(0.05,0.2,0.05)) {
      for (max_dep in seq(3,10,1)) {
        for (colsamp in seq(0.1,1,0.1)) {
          for(min_chid in seq(1,5,1)) {
            
            indsel<-mean_temp[mean_temp$type==paste(accu4_2$eta_s[i],accu4_2$ss[i],accu4_2$tims[i],sep = "_"),]
            
            set.seed(indsel$tims[1])
            mat1_pat<-unique(mat1$pat_id)
            mat2_pat<-unique(mat2$pat_id)
            train_0<-sample(mat1_pat,48)
            train_1<-sample(mat2_pat,51)
            train_mat<-rbind(train_prot[train_prot$pat_id%in%train_0,-ncol(train_prot)],train_prot[train_prot$pat_id%in%train_1,-ncol(train_prot)])
            test_mat<-train_prot[!row.names(train_prot)%in%row.names(train_mat),-ncol(train_prot)]
            
            eta_s<-accu4_2$eta_s[i]
            ss<-accu4_2$ss[i]
            feature_choose<-indsel$feature[order(indsel$Frequency,decreasing=T)][1:accu4_2$fea_selt[i]]
              # write.csv(feature_choose,"20221123_zzov_xgboost_best_search_plasma_MRM_top3_model_feature_final.csv",row.names = F)
              # feature_choose<-importance$Feature[importance$Gain>fea_selt]
              train_mat2<-train_mat[,colnames(train_mat)%in%feature_choose]
              train_mat2$label<-train_mat$label
              test2<-test_mat[,colnames(test_mat)%in%feature_choose]
              test2$label<-test_mat$label
              vali2<-vali_mat2[,colnames(vali_mat2)%in%feature_choose]
              vali2$label<-vali_mat2$label
              
              
              train_matrix2 <- sparse.model.matrix(label ~ ., data = train_mat2)
              test_matrix2 <- sparse.model.matrix(label ~ ., data =test2)
              vali_matrix2 <- sparse.model.matrix(label ~ ., data =vali2)
              
              train_label2 <- as.numeric(train_mat2$label)
              test_label2 <-  as.numeric(test2$label)
              vali_label2 <-  as.numeric(vali2$label)
              
              train_fin2 <- list(data=train_matrix2,label=train_label2) 
              test_fin2 <- list(data=test_matrix2,label=test_label2) 
              vali_fin2 <- list(data=vali_matrix2,label=vali_label2) 
              
              dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
              dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)
              dvali2 <- xgb.DMatrix(data = vali_fin2$data, label = vali_fin2$label)
              
              
              xgb2 <- xgboost(data = dtrain2, eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
                              gamma = gam,max_depth =max_dep,colsample_bytree =colsamp,min_child_weight =min_chid,
                              scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
              
              
              
              
              
              ########################################
              #train
              
              
              pre_xgb = round(predict(xgb2,newdata = dtrain2))
              pre_xgb_train = predict(xgb2,newdata = dtrain2)
              aaa<-data.frame(table(train_label2,pre_xgb,dnn=c("true","pre")))
              acc_train<-sum(pre_xgb==train_label2)/length(train_label2)
              
              #ROC
              library(pROC)
              xgboost_roc_train <- roc(train_label2,as.numeric(pre_xgb_train))
              
              
              auc_train<-xgboost_roc_train[["auc"]]
              ##################################################
              
              ########################################
              #test
             
              
              pre_xgb = round(predict(xgb2,newdata = dtest2))
              pre_xgb_test = predict(xgb2,newdata = dtest2)
              aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
              acc_test<-sum(pre_xgb==test_label2)/length(test_label2)
              
              #ROC曲线
              library(pROC)
              xgboost_roc_test <- roc(test_label2,as.numeric(pre_xgb_test))
              
             
              auc_test<-xgboost_roc_test[["auc"]]
              ##################################################
              #vali
             
              
              pre_xgb = round(predict(xgb2,newdata =dvali2))
              pre_xgb_vali = predict(xgb2,newdata =dvali2)
              aaa_vali<-data.frame(table(vali_label2,pre_xgb_vali,dnn=c("true","pre")))
              acc_vali<-sum(pre_xgb==vali_label2)/length(vali_label2)
              
              #ROC
              library(pROC)
              xgboost_roc_vali <- roc(vali_label2,as.numeric(pre_xgb_vali))
              
             
              auc_vali<-xgboost_roc_vali[["auc"]]
              
              accu5 <- rbind(accu5,c(row.names(accu4_2)[i],accu4_2$tims[i],eta_s,ss,gam,max_dep,colsamp,min_chid,accu4_2$fea_selt[i],acc_train,auc_train,acc_test,auc_test,acc_vali,auc_vali))
              }
            }
          }
          
        }
        
      }
      
    
    
  
  



accu5<-data.frame(accu5)
colnames(accu5)<-c("model","tims","eta_s","ss",'gam','max_dep','colsamp','min_chid',"fea_selt","acc_train","auc_train","acc_test","auc_test","acc_vali","auc_vali")

write.csv(accu5,"20230207_zzov_xgboost_plasma_MRM_step3_3_2.csv",row.names = F)




##########################################################

 mean_temp<-read.csv("20230201_zzov_plasma_MRM_xdboost_step1_mean.csv",header = T)
 accu3<-read.csv("20230201_zzova_xgboost_plasma_MRM_step2_mean.csv",header = T)
 unique(accu3$tims)[3]
 #1789
 eta_s=0.26;ss=0.6;gam=0.2;max_dep=5;colsamp=0.7;min_chid=1 ;fea_selt=11
          
         indsel<-mean_temp[mean_temp$type==paste(accu3$eta_s[1789],accu3$ss[1789],accu3$tims[1789],sep = "_"),]
          
          set.seed(unique(accu3$tims)[3])
          mat1_pat<-unique(mat1$pat_id)
          mat2_pat<-unique(mat2$pat_id)
          train_0<-sample(mat1_pat,48)
          train_1<-sample(mat2_pat,51)
          train_mat<-rbind(train_prot[train_prot$pat_id%in%train_0,-ncol(train_prot)],train_prot[train_prot$pat_id%in%train_1,-ncol(train_prot)])
          test_mat<-train_prot[!row.names(train_prot)%in%row.names(train_mat),-ncol(train_prot)]
          
          feature_choose<-indsel$feature[order(indsel$Frequency,decreasing=T)][1:fea_selt]
          write.csv(feature_choose,paste0("20230216_zzov_xgboost_plasma_MRM_feature_accu3_mode1789.csv"),row.names = F)
          # feature_choose<-importance$Feature[importance$Gain>fea_selt]
          
          train_mat2<-train_mat[,colnames(train_mat)%in%feature_choose]
          train_mat2$label<-train_mat$label
          test2<-test_mat[,colnames(test_mat)%in%feature_choose]
          test2$label<-test_mat$label
          vali2<-vali_mat2[,colnames(vali_mat2)%in%feature_choose]
          vali2$label<-vali_mat2$label
          
          
          train_matrix2 <- sparse.model.matrix(label ~ ., data = train_mat2)
          test_matrix2 <- sparse.model.matrix(label ~ ., data =test2)
          vali_matrix2 <- sparse.model.matrix(label ~ ., data =vali2)
          
          train_label2 <- as.numeric(train_mat2$label)
          test_label2 <-  as.numeric(test2$label)
          vali_label2 <-  as.numeric(vali2$label)
          
          train_fin2 <- list(data=train_matrix2,label=train_label2) 
          test_fin2 <- list(data=test_matrix2,label=test_label2) 
          vali_fin2 <- list(data=vali_matrix2,label=vali_label2) 
          
          dtrain2 <- xgb.DMatrix(data = train_fin2$data, label = train_fin2$label) 
          dtest2 <- xgb.DMatrix(data = test_fin2$data, label = test_fin2$label)
          dvali2 <- xgb.DMatrix(data = vali_fin2$data, label = vali_fin2$label)
          
          
          xgb2 <- xgboost(data = dtrain2, eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
                          gamma = gam,max_depth =max_dep,colsample_bytree =colsamp,min_child_weight =min_chid,
                          scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
          
          
         
          
          ########################################
          #train
         
          
          pre_xgb = round(predict(xgb2,newdata = dtrain2))
          pre_xgb_train<-predict(xgb2,newdata = dtrain2)
          aaa<-data.frame(table(train_label2,pre_xgb,dnn=c("true","pre")))
          acc_train<-sum(pre_xgb==train_label2)/length(train_label2)
          
          #ROC
          library(pROC)
          xgboost_roc_train <- roc(train_label2,as.numeric(pre_xgb_train))
          
          
          auc_train<-xgboost_roc_train[["auc"]]
          
          #predict 
          
          predicted1_1<-cbind(train_mat2,pre_xgb_train )
          write.csv(predicted1_1,"20230509_plasma_mrm_train_pred_mode1789.csv",row.names = T)
          predicted1_1$mean<-apply(predicted1_1[,1:11],1,mean)
          
          
          
          #
          library(plyr)
          library(magrittr)
          library(rlist)
          library('SHAPforxgboost')
          train_mat3<-train_mat[,colnames(train_mat)%in%feature_choose]
          train_label<-train_mat$label
          train_mat3<-as.matrix(train_mat3)
          xgb3 <- xgboost(data = train_mat3, label = train_label,eta=eta_s,objective='binary:logistic', nround=50,subsample=ss,
                          gamma = gam,max_depth =max_dep,colsample_bytree =colsamp,min_child_weight =min_chid,
                          scale_pos_weight=1,set.seed(1234),eval_metric='logloss')
          X2<-as.matrix(train_mat2[,-ncol(train_mat2)])
          shap_values <- shap.prep(xgb_model = xgb3, X_train = X2)
          
          pdf("20230216_zzov_xgboost_plasma_MRM_train_feature_final_accu3_mode1789.pdf",width = 8,height = 8)
          shap.plot.summary(shap_values,dilute = 1)
          dev.off()
          
          #kmplot
          df.score<-data.frame(row.names(train_mat2),as.numeric(pre_xgb_train))
          names(df.score)<-c("Sample_name","score")
          df.score$label<-MRM_mat$label[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_recurrence<-MRM_mat$cox_recurrence[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_RFS<-MRM_mat$cox_RFS[match(df.score$Sample_name,MRM_mat$batch)]
          
          df.score$score<-ifelse(df.score$score<0.5,"0","1")
          df.score$score<-as.numeric(df.score$score)
          write.csv(df.score,"20230216_zzov_xgboost_plasma_MRM_train_feature_final_pred_accu3_mode1789.csv",row.names = T)
          df.score$cox_recurrence<-ifelse(df.score$cox_recurrence=="Yes","1","0")
          df.score$cox_recurrence<-as.numeric(df.score$cox_recurrence)
          tmp.cox <-  survfit(Surv(cox_RFS,cox_recurrence)~score, data=df.score)
          
          pdf("20230216_zzov_xgboost_plasma_MRM_train_feature_final_kmplot_accu3_mode1789.pdf", width = 8, height = 6)
          ggsurvplot(tmp.cox,
                     pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
          )
          dev.off()
          
          
          ##################################################
          
          ########################################
          #test
          #混淆矩阵
          
          pre_xgb = round(predict(xgb2,newdata = dtest2))
          pre_xgb_test = predict(xgb2,newdata = dtest2)
          aaa<-data.frame(table(test_label2,pre_xgb,dnn=c("true","pre")))
          acc_test<-sum(pre_xgb==test_label2)/length(test_label2)
          
          #ROC曲线
          library(pROC)
          xgboost_roc_test <- roc(test_label2,as.numeric(pre_xgb_test))
          
          auc_test<-xgboost_roc_test[["auc"]]
          
          
          
          #蜜蜂图
          library(plyr)
          library(magrittr)
          library(rlist)
          library('SHAPforxgboost')
          
          X2<-as.matrix(test2[,-ncol(test2)])
          shap_values <- shap.prep(xgb_model = xgb3, X_train = X2)
          
          pdf("20230216_zzov_xgboost_plasma_MRM_test_feature_final_accu3_mode1789.pdf",width = 8,height = 8)
          shap.plot.summary(shap_values,dilute = 1)
          dev.off()
          
          #kmplot
          df.score<-data.frame(row.names(test2),as.numeric(pre_xgb_test))
          names(df.score)<-c("Sample_name","score")
          df.score$label<-MRM_mat$label[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_recurrence<-MRM_mat$cox_recurrence[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_RFS<-MRM_mat$cox_RFS[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$score<-ifelse(df.score$score<0.5,"0","1")
          df.score$score<-as.numeric(df.score$score)
          write.csv(df.score,"20230216_zzov_xgboost_plasma_MRM_test_feature_final_pred_accu3_mode1789.csv",row.names = T)
          df.score$cox_recurrence<-ifelse(df.score$cox_recurrence=="Yes","1","0")
          df.score$cox_recurrence<-as.numeric(df.score$cox_recurrence)
          tmp.cox <-  survfit(Surv(cox_RFS,cox_recurrence)~score, data=df.score)
          
          pdf("20230216_zzov_xgboost_plasma_MRM_test_feature_final_kmplot_accu3_mode1789.pdf", width = 8, height = 6)
          ggsurvplot(tmp.cox,
                     pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
          )
          dev.off()
          
          ##################################################
          #vali
         
          
          pre_xgb = round(predict(xgb2,newdata =dvali2))
          pre_xgb_vali = predict(xgb2,newdata =dvali2)
          aaa_vali<-data.frame(table(vali_label2,pre_xgb_vali,dnn=c("true","pre")))
          acc_vali<-sum(pre_xgb==vali_label2)/length(vali_label2)
          
          #ROC
          library(pROC)
          xgboost_roc_vali <- roc(vali_label2,as.numeric(pre_xgb_vali))
         res<- roc(vali_label2,as.numeric(pre_xgb_vali),
              smooth=TRUE,)
         res$auc
        
          auc_vali<-xgboost_roc_vali[["auc"]]
          
          
          
          #蜜蜂图
          library(plyr)
          library(magrittr)
          library(rlist)
          library('SHAPforxgboost')
          
          X2<-as.matrix(vali2[,-ncol(vali2)])
          shap_values <- shap.prep(xgb_model = xgb3, X_train = X2)
          
          pdf("20230407_zzov_xgboost_plasma_MRM_vali_feature_final_accu3_mode1789.pdf",width = 5,height = 5)
          shap.plot.summary(shap_values,dilute = 1)
          dev.off()
          
          
          #kmplot
          df.score<-data.frame(row.names(vali2),as.numeric(pre_xgb_vali))
          names(df.score)<-c("Sample_name","score")
          df.score$label<-MRM_mat$label[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_recurrence<-MRM_mat$cox_recurrence[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$cox_RFS<-MRM_mat$cox_RFS[match(df.score$Sample_name,MRM_mat$batch)]
          df.score$score<-ifelse(df.score$score<0.5,"0","1")
          df.score$score<-as.numeric(df.score$score)
          write.csv(df.score,"20230216_zzov_xgboost_plasma_MRM_vali_feature_final_pred_accu3_mode1789.csv",row.names = T)
          df.score$cox_recurrence<-ifelse(df.score$cox_recurrence=="Yes","1","0")
          df.score$cox_recurrence<-as.numeric(df.score$cox_recurrence)
          tmp.cox <-  survfit(Surv(cox_RFS,cox_recurrence)~score, data=df.score)
          
          pdf("20230216_zzov_xgboost_plasma_MRM_vali_feature_final_kmplot_accu3_mode1789.pdf", width = 8, height = 6)
          ggsurvplot(tmp.cox,
                     pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF")
          )
          dev.off()
          
          
          #####################
          

  
  
  
  
  