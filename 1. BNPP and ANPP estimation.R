rm(list=ls())
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)


############################################################################################
###1 BNPP model
############################################################################################

#1.1 read data
bnpp=read.csv("./NPP and fBNPP Data.csv",header = T,as.is = T)%>%
  select(Code,Site,Longititude,Latitude,BNPP_unit_unified,Biomes,Order,
         BIO1,BIO2,BIO8,BIO12,BIO14,BIO15,BIO18,
         CECc,CFRAG,CLPC,CNrt,ELCO,ESP,GYPS,BULK,PHAQ,SDTO,TCEQ,MNPP)%>%
  mutate(Biomes=as.factor(Biomes),Order=factor(Order,levels=c("Alfisols","Andisols" ,"Aridisols","Entisols","Gelisols",
                                                              "Histosols","Ice/Glacier","Inceptisols","Mollisols",
                                                              "Oxisols","Spodosols","Ultisols","Vertisols","Water",
                                                              "Rocky Land","Shifting Sands","Salt"),
                                               labels = c(1:14,4,4,4)))%>%
  mutate(BNPP_unit_unified=log(BNPP_unit_unified*1000))  ##Mg ha-1 --> log(g m-2)

##remove outliers
bnpp_fine=bnpp%>%filter(BNPP_unit_unified>-10)%>%filter(!Code%in%c(6,7))


#1.2 datasplit
SplitFun=function(all.data_fine,train_ratio=0.8,seed=150){
  dataout=list()
  Sites=unique(all.data_fine$Code)
  SitesNO=length(Sites)
  set.seed(seed)
  train_ix=sample(1:SitesNO, SitesNO*train_ratio)
  train_Sites=Sites[train_ix]
  test_Sites=Sites[-train_ix]
  train_ix_alldata=which(all.data_fine$Code%in%train_Sites)
  dataout$train_data=all.data_fine[train_ix_alldata,]
  dataout$test_data=all.data_fine[-train_ix_alldata,]
  dataout$vncol=ncol(dataout$train_data)
  print(paste0("seed=",seed,", train_ratio=",round(nrow(dataout$train_data)/nrow(all.data_fine),2)))
  return(dataout)
}

# datasplit 80% for train, 20% for test
Datasplit=SplitFun(all.data_fine=bnpp_fine,train_ratio=0.8,seed=150)

# select predictors
input_ix=c(5,6:Datasplit$vncol)
colnames(Datasplit$train_data)[input_ix]

##1.3 Machine Learning model building
ModelML=function(Datasplit,input_ix,mod="RF"){
  Datasplit=Datasplit
  #.1 data split
  Datasplit$train_data=Datasplit$train_data[,input_ix]
  Datasplit$test_data=Datasplit$test_data[,input_ix]
  
  #.2 ML formula
  mod.formula <- as.formula(paste0(colnames(Datasplit$train_data[1]),"~."))
  
  #.3 Corss-Validation
  fitControl <- trainControl(method="repeatedcv", number=10, repeats = 10,allowParallel=T) 
  
  #.4 preset
  MethodList=data.frame(mod=c("RF","XG","CU","SVM","MS","BN","BR","BL","BG","RR","MNN"),
                        method=c("ranger","xgbTree","cubist","svmRadial","bagEarth",
                                 "brnn","bridge","blassoAveraged","bayesglm","foba","mlpWeightDecayML"),
                        Unique=1:11)
  TuneGrid=list(RF=expand.grid(mtry = c(2,4,6),splitrule = "variance",min.node.size=c(3,5)),
                XG=expand.grid(eta = c(0.3,0.5), nrounds = c(150), 
                               max_depth = c(3,5), gamma = 0, colsample_bytree = 0.8, 
                               min_child_weight = 1, subsample=1),
                CU=expand.grid(committees=c(10,20), 
                               neighbors=c(9,18)),   
                SVM=expand.grid(sigma=c(0.02,0.04), C=1),
                MS=expand.grid(degree = c(1,2),nprune = c(40,50)),
                BN=expand.grid(neurons=3))
  eval(parse(text=paste0("TuneGridselect=TuneGrid$",mod)))
  
  #.5 machine training
  if(mod=="MS"|mod=="SVM"|mod=="BN"){
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"),
                          tuneGrid = TuneGridselect,
                          metric = "RMSE")
  }else if(mod=="BR"|mod=="BL"|mod=="BG"|mod=="RR"|mod=="MNN"){
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"),
                          metric = "RMSE")
  }else{
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"), 
                          tuneGrid = TuneGridselect,
                          importance = "impurity",
                          metric = "RMSE")
  }
  
  
  #.6 predict
  Datasplit$train_data$Pred= predict(Datasplit$mod,newdata = Datasplit$train_data)
  Datasplit$test_data$Pred = predict(Datasplit$mod,newdata = Datasplit$test_data)
  Datasplit$train_data$Model=mod
  Datasplit$test_data$Model=mod
  
  #.7 evaluation
  Train_validation=postResample(pred =  Datasplit$train_data$Pred, obs = Datasplit$train_data[,1])
  Train_validation=data.frame(RMSE=Train_validation[1],R2=Train_validation[2]);Train_validation$Type="Train"
  Test_validation=postResample(pred = Datasplit$test_data$Pred, obs = Datasplit$test_data[,1])
  Test_validation=data.frame(RMSE=Test_validation[1],R2=Test_validation[2]);Test_validation$Type="Test"
  Datasplit$Evaluation=rbind(Train_validation,Test_validation);rownames(Datasplit$Evaluation)=1:nrow(Datasplit$Evaluation)
  
  
  #.8 importance
  Datasplit$MLImp=varImp(Datasplit$mod,scale=F)
  
  return(Datasplit)
}


## 1.4 Train model
TrainoutRF=ModelML(Datasplit,input_ix,mod="RF")


##1.5 compile results
##Prediction for trainning data
dataF=data.frame(bnpp=TrainoutRF$train_data$BNPP_unit_unified)
dataF$RF=TrainoutRF$train_data$Pred

##prediction for testing data
dataT=data.frame(bnpp=TrainoutRF$test_data$BNPP_unit_unified)
dataT$RF=TrainoutRF$test_data$Pred

dataF$Type="Calibration"
dataT$Type="Validation"

AllPre=rbind(dataF,dataT)
AllPre=AllPre%>%gather(key = Model,value = Pred,-bnpp,-Type)

RMSEFun=function(simvalue,obsvalue) {round(sqrt(mean((simvalue-obsvalue)^2)),2)}
R2Fun=function(simvalue,obsvalue){round(summary(lm(simvalue~1+obsvalue))$r.squared,2)}


Evaluation=AllPre%>%dplyr::group_by(Type,Model)%>%
  dplyr::summarise(RMSE=RMSEFun(Pred,bnpp),R2=R2Fun(Pred,bnpp))%>%filter(Type=="Validation")

##1.6 save results
saveRDS(Evaluation,"BNPP_Evaluation.rds")
saveRDS(AllPre,"BNPP_Prediction.rds")



############################################################################################
###2 ANPP model
############################################################################################

#2.1 read data
Anpp=read.csv("./fBNPP_all_Year_EVP_12.30.csv",header = T,as.is = T)%>%
  select(Code,Site,Longititude,Latitude,ANPP_unit_unified,Biomes,Order,
         BIO1,BIO2,BIO8,BIO12,BIO14,BIO15,BIO18,
         CECc,CLPC,CNrt,ELCO,ESP,GYPS,PHAQ,SDTO,TCEQ,MNPP)%>%
  mutate(Biomes=as.factor(Biomes),Order=factor(Order,levels=c("Alfisols","Andisols" ,"Aridisols","Entisols","Gelisols",
                                                              "Histosols","Ice/Glacier","Inceptisols","Mollisols",
                                                              "Oxisols","Spodosols","Ultisols","Vertisols","Water",
                                                              "Rocky Land","Shifting Sands","Salt"),
                                               labels = c(1:14,4,4,4)))%>%
  mutate(ANPP_unit_unified=log(ANPP_unit_unified*1000))  ##Mg ha-1 --> log(g m-2)



#2.2 datasplit
SplitFun=function(all.data_fine,train_ratio=0.8,seed=150){
  dataout=list()
  Sites=unique(all.data_fine$Code)
  SitesNO=length(Sites)
  set.seed(seed)
  train_ix=sample(1:SitesNO, SitesNO*train_ratio)
  train_Sites=Sites[train_ix]
  test_Sites=Sites[-train_ix]
  train_ix_alldata=which(all.data_fine$Code%in%train_Sites)
  dataout$train_data=all.data_fine[train_ix_alldata,]
  dataout$test_data=all.data_fine[-train_ix_alldata,]
  dataout$vncol=ncol(dataout$train_data)
  print(paste0("seed=",seed,", train_ratio=",round(nrow(dataout$train_data)/nrow(all.data_fine),2)))
  return(dataout)
}

# datasplit 80% for train, 20% for test
Datasplit=SplitFun(all.data_fine=Anpp,train_ratio=0.8,seed=140)

# select predictors
input_ix=c(5,6:Datasplit$vncol)
colnames(Datasplit$train_data)[input_ix]

##2.3 Machine Learning model building
ModelML=function(Datasplit,input_ix,mod="RF"){
  Datasplit=Datasplit
  #.1 data split
  Datasplit$train_data=Datasplit$train_data[,input_ix]
  Datasplit$test_data=Datasplit$test_data[,input_ix]
  
  #.2 ML formula
  mod.formula <- as.formula(paste0(colnames(Datasplit$train_data[1]),"~."))
  
  #.3 Corss-Validation
  fitControl <- trainControl(method="repeatedcv", number=10, repeats = 10,allowParallel=T) 
  
  #.4 preset
  MethodList=data.frame(mod=c("RF","XG","CU","SVM","MS","BN","BR","BL","BG","RR","MNN"),
                        method=c("ranger","xgbTree","cubist","svmRadial","bagEarth",
                                 "brnn","bridge","blassoAveraged","bayesglm","foba","mlpWeightDecayML"),
                        Unique=1:11)
  TuneGrid=list(RF=expand.grid(mtry = c(2,4,6),splitrule = "variance",min.node.size=c(3,5)),
                XG=expand.grid(eta = c(0.3,0.5), nrounds = c(150), 
                               max_depth = c(3,5), gamma = 0, colsample_bytree = 0.8, 
                               min_child_weight = 1, subsample=1),
                CU=expand.grid(committees=c(10,20), 
                               neighbors=c(9,18)),   
                SVM=expand.grid(sigma=c(0.02,0.04), C=1),
                MS=expand.grid(degree = c(1,2),nprune = c(40,50)),
                BN=expand.grid(neurons=3))
  eval(parse(text=paste0("TuneGridselect=TuneGrid$",mod)))
  
  #.5 machine training
  if(mod=="MS"|mod=="SVM"|mod=="BN"){
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"),
                          tuneGrid = TuneGridselect,
                          metric = "RMSE")
  }else if(mod=="BR"|mod=="BL"|mod=="BG"|mod=="RR"|mod=="MNN"){
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"),
                          metric = "RMSE")
  }else{
    Datasplit$mod = train(mod.formula, 
                          data=Datasplit$train_data,
                          method = MethodList[MethodList$mod==mod,"method"],
                          trControl = fitControl,
                          preProcess=c("center", "scale", "YeoJohnson", "nzv"), 
                          tuneGrid = TuneGridselect,
                          importance = "impurity",
                          metric = "RMSE")
  }
  
  
  #.6 predict
  Datasplit$train_data$Pred= predict(Datasplit$mod,newdata = Datasplit$train_data)
  Datasplit$test_data$Pred = predict(Datasplit$mod,newdata = Datasplit$test_data)
  Datasplit$train_data$Model=mod
  Datasplit$test_data$Model=mod
  
  #.7 evaluation
  Train_validation=postResample(pred =  Datasplit$train_data$Pred, obs = Datasplit$train_data[,1])
  Train_validation=data.frame(RMSE=Train_validation[1],R2=Train_validation[2]);Train_validation$Type="Train"
  Test_validation=postResample(pred = Datasplit$test_data$Pred, obs = Datasplit$test_data[,1])
  Test_validation=data.frame(RMSE=Test_validation[1],R2=Test_validation[2]);Test_validation$Type="Test"
  Datasplit$Evaluation=rbind(Train_validation,Test_validation);rownames(Datasplit$Evaluation)=1:nrow(Datasplit$Evaluation)
  
  
  #.8 importance
  Datasplit$MLImp=varImp(Datasplit$mod,scale=F)
  
  return(Datasplit)
}


## 2.4 Train model
TrainoutRF=ModelML(Datasplit,input_ix,mod="RF")


##2.5 compile results
##Prediction for trainning data
dataF=data.frame(ANPP=TrainoutRF$train_data$ANPP_unit_unified)
dataF$RF=TrainoutRF$train_data$Pred

##prediction for testing data
dataT=data.frame(ANPP=TrainoutRF$test_data$ANPP_unit_unified)
dataT$RF=TrainoutRF$test_data$Pred

dataF$Type="Calibration"
dataT$Type="Validation"

AllPre=rbind(dataF,dataT)
AllPre=AllPre%>%gather(key = Model,value = Pred,-ANPP,-Type)

RMSEFun=function(simvalue,obsvalue) {round(sqrt(mean((simvalue-obsvalue)^2)),2)}
R2Fun=function(simvalue,obsvalue){round(summary(lm(simvalue~1+obsvalue))$r.squared,2)}


Evaluation=AllPre%>%dplyr::group_by(Type,Model)%>%
  dplyr::summarise(RMSE=RMSEFun(Pred,ANPP),R2=R2Fun(Pred,ANPP))%>%filter(Type=="Validation")

##2.6 save results
saveRDS(Evaluation,"ANPP_Evaluation.rds")
saveRDS(AllPre,"ANPP_Prediction.rds")

