rm(list=ls())
library(MASS)
library("raster")
library(skimr)
library(rgdal)


###############################################################
##01 Cage estimation
###############################################################

### Obtain 14C data
###
#install.packages("ISRaD")
# library("ISRaD")
# setwd("E:/Wanggc/RadioCarbon")
# getwd()
# # set 'directory' argument to local path
# ISRaD.getdata(directory = "~", dataset= "full", extra = "TRUE", force_download = "TRUE")
# ISRaD.getdata(directory = "~", dataset = "full", force_download = "TRUE")
# #install.packages("soilcarbon")
# #library("soilcarbon")
# #radiocarbondata = soilcarbon_database
# #write.csv(radiocarbondata,file="radiocarbondata1.csv",row.names=F)
# 
# ###
# ### age curves
# ###
# AtmC14 = read.csv("./AtmC14_1.csv",header = T,as.is = T)
# AtmC14
# k_all = 1/seq(1,50000,by=1) ## Possabile k
# Ainput_t_N = AtmC14$deltaC14[AtmC14$Hemisphere=="North"] ## North hemisphere C14 in carbon input
# Ainput_t_S = AtmC14$deltaC14[AtmC14$Hemisphere=="South"] ## South hemisphere C14 in carbon input
# Years_sim_N = AtmC14$Year[AtmC14$Hemisphere=="North"]
# Years_sim_S = AtmC14$Year[AtmC14$Hemisphere=="South"]
# 
# Asoil_t_N = data.frame(array(NA,dim=c(length(Years_sim_N),length(k_all))))               ## North hemisphere C14 in soil
# Asoil_t_S = data.frame(array(NA,dim=c(length(Years_sim_S),length(k_all))))               ## South hemisphere C14 in soil
# colnames(Asoil_t_N) = 1:50000
# colnames(Asoil_t_S) = 1:50000
# rownames(Asoil_t_N) = Years_sim_N
# rownames(Asoil_t_S) = Years_sim_S
# dim(Asoil_t_N)
# dim(Asoil_t_S)
# 
# 
# ### North Hemisphere
# for (ii in 1:length(Years_sim_N)) {
#   if (ii == 1) {
#       Asoil_t_N[ii,] = 1000*(k_all/(k_all+1/8276)*(Ainput_t_N[ii]/1000+1) - 1)
#       
#     } else {
#       Asoil_t_N[ii,] = 1000*(k_all*(Ainput_t_N[ii]/1000+1)+(1-k_all-1/8276)*(Asoil_t_N[ii-1,]/1000+1) - 1)
#     }
# }
# save(Asoil_t_N,file="Asoil_t_N.Rdata")
# #windows()
# #par(mar=c(4,4,1,1))
# #filled.contour(x=Years_sim_N,y=log(1/k_all),z=as.matrix(Asoil_t_N))
# 
# ### South Hemisphere
# for (ii in 1:length(Years_sim_S)) {
#   if (ii == 1) {
#     Asoil_t_S[ii,] = 1000*(k_all/(k_all+1/8276)*(Ainput_t_S[ii]/1000+1) - 1)
#     
#   } else {
#     Asoil_t_S[ii,] = 1000*(k_all*(Ainput_t_S[ii]/1000+1)+(1-k_all-1/8276)*(Asoil_t_S[ii-1,]/1000+1) - 1)
#   }
# }
# save(Asoil_t_S,file="Asoil_t_S.Rdata")
# #filled.contour(x=Years_sim_S,y=rev(log(k_all)),z=as.matrix(Asoil_t_N[,2000:1]),log="y")
# #Asoil_t_S[,1:3]

###
### TAU: 14C-data-estimated 
###

radiocarbondata = read.csv("./ISRaD_database_files/ISRaD_extra_flat_layer_v 1.8.9.2021-04-13.csv",header = T,as.is = T)#,colClasses = "character")
names(radiocarbondata)
radiocarbondata$lat <- as.numeric(radiocarbondata$pro_lat)
radiocarbondata$X14c_sigma <- as.numeric(radiocarbondata$lyr_14c_sigma)
radiocarbondata$SOCD <- as.numeric(radiocarbondata$lyr_soc_fill_extra)*100 # g cm-2 to Mg ha-1
radiocarbondata=radiocarbondata[is.na(radiocarbondata$lat)==FALSE&is.na(radiocarbondata$lyr_14c)==FALSE&is.na(radiocarbondata$lyr_obs_date_y)==FALSE&radiocarbondata$lyr_obs_date_y < 2014&radiocarbondata$lyr_top > 0,] #atmospheric data is only available before 2014
length(radiocarbondata[,1])
Xsite =  radiocarbondata$site_name
Xprofile = radiocarbondata$pro_name
X14C = as.numeric(radiocarbondata$lyr_14c)
X14C_SD = as.numeric(radiocarbondata$X14c_sigma)
Xlayer_top = as.numeric(radiocarbondata$lyr_top)
Xlayer_bot = as.numeric(radiocarbondata$lyr_bot)
XYear = as.numeric(as.character(radiocarbondata$lyr_obs_date_y))
XLat = as.numeric(radiocarbondata$lat)
XLong = as.numeric(radiocarbondata$pro_long)
XOC = as.numeric(radiocarbondata$SOCD)

new.radiocarbondata <- cbind(Xsite,Xprofile,XLat,XLong,X14C,X14C_SD,Xlayer_top,Xlayer_bot,XYear,XOC)
summary(Xlayer_top)



new.radiocarbondata <- read.csv("./new.radiocarbondata.csv")
############################
#data pre-processing
############################
mean_sigma <- as.data.frame(cbind(new.radiocarbondata$X14C,new.radiocarbondata$X14C_SD))
#mean_sigma <- df[is.finite(rowSums(mean_sigma)),]
mean_sigma <- mean_sigma[complete.cases(mean_sigma),]#;length(mean_sigma[,1])
y=mean_sigma[,2];x=mean_sigma[,1]#;summary(mean_sigma)
fit_mean_sigma <- lm(y~x);summary(fit_mean_sigma)
#predict(fit_mean_sigma,newdata=data.frame(x=-200))

#sigma filling
new.radiocarbondata[is.na(new.radiocarbondata$X14C_SD)==TRUE,"X14C_SD"] <- predict(fit_mean_sigma,newdata=data.frame(x=as.numeric(new.radiocarbondata[is.na(new.radiocarbondata$X14C_SD)==TRUE,"X14C"])))

summary(cbind(new.radiocarbondata$X14C,new.radiocarbondata$X14C_SD))
new.radiocarbondata <- new.radiocarbondata[new.radiocarbondata$X14C < 1000,]

########################################
names(new.radiocarbondata)
Xsite =  new.radiocarbondata$Xsite
Xprofile = new.radiocarbondata$Xprofile
X14C = as.numeric(new.radiocarbondata$X14C)
X14C_SD = as.numeric(new.radiocarbondata$X14C_SD)
Xlayer_top = as.numeric(new.radiocarbondata$Xlayer_top)
Xlayer_bot = as.numeric(new.radiocarbondata$Xlayer_bot)
XYear = as.numeric(as.character(new.radiocarbondata$XYear))
XLat = as.numeric(new.radiocarbondata$XLat)
XLong = as.numeric(new.radiocarbondata$XLong)
XOC = as.numeric(new.radiocarbondata$XOC)
#####################################################
## Identify the tau in the possible 100000 values
Asoil_t_N<-get(load("Asoil_t_N.Rdata"))
Asoil_t_S<-get(load("Asoil_t_S.Rdata"))
AtmC14 = read.csv("./AtmC14.csv",header = T,as.is = T)
AtmC14
Ainput_t_N = AtmC14$deltaC14[AtmC14$Hemisphere=="North"] ## North hemisphere C14 in carbon input
Ainput_t_S = AtmC14$deltaC14[AtmC14$Hemisphere=="South"] ## South hemisphere C14 in carbon input
Years_sim_N = AtmC14$Year[AtmC14$Hemisphere=="North"]
Years_sim_S = AtmC14$Year[AtmC14$Hemisphere=="South"]

tau.finder <- function(x) {
  
  if (XLat[x]>=0) {
    iDATA = as.numeric(Asoil_t_N[which(Years_sim_N==max(c(1950,XYear[x]))),])
  } else {
    iDATA = as.numeric(Asoil_t_S[which(Years_sim_S==max(c(1950,XYear[x]))),])
  }
  
  tau.mean = which.min((iDATA-X14C[x])^2)
  tau.low = which.min((iDATA-X14C[x]-1.96*X14C_SD[x])^2)
  tau.up = which.min((iDATA-X14C[x]+1.96*X14C_SD[x])^2)
  return(c(tau.mean,tau.up,tau.low))
}

All.TAU = data.frame(array(NA,dim=c(length(XLong),3)))
for (ii in 1:length(XLong)) {
  All.TAU[ii,] = tau.finder(ii)
  print(c(ii,All.TAU[ii,],All.TAU[ii,2]-All.TAU[ii,1]))
}

All.TAU = t(sapply(1:length(XLong),tau.finder))
save(All.TAU,file="All.TAU.Rdata")

All.TAU1 = data.frame(array(NA,dim=c(length(XLong),3)))

for (i in c(1:length(XLong))) {
  All.TAU1[i,] <- unlist(All.TAU[i,])
  #print(c(i,unlist(All.TAU[,i])))
}
colnames(All.TAU1) = c("tau_mean","tau_up","tau_low")
summary(All.TAU1)
save(All.TAU1,file="All.TAU.DataFrame.Rdata")

All.TAU1 <- get(load("All.TAU.DataFrame.Rdata"))
Tau_14C = cbind(Xsite,Xprofile,XLat,XLong,XYear,Xlayer_top,Xlayer_bot,X14C,X14C_SD,XOC,XBD,XCorase,All.TAU1)
save(Tau_14C,file="Tau_14C.Rdata")


Tau_14C <- get(load("Tau_14C.Rdata"))
summary(Tau_14C)
length(All.TAU1[,1])


########################################################
##02 obtain environmental variables
############################################################
#########
############ Complete XOC, XBD and XCorase, and add other soil variables using WISE30 data
#########
XY <- data.frame(Tau_14C$XLong,Tau_14C$XLat) #in extracting WISE data, must be long,lat format
Uni_Pro <- match(unique(Tau_14C$Xprofile),Tau_14C$Xprofile)
Long_Lat <- data.frame(Tau_14C$Xsite[Uni_Pro],Tau_14C$Xprofile[Uni_Pro],Tau_14C$XLong[Uni_Pro],Tau_14C$XLat[Uni_Pro])
colnames(Long_Lat) <- c("Site","Profile","Longitude","Latitude")
write.table(Long_Lat,file="./Longs_Lats.csv",sep=",",col.names = TRUE,row.names = F)

# the WISE database has 7 consistent layers from 0 to 200cm
WISE30sec0 = read.table("./HW30s_wDi.txt",sep=",",header=T) # the actual data
names(WISE30sec0)
#use only a few variables
WISE30sec = WISE30sec0[,c("NEWSUID","TopDep","BotDep","CFRAG","SDTO","STPC","CLPC","BULK","TAWC","ORGC","TOTN","CNrt","PHAQ","CECS","ECEC","CECc","TEB","BSAT","ESP","ALSA","TCEQ","GYPS","ELCO")]
#select the useful rows
WISE30sec = WISE30sec[WISE30sec$BULK>0&WISE30sec$ORGC>0&WISE30sec$CFRAG>0,]
### input the gis information (long and lat) from WISE database:
AllData = raster("./w001000.adf")
profilecells = extract(AllData,XY) #extracting
AllAttributes = levels(AllData)[[1]]
profileNEWSUID = as.character(AllAttributes$NEWSUID[match(profilecells,AllAttributes$ID)])

tt <- 0
for (i in c(1:length(XLong))) {
  tt[i] <- match(profileNEWSUID[i],WISE30sec$NEWSUID)
}  
nas = which(is.na(tt)==TRUE)
nas
Tau_14C <- Tau_14C[-nas,]
XY <- data.frame(Tau_14C$XLong,Tau_14C$XLat)
profilecells = extract(AllData,XY) #extracting
AllAttributes = levels(AllData)[[1]]
profileNEWSUID = as.character(AllAttributes$NEWSUID[match(profilecells,AllAttributes$ID)])

Tau_14C <- Tau_14C[Tau_14C$Xlayer_bot>0&Tau_14C$Xlayer_top>0&is.na(Tau_14C$Xlayer_top)==F&Tau_14C$Xlayer_bot<1000&Tau_14C$Xlayer_top<1000,]
length(Tau_14C[,1])
summary(Tau_14C)
### Data harmonization 
source("ea_spline.R")
vars <- c("CFRAG","SDTO","STPC","CLPC","BULK","TAWC","ORGC","TOTN","CNrt","PHAQ","CECS","ECEC","CECc","TEB","BSAT","ESP","ALSA","TCEQ","GYPS","ELCO")
WISE_estimated = data.frame(array(NA,dim=c(length(Tau_14C$XLong),length(vars))))
colnames(WISE_estimated) <- c("CFRAG","SDTO","STPC","CLPC","BULK","TAWC","ORGC","TOTN","CNrt","PHAQ","CECS","ECEC","CECc","TEB","BSAT","ESP","ALSA","TCEQ","GYPS","ELCO")

for (i in c(1:length(Tau_14C$XLong))) {
  print(i)
  data <- WISE30sec[WISE30sec$NEWSUID==profileNEWSUID[i],]
  if (length(data[,1])==7) {
    print(as.character(data$NEWSUID)[1])
    SoilID <- as.character(data$NEWSUID)[1]
    Upper <- data$TopDep
    Lower <- data$BotDep
    # if (Tau_14C$Xlayer_bot[i] < 200) {
    for (j in vars) {
      value <- data[,j]
      #names(value) <- "value"
      obj = as.data.frame(cbind(SoilID,Upper,Lower,value))
      obj = obj[complete.cases(obj),]
      rownames(obj) <- NULL
      Xtop <- Tau_14C$Xlayer_top[i];Xbot <- Tau_14C$Xlayer_bot[i]
      if ((Xbot-Xtop) < 1) {Xbot <- Xtop+1}  #the spline function requires so!!!!!
      sp.fit<- ea_spline(obj, var.name="value",d=t(c(Xtop,Xbot)))
      WISE_estimated[i,j] <- sp.fit$harmonised[2]
    }
    #} else {
    #      WISE_estimated[i,] <- data[7,c(2:21)]
    #  }
  } else {
    WISE_estimated[i,] <- rep(NA,length(vars))
  }
}
Tau_14C[1001,]
Tau_Wise <- cbind(Tau_14C,WISE_estimated)
rownames(Tau_Wise) <- NULL
summary(Tau_Wise)


## Make the Tau_14C$XOC complete, by using ORGC from WISE if it is NA
Tau_Wise[is.na(Tau_Wise$XOC)==T,"XOC"] <- Tau_Wise[is.na(Tau_Wise$XOC)==T,"ORGC"]/10 # He's data is %; while WISE and WOSIS are g/kg
Tau_Wise[is.na(Tau_Wise$XBD)==T,"XBD"] <- Tau_Wise[is.na(Tau_Wise$XBD)==T,"BULK"] 
Tau_Wise[is.na(Tau_Wise$XCorase)==T,"XCorase"] <- Tau_Wise[is.na(Tau_Wise$XCorase)==T,"CFRAG"] 
summary(Tau_Wise)
save(Tau_Wise ,file="Tau_Wise.Rdata")  
Tau_Wise_R = Tau_Wise[is.na(Tau_Wise$CFRAG)==F,]
summary(Tau_Wise_R)
save(Tau_Wise_R,file="Tau_Wise_R.Rdata")  


### BIOMEs extracting #### 
Tau_Wise <- get(load("Tau_Wise_R.Rdata")) # start with this new data.frame
XY <- data.frame(Tau_Wise$XLong,Tau_Wise$XLat)
Biomes <- raster("E:/Zhongkui/Global_SOC_Warming/Data/Biomes/biomes/w001000.adf")
Extracted_Biomes = extract(Biomes,XY)### created around 136 NAs, complete them by using the nearest values #
nas = which(is.na(Extracted_Biomes)==TRUE)
Extracted_Biomes[nas]=extract(Biomes,XY[nas,],buffer=10000,fun = function(x) median(x,na.rm=T)) #Biome is unique values, can not be interpolated or averaged
unique(Extracted_Biomes) #avoid strange values and NA values
length(Extracted_Biomes)
Tau_Wise_Biomes <- cbind(Tau_Wise,Extracted_Biomes)
names(Tau_Wise_Biomes)[length(names(Tau_Wise_Biomes))] <- "Biomes"
summary(Tau_Wise_Biomes)

### Bioclimatic variables extracting #### 
ClimVars = as.data.frame(array(NA,dim=c(length(Tau_Wise_Biomes[,1]),19)))
colnames(ClimVars) = paste("BIO",1:19,sep="")
index = 0
for (ii in c(1:19)) {
  index = index + 1
  BIOx = raster(paste("./wc2.1_30s_bio/wc2.1_30s_bio_",ii,".tif",sep=""))
  BIOy= extract(BIOx,XY)
  nas = which(is.na(BIOy)==TRUE)
  BIOy[nas]=extract(BIOx,XY[nas,],buffer=10000,fun = function(x) mean(x,na.rm=T)) #buffer=5000 creats NA, so increase it to 10000
  ClimVars[,index] = BIOy	
}
#save(ClimVars,file="./ClimVars.Rdata")
## Put Tau, Wise, Biome and Climate into a single data frame
Tau_Wise_Biomes_Climate = as.data.frame(cbind(Tau_Wise_Biomes,ClimVars))
summary(Tau_Wise_Biomes_Climate)
save(Tau_Wise_Biomes_Climate,file="./Tau_Wise_Biomes_Climate.Rdata")

#Calculate C input and combine it to the big dataset
Tau_Wise_Biomes_Climate <- get(load("Tau_Wise_Biomes_Climate.Rdata"))
Tau_Wise_Biomes_Climate$Input_mean <- Tau_Wise_Biomes_Climate$XOC*Tau_Wise_Biomes_Climate$XBD*(Tau_Wise_Biomes_Climate$Xlayer_bot-Tau_Wise_Biomes_Climate$Xlayer_top)*(1-Tau_Wise_Biomes_Climate$XCorase/100)/Tau_Wise_Biomes_Climate$tau_mean

#Tau_Wise_Biomes_Climate$Input_up <- Tau_Wise_Biomes_Climate$XOC*Tau_Wise_Biomes_Climate$XBD*(Tau_Wise_Biomes_Climate$Xlayer_bot-Tau_Wise_Biomes_Climate$Xlayer_top)*(1-Tau_Wise_Biomes_Climate$XCorase/100)/Tau_Wise_Biomes_Climate$tau_low

#Tau_Wise_Biomes_Climate$Input_low <- Tau_Wise_Biomes_Climate$XOC*Tau_Wise_Biomes_Climate$XBD*(Tau_Wise_Biomes_Climate$Xlayer_bot-Tau_Wise_Biomes_Climate$Xlayer_top)*(1-Tau_Wise_Biomes_Climate$XCorase/100)/Tau_Wise_Biomes_Climate$tau_up
#The following file can be transferred to Xiao LJ
save(Tau_Wise_Biomes_Climate,file="./Tau_Wise_Biomes_Climate_Input.Rdata")

summary(Tau_Wise_Biomes_Climate)
length(Tau_Wise_Biomes_Climate[,1])


LandCover = raster("/MOD12Q1_UMD.tif")
#Biomes = Biomes_Paper[extract(Biomeslayer,xy)]
#Biomes.na = which(is.na(Biomes)==TRUE)
#Biomes[Biomes.na] = Biomes_Paper[extract(Biomeslayer,xy[Biomes.na,],buffer=5000,fun = function(x) (getmode(v=x)))]
#Biomes[which(extract(LandCover,xy)==12)] = 9
## Land cover
iLandCover = extract(LandCover,XY)
df <- Tau_Wise_Biomes_Climate
XY <- data.frame(df$XLong,df$XLat)
length(XY[,1])
iLandCover = extract(LandCover,XY)
summary(iLandCover)
LJ <- cbind(df,iLandCover)
LJ$Biomes[LJ$Biomes<=3]=1
LJ$Biomes[LJ$Biomes>3 & LJ$Biomes<5]=3
LJ$Biomes[LJ$Biomes==7]=2
LJ$Biomes[LJ$Biomes==8]=4
LJ$Biomes[LJ$Biomes==9]=2
LJ$Biomes[LJ$Biomes==10]=5
LJ$Biomes[LJ$Biomes==11]=7
LJ$Biomes[LJ$Biomes==12]=5
LJ$Biomes[LJ$Biomes==13]=8
LJ$Biomes[LJ$Biomes==14]=1
LJ$Biomes[LJ$iLandCover==12]=9


############################################################################################
###3 Train Cage model
############################################################################################

#3.1 read data
bnpp=read.csv("./14C Data.csv",header = T,as.is = T)%>%
  select(Code,XLong,XLat,tau_mean,Biomes,Order,Xlayer_top,Xlayer_bot,
         BIO1,BIO2,BIO8,BIO12,BIO14,BIO15,BIO18,
         CECc,CLPC,CNrt,ELCO,ESP,GYPS,PHAQ,SDTO,TCEQ,MNPP)%>%
  mutate(Biomes=as.factor(Biomes),Order=factor(Order,levels=c("Alfisols","Andisols" ,"Aridisols","Entisols","Gelisols",
                                                              "Histosols","Ice/Glacier","Inceptisols","Mollisols",
                                                              "Oxisols","Spodosols","Ultisols","Vertisols","Water",
                                                              "Rocky Land","Shifting Sands","Salt"),
                                               labels = c(1:14,4,4,4)))%>%
  filter(tau_mean<50000)%>%mutate(tau_mean=log(tau_mean))


#3.2 datasplit
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
Datasplit=SplitFun(all.data,train_ratio=0.8,seed=31) 

# select predictors
input_ix=c(4,5:(Datasplit$vncol))
colnames(Datasplit$train_data)[input_ix]

##3.3 Machine Learning model building
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


## 3.4 Train model
TrainoutRF=ModelML(Datasplit,input_ix,mod="RF")


##1.5 compile results
##Prediction for trainning data
ddataF=data.frame(tau_mean=TrainoutRF$train_data$tau_mean)
dataF$RF=TrainoutRF$train_data$Pred

##prediction for testing data
dataT=data.frame(tau_mean=TrainoutRF$test_data$tau_mean)
dataT$RF=TrainoutRF$test_data$Pred

dataF$Type="Calibration"
dataT$Type="Validation"

AllPre=rbind(dataF,dataT)
AllPre=AllPre%>%gather(key = Model,value = Pred,-tau_mean,-Type)

RMSEFun=function(simvalue,obsvalue) {round(sqrt(mean((simvalue-obsvalue)^2)),2)}
R2Fun=function(simvalue,obsvalue){round(summary(lm(simvalue~1+obsvalue))$r.squared,2)}


Evaluation=AllPre%>%dplyr::group_by(Type,Model)%>%
  dplyr::summarise(RMSE=RMSEFun(Pred,tau_mean),R2=R2Fun(Pred,tau_mean))%>%filter(Type=="Validation")

##3.6 save results
saveRDS(Evaluation,"Cage_Evaluation.rds")
saveRDS(AllPre,"Cage_Prediction.rds")




############################################################################################
###4 Model evaluation and plot
############################################################################################

ANPPEvaluation=readRDS("ANPP_Evaluation.rds")
ANPPpred=readRDS("ANPP_Prediction.rds")%>%select(ANPP:Pred)%>%rename(Obs=ANPP)
ANPPpred$Index="ANPP"
ANPPEvaluation$Index="ANPP"

CageEvaluation=readRDS("Cage_Evaluation.rds")
Cagepred=readRDS("Cage_Prediction.rds")%>%select(tau_mean:Pred)%>%rename(Obs=tau_mean)
Cagepred$Index="Cage"
CageEvaluation$Index="Cage"


BNPPEvaluation=readRDS("./NPP/BNPP_Evaluation.rds")
BNPPpred=readRDS("./NPP/BNPP_Prediction.rds")%>%select(bnpp:Pred)%>%rename(Obs=bnpp)
BNPPpred$Index="BNPP"
BNPPEvaluation$Index="BNPP"

AllPre=rbind(Cagepred,BNPPpred,ANPPpred)
Evaluation=rbind(CageEvaluation,BNPPEvaluation,ANPPEvaluation)

AllPre=AllPre%>%filter(Model=="RF")%>%
  mutate(Index=factor(Index,levels = c("Cage","BNPP","ANPP")))
Evaluation=Evaluation%>%filter(Model=="RF")%>%
  mutate(Index=factor(Index,levels = c("Cage","BNPP","ANPP")))

Cols = pal_npg("nrc", alpha = 0.8)(9)[c(1,2)]
library("scales")
show_col(Cols)

windows(12,5)
AllPre%>%
  ggplot(.,aes(x=Obs,y=Pred,group=Type,col=Type))+
  geom_point(aes(fill=Type),col="black",shape=21,size=4,stroke=1.1,alpha=0.45)+
  scale_fill_manual(values=Cols)+
  scale_color_manual(values=Cols)+
  scale_alpha(range = c(.4, .5),guide = FALSE)+
  geom_abline(intercept = 0, slope = 1,size=0.5,linetype=2)+
  scale_x_continuous(limits = c(-4,12),breaks = seq(-4,10,2))+
  scale_y_continuous(limits = c(-4,12),breaks = seq(-4,10,2))+
  
  facet_wrap(~Index, scales="free")+
  labs(x=expression(paste("Log"[e], " observed Cage (g/m2)",sep="")),
       y=expression(paste("Log"[e], " simulated Cage (g/m2)",sep="")))+
  theme(axis.ticks.length = unit(0.2, "cm"))+
  theme(legend.title=element_blank(),
        legend.justification=c(1,0),legend.position=c(1,0),
        strip.text = element_blank(),
        legend.text=element_text(size=18,colour="black"))+
  geom_text(
    data    = Evaluation,
    mapping = aes(x = -Inf, y = Inf, label = Index),size=5,
    hjust   = -0.02, vjust   = 1.3,parse = F, inherit.aes=FALSE
  )+
  geom_text(
    data    = Evaluation,
    mapping = aes(x = -Inf, y = Inf, label = Label),size=5,
    hjust   = -0.02, vjust   = 2,parse = T, inherit.aes=FALSE
  )






















