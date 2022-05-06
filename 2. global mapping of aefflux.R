rm(list = ls())
library(raster)
library(dplyr)
library(tidyr)
library(caret)
gc(reset = T)


##################################################################################
####1 global BNPP
##################################################################################
fileNO=sort(as.numeric(unlist(strsplit(list.files(paste0("./inputdata/UnLayer",1),pattern = ".rds$"),".rds"))))
PreF=function(j,i,mod_BNPP){
  library(ranger)
  input=readRDS(paste0("./inputdata/UnLayer",i,"/",j,".rds"))
  BNPPpred=predict(mod_BNPP,input,predict.all=T)$predictions
  BNPPpred.mean=apply(BNPPpred,1,mean)
  BNPPpred.sd=apply(BNPPpred,1,sd)
  rm(BNPPpred);
  output=data.frame(NO=input[,1],BNPP.mean=BNPPpred.mean,BNPP.sd=BNPPpred.sd,
                    MODIS.NPP=input$BNPP,MODIS.sd=input$MNPPSD)
  rm(input);rm(BNPPpred.mean);rm(BNPPpred.sd)
  gc(reset = T)
  return(output)
}

###parallel running
library(snowfall)
library(parallel)
sfInit(parallel = T,cpus = 70)
tictoc::tic()
results2=sfLapply(fileNO,PreF,i=1,mod_BNPP=mod_BNPP)
tictoc::toc()
results2=do.call("rbind",results2)
saveRDS(results2,paste0("./Out/BNPP_Layerall.rds"))
rm(results)
sfStop()
gc(reset = T)





##################################################################################
#### 2 global ANPP
##################################################################################
fileNO=sort(as.numeric(unlist(strsplit(list.files(paste0("./inputdata/UnLayer",1),pattern = ".rds$"),".rds"))))
PreF=function(j,i,mod_ANPP){
  library(ranger)
  input=readRDS(paste0("./inputdata/UnLayer",i,"/",j,".rds"))
  ANPPpred=predict(mod_ANPP,input,predict.all=T)$predictions
  ANPPpred.mean=apply(ANPPpred,1,mean)
  ANPPpred.sd=apply(ANPPpred,1,sd)
  rm(ANPPpred);
  output=data.frame(NO=input[,1],ANPP.mean=ANPPpred.mean,ANPP.sd=ANPPpred.sd)
  rm(input);rm(ANPPpred.mean);rm(ANPPpred.sd)
  gc(reset = T)
  return(output)
}

###parallel running
library(snowfall)
library(parallel)
sfInit(parallel = T,cpus = 70)
tictoc::tic()
results2=sfLapply(fileNO,PreF,i=1,mod_ANPP=mod_ANPP)
tictoc::toc()
results2=do.call("rbind",results2)
saveRDS(results2,paste0("./Out/ANPP_Layerall.rds"))
rm(results)
sfStop()
gc(reset = T)



##3 Convert gridcell dataframe to map
for(i in 1:7){
  Cage1=readRDS(paste0("./Out/Cage_Layer",i,".rds"))
  Cage_Layer_1=raster(res=1/120)
  Cage_Layer_1[Cage1$NO]=round(exp(Cage1$Cage_mean),2)
  plot(Cage_Layer_1)
  CV_1=raster(res=1/120)
  CV_1[Cage1$NO]=round((Cage1$Cage_sd/Cage1$Cage_mean),2)
  plot(CV_1)
  writeRaster(Cage_Layer_1,paste0("F:/Paper1/BNPPproducts/Cage",i,".tif"))
  writeRaster(CV_1,paste0("F:/Paper1/BNPPproducts/CV_Cage",i,".tif"))
  gc(reset = T)
}

ANPP=readRDS(paste0("./Out/ANPP_Layerall.rds"))
ANPP_Layer=raster(res=1/120)
ANPP_Layer[ANPP$NO]=round(exp(ANPP$ANPP.mean)/100,2)
CV_ANPP_Layer=raster(res=1/120)
CV_ANPP_Layer[ANPP$NO]=round((ANPP$ANPP.sd/ANPP$ANPP.mean),2)
writeRaster(ANPP_Layer,paste0("./BNPPproducts/ANPP.tif"))
writeRaster(CV_ANPP_Layer,paste0("./BNPPproducts/CV_ANPP.tif"))




###################################################
#4 SOC stock
#4.1 SOC stock from WISE % to Mg ha-1
ReadSOC=function(j,layer){
  Soil1=readRDS(paste0("./inputdata/SoilLayer",layer,"/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))
  return(Soil1)
}



##Save Layer1-7 SOC stock
T1=Sys.time()
for(layer in 1:7){
  inputfiles=list.files("./inputdata/SoilLayer",layer)
  print(layer)
  tictoc::tic()
  library(snowfall)
  sfInit(parallel = T,cpus = 70)
  sfLibrary(dplyr)
  SOC1=sfLapply(1:length(inputfiles),ReadSOC,layer)
  sfStop()
  
  SOC1=do.call(rbind,SOC1)
  SOCraster=raster(res=1/120)
  SOCraster[SOC1$NO]=SOC1$SOC
  plot(SOCraster)
  summary(SOC1$SOC)
  writeRaster(SOCraster,paste0("./Products/SOC",layer,".tif"))
  rm(SOCraster)
  gc(reset = T)
  tictoc::toc()
  print(Sys.time()-T1)
}

#SOC stock for ALL layers
library(raster)
SOCfiles=paste0("./Products/SOC",1:7,".tif")
rast_stack <- stack(SOCfiles)
plot(rast_stack)
SOCraster_all<- calc(rast_stack,sum,na.rm = TRUE)
values(SOCraster_all)[values(SOCraster_all)==0]=NA
sum(!is.na(values(SOCraster_all)))
plot(SOCraster_all)
writeRaster(SOCraster_all,paste0("./Products/SOC_All.tif"))

##4.2 Combing SOC of all layers to one dataframe
library(snowfall)
T1=Sys.time()
#Table 
for(layer in 1:7){
  if(layer==1){
    print(layer)
    tictoc::tic()
    sfInit(parallel = T,cpus = 70)
    sfLibrary(dplyr)
    inputfiles=list.files("./inputdata/SoilLayer",layer)
    SOC1=sfLapply(1:length(inputfiles),ReadSOC,layer)
    sfStop()
    SOC1=do.call(rbind,SOC1)%>%select(NO,SOC)
    colnames(SOC1)[2]=paste0("SOC",layer)
    tictoc::toc()
  }else{
    tictoc::tic()
    print(layer)
    tictoc::tic()
    inputfiles=list.files("./inputdata/SoilLayer",layer)
    sfInit(parallel = T,cpus = 70)
    sfLibrary(dplyr)
    SOC2=sfLapply(1:length(inputfiles),ReadSOC,layer)
    sfStop()
    SOC2=do.call(rbind,SOC2)%>%select(NO,SOC)
    colnames(SOC2)[2]=paste0("SOC",layer)
    SOC1=SOC1%>%left_join(.,SOC2)
    tictoc::toc()
  }
  print(Sys.time()-T1)
}
SOC1=SOC1%>%mutate_at(vars(SOC1:SOC7),.funs =~round(.,3) )
saveRDS(SOC1,paste0("./Out/SOC_7Layer.rds"))



##4.3 Caculating vertical carbon transport
ReadSOC=function(j){
  Soil1=readRDS(paste0("./inputdata/SoilLayer1/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil2=readRDS(paste0("./inputdata/SoilLayer2/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil3=readRDS(paste0("./inputdata/SoilLayer3/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil4=readRDS(paste0("./inputdata/SoilLayer4/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil5=readRDS(paste0("./inputdata/SoilLayer5/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil6=readRDS(paste0("./inputdata/SoilLayer6/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  Soil7=readRDS(paste0("./inputdata/SoilLayer7/",j,".rds"))%>%
    select(NO,Xlayer_top,Xlayer_bot,ORGC,BULK,CFRAG)%>%
    mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))%>%select(NO,SOC)
  
  Soil=cbind(Soil1,Soil2[,2],Soil3[,2],Soil4[,2],Soil5[,2],Soil6[,2],Soil7[,2])
  colnames(Soil)[2:8]=paste0("SOC",1:7)
  rm(Soil1,Soil2,Soil3,Soil4,Soil5,Soil6,Soil7)
  
  
  ##Vertical transport
  Vertical_T <- function(i, dx=c(0.2,0.2,0.2,0.2,0.2,0.5,0.5),  Df=2.74*10^-7*365){
    # LA = 0.1 * anpp
    SOC=unlist(c(Soil[i,2:8]))/dx/10
    SOC[is.na(SOC)] <- 0
    Ad=Df
    N=7
    D =  diff(Df * diff(c(SOC[1],  SOC, 0.8*SOC[N])) / c(dx[1], dx)) / dx 
    # A =  Ad * diff(c(SOC, SOC[N])) / dx
    VT = round(D[2:7],4)*10 #-A
    # VT[1] = LA
    return(VT)
  }
  
  # tictoc::tic()
  results=t(sapply(1:nrow(Soil),Vertical_T))
  results=as.data.frame(results)
  colnames(results)=paste0("SOC",2:7)
  results$NO=Soil$NO
  # tictoc::toc()
  saveRDS(results,paste0("./inputdata/SoilVertical/",j,".rds"))
}

inputfiles=list.files("./inputdata/SoilLayer1")
library(snowfall)
sfInit(parallel = T,cpus = 75)
sfLibrary(dplyr)
sfLapply(1:length(inputfiles),ReadSOC)
sfStop()

##04 Combing SOC stock with vertical carbon transport
ReadVer=function(i){
  readRDS(paste0("./inputdata/SoilVertical/",i,".rds"))
}

sfInit(parallel = T,cpus = 70)
SoilVertical=sfLapply(1:length(inputfiles),ReadVer)
sfStop()

SoilVertical=do.call(rbind,SoilVertical)

##0.1*ANPP will transport to SOC of the first layer
ANPP=readRDS("./ANPP_Layerall.rds")
head(ANPP)
SoilVertical$SOC1=ANPP$ANPP.mean *0.1
VerRaster=raster(res=1/120)
VerRaster[SoilVertical$NO]=SoilVertical$SOC1

saveRDS(SoilVertical,"./Out/Vertical_TransportC.rds")



###################################################
#5 transit time aefflux
#5.1 aefflux=SOC/BNPP
# aefflux without vertical transport
TransitT=BNPP%>%select(NO,XLat,Biomes)
TransitT$Transit_All=SOC$SOC/BNPP$BNPP
TransitT$Transit1=SOC$SOC1/BNPP$BNPP1
TransitT$Transit2=SOC$SOC2/BNPP$BNPP2
TransitT$Transit3=SOC$SOC3/BNPP$BNPP3
TransitT$Transit4=SOC$SOC4/BNPP$BNPP4
TransitT$Transit5=SOC$SOC5/BNPP$BNPP5
TransitT$Transit6=SOC$SOC6/BNPP$BNPP6
TransitT$Transit7=SOC$SOC7/BNPP$BNPP7
TransitT=TransitT%>%mutate_at(vars(Transit_All:Transit7),.funs = ~round(.,2))

#5.2 aefflux with vertical transport
# aefflux=SOC/(BNPP+vertical C)

SoilVertical=readRDS("./Vertical_TransportC.rds")
SoilVertical=SoilVertical%>%mutate(SOC = rowSums(.[grep("SOC", names(.))], na.rm = TRUE))
TransitT_vert=BNPP%>%select(NO,XLat,Biomes)
TransitT_vert$TransitT_vert_All=SOC$SOC/(BNPP$BNPP+SoilVertical$SOC)
TransitT_vert$TransitT_vert1=SOC$SOC1/(BNPP$BNPP1+SoilVertical$SOC1)
TransitT_vert$TransitT_vert2=SOC$SOC2/(BNPP$BNPP2+SoilVertical$SOC2)
TransitT_vert$TransitT_vert3=SOC$SOC3/(BNPP$BNPP3+SoilVertical$SOC3)
TransitT_vert$TransitT_vert4=SOC$SOC4/(BNPP$BNPP4+SoilVertical$SOC4)
TransitT_vert$TransitT_vert5=SOC$SOC5/(BNPP$BNPP5+SoilVertical$SOC5)
TransitT_vert$TransitT_vert6=SOC$SOC6/(BNPP$BNPP6+SoilVertical$SOC6)
TransitT_vert$TransitT_vert7=SOC$SOC7/(BNPP$BNPP7+SoilVertical$SOC7)

TransitT_vert=TransitT_vert%>%mutate_at(vars(TransitT_vert_All:TransitT_vert7),.funs = ~round(.,2))
saveRDS(TransitT_vert,"./Out/TransitT_vert_MNPP_7layers.rds")

for(i in 1:8){
  Traster=raster(res=1/120)
  Traster[TransitT_vert$NO]=TransitT_vert[,3+i]
  plot(Traster,zlim=c(0,1000))
  outName=paste0("./Products/",colnames(TransitT_vert)[3+i],".tif")
  writeRaster(Traster,outName,overwrite=TRUE)
}


###################################################
#6 mapping
world=readRDS("./Map/Worldcountry.rds")
world=world[world@data$NAME!="Antarctica",]
world=st_as_sf(world)
worldCombine=read_sf("./Map/World_combine.shp")
#Without vertical transport
Transitfiles=paste0("./Transit",1:7,".tif")
TransitRas= stack(Transitfiles)
names(TransitRas) <- c(letters[1:7])
#With vertical transport
Transitfiles=c(paste0("./TransitT_vert",1:7,".tif"),"./Products/TransitT_vert_All.tif")
TransitRas= stack(Transitfiles)
names(TransitRas) <- c(letters[1:8])


##6.1 aefflux
#################################
##1 layer
cols=c("#FDE725","#CFE11C","#A8DB34","#38BA76","#1F978A","#31688E","#404687")
PT1=tm_shape(shp = TransitRas[[1]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA, 
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


##2 layer
PT2=tm_shape(shp = TransitRas[[2]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA,  
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


##3 layer
PT3=tm_shape(shp = TransitRas[[3]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA,  
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


##4 layer
PT4=tm_shape(shp = TransitRas[[4]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA,  
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )

windows(12,6)
PT4


##5 layer
PT5=tm_shape(shp = TransitRas[[5]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA,  
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )




#6 layer
PT6=tm_shape(shp = TransitRas[[6]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.show = FALSE,
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA, 
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


#layer 7
PT7=tm_shape(shp = TransitRas[[7]],ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=12,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,30,50,100,200,500,1000,5000,10000),
    stretch.palette = F,
    style = "cont",
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  tm_layout(
    scale = 1,
    legend.position = c("left","bottom"),
    bg.color = NA,  
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.hist.size = 1,
    legend.text.size = 0.8,
    legend.hist.width = .25,
    legend.frame = FALSE,
    legend.outside = T,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )

library(grid)
windows(24,22)
print(PT1, vp = grid::viewport(0.26, 0.88, width = 0.5, height = 0.22))
print(PT2, vp = grid::viewport(0.7, 0.88, width = 0.5, height = 0.22))
print(PT3, vp = grid::viewport(0.26, 0.66, width = 0.45, height = 0.22))
print(PT4, vp = grid::viewport(0.7, 0.66, width = 0.45, height = 0.22))
print(PT5, vp = grid::viewport(0.26, 0.44, width = 0.45, height = 0.22))
print(PT6, vp = grid::viewport(0.7, 0.44, width = 0.45, height = 0.22))
print(PT7, vp = grid::viewport(0.358, 0.22, width = 0.65, height = 0.22))















