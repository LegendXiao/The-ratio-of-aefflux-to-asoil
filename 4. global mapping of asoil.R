rm(list = ls())
library(raster)
library(dplyr)
library(tidyr)
library(caret)
gc(reset = T)

##################################################################################
#### 1 Cage
##################################################################################

##1 Cage global prediction
time1=Sys.time()
for(i in 1:7){
  fileNO=sort(as.numeric(unlist(strsplit(list.files(paste0("./inputdata/UnLayer",i),pattern = ".rds$"),".rds"))))
  PreF=function(j,i,mod){
    library(ranger)
    input=readRDS(paste0("./inputdata/UnLayer",i,"/",j,".rds"))
    
    NR=nrow(input)
    
    ##missing value
    input_NA_NO=which(is.na(input$BIO1))
    input_NONA_NO=which(!is.na(input$BIO1))
    output=data.frame(NO=input[,1],Cage_mean=rep(0,NR),Cage_sd=rep(0,NR),Layer=i,MissNO=rep(0,NR))
    
    Cagepred=predict(mod,input[input_NONA_NO,],predict.all=T)$predictions
    output$Cage_mean[input_NONA_NO]=round(apply(Cagepred,1,mean),4)
    output$Cage_sd[input_NONA_NO]=round(apply(Cagepred,1,sd),4)
    output$MissNO[input_NA_NO]=1
    
    rm(Cagepred);
    rm(input);
    gc(reset = T)
    return(output)
  }
  
  ###parallel running
  library(snowfall)
  sfInit(parallel = T,cpus = 70)
  sfLibrary(dplyr)
  tictoc::tic()
  results2=sfLapply(fileNO,PreF,i=i,mod=mod_CAge)
  tictoc::toc()
  results2=do.call("rbind",results2)
  saveRDS(results2,paste0("./Out/Cage_Layer",i,".rds"))
  rm(results2)
  sfStop()
  gc(reset = T)
}
print(Sys.time()-time1)


##################################################################################
#### 2 mapping
##################################################################################
##3.2 asoil = Cage
############################################
Cages=paste0("./Products/C14_Cage",1:7,".tif")
Cagesras= stack(Cages)
names(Cagesras) <- c(letters[1:7])
cols=c("#FDE725","#A8DB34","#38BA76","#1F978A","#31688E","#404687")
CageF=function(Rshape=Cagesras[[1]]){
  PCage1=tm_shape(shp = Rshape,ylim = c(-70,90)) +
    tm_raster(
      title = "asoil",n=9,
      palette =  colorRampPalette(cols)(100),
      breaks = c(0,300,600,1500,2000,3000,5000,8000,10000,100000),
      # breaks = c(0,300,600,1000,1500,2000,3000,4000,5000),
      stretch.palette = F,
      legend.show = FALSE,
      style = "cont",
      legend.hist = F
    ) +
    tm_shape(shp = worldCombine) + 
    tm_borders(col = "grey", lwd = 0.7) +
    tm_layout(
      scale = 1,
      legend.position = c("left","bottom"),
      bg.color = NA,  ##背景色
      legend.title.size = 1,
      legend.title.fontface = 2,
      legend.hist.size = 1,
      legend.text.size = 0.8,
      legend.hist.width = .3,
      legend.frame = FALSE,
      legend.outside = F,
      legend.bg.color = NA,
      legend.bg.alpha = 1
    )
  return(PCage1)
}
PCage1=CageF(Cagesras[[1]])
PCage2=CageF(Cagesras[[2]])
PCage3=CageF(Cagesras[[3]])
PCage4=CageF(Cagesras[[4]])
PCage5=CageF(Cagesras[[5]])
PCage6=CageF(Cagesras[[6]])
PCage7=tm_shape(shp = Cagesras[[7]],ylim = c(-70,90)) +
  tm_raster(
    title = "asoil",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,300,600,1500,2000,3000,5000,8000,10000,100000),
    # breaks = c(0,300,600,1000,1500,2000,3000,4000,5000),
    stretch.palette = F,
    legend.is.portrait = F,
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
    legend.hist.width = .3,
    legend.frame = FALSE,
    legend.outside = T,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )

library(grid)
windows(24,22)
print(PCage1, vp = grid::viewport(0.26, 0.88, width = 0.5, height = 0.22))
print(PCage2, vp = grid::viewport(0.7, 0.88, width = 0.5, height = 0.22))
print(PCage3, vp = grid::viewport(0.26, 0.66, width = 0.45, height = 0.22))
print(PCage4, vp = grid::viewport(0.7, 0.66, width = 0.45, height = 0.22))
print(PCage5, vp = grid::viewport(0.26, 0.44, width = 0.45, height = 0.22))
print(PCage6, vp = grid::viewport(0.7, 0.44, width = 0.45, height = 0.22))
print(PCage7, vp = grid::viewport(0.358, 0.22, width = 0.65, height = 0.22))
