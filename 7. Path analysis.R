rm(list = ls())
gc(reset = T)
library(raster)
library(dplyr)
library(tidyr)

C14data=read.csv("./C14DataforSEM.csv")
C14data=C14data%>%mutate(SOC=ORGC/10*BULK*(Xlayer_bot-Xlayer_top)*(1-CFRAG/100))

#Load plspm package
library(plspm)

#latent variable 
dat_blocks <- list(
  Topography=c("LandForm","TWI","mrVBF"),
  climate = c('BIO1', 'BIO12'), 
  soil = c( "Xlayer_top","Order"), 
  Cinput = c('MNPP',"Biomes"), 
  asoil ="asoil",
  aefflux = 'aefflux',
  Ratio="Ratio"
  
)
dat_blocks


#relationship for latent variable  
Topography <- c(0, 0, 0, 0, 0,0,0);length(Topography)
climate <- c(1, 0, 0, 0, 0,0,0);length(climate)
soil <- c(1, 1, 0, 0, 0,0,0);length(soil)
Cinput <- c(1, 1, 1, 0, 0,0,0);length(Cinput)
asoil <- c(1, 1, 1, 1, 0,0,0);length(asoil)
aefflux <- c(0, 1, 1, 1, 1,0,0);length(aefflux)
Ratio <- c(1, 1, 1, 1, 1,1,0);length(Ratio)

dat_path <- rbind(Topography,climate, soil,Cinput,asoil,aefflux,Ratio)
colnames(dat_path) <- rownames(dat_path)
dat_path

dat_modes <- rep('A', 7)
dat_modes

dat=C14data[,c("x","y","Elevdata","LandForm",'BIO1', 'BIO12', 'PHAQ',"BULK","ORGC",'SDTO', "SOC",'CLPC', 'STPC','TOTN' ,'CNrt',"TEB","TAWC",
               'CECc',"ELCO","TCEQ",'CFRAG','MNPP',"Biomes", "Xlayer_top","Xlayer_bot","Order",'aefflux',"asoil","Ratio","TWI","mrVBF","SP","fm","Aspect")]
  
##PLS-PM
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes,scaled = F)
dat_pls
summary(dat_pls)
dat_pls$gof
plot(dat_pls,what="all", arr.pos=.4, box.prop=.4, cex.txt=.8)

str(C14data)


#Parameter Estimate
dat_pls$path_coefs
dat_pls$inner_model

innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

dat_pls$inner_summary

dat_pls$effects

dat_pls$outer_model
outerplot(dat_pls, what = 'loadings', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')
outerplot(dat_pls, what = 'weights', arr.width = 0.1, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray')

#goodness-of-fit
dat_pls$gof

#scores
dat_pls$scores



##
C14data=C14data%>%mutate(Layer=(Xlayer_top+Xlayer_bot)/2)

model<-psem(lme(BIO1~LandForm+Elevdata,random=~1|Biomes,data=C14data, method = "ML"),
            lme(BIO12~Elevdata,random=~1|Biomes,data=C14data, method = "ML"),
            lme(SOC~BIO1+BIO12+Layer+LandForm+Elevdata,random=~1|Biomes,data=C14data, method = "ML"),
            lme(MNPP~BIO12+BIO1,random=~1|Biomes,data=C14data, method = "ML"),
            lme(asoil~BIO1+BIO12+Layer,random=~1|Biomes,data=C14data, method = "ML"),
            lme(aefflux~asoil+BIO1+BIO12+Layer+SOC,random=~1|Biomes,data=C14data, method = "ML"),
            lme(Ratio~asoil+aefflux+SOC+Layer+BIO1+BIO12+MNPP,random=~1|Biomes,data=C14data, method = "ML"),
            data = C14data)
new.summary<-summary(model, .progressBar = F)
new.summary$dTable


