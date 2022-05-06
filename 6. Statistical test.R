library(raster)
library(dplyr)
library(tidyr)
require(maptools)
library(maps)
library(dplyr)
library(RColorBrewer)
require(rasterVis)
library(ggsci)
library(ggrepel)
library(scales)


#01 aefflux
TransitT=readRDS("./Out/TransitT_vert_MNPP_7layers.rds")

set.seed(120)
TransitT=TransitT%>%gather(layer,value,TransitT_vert_All:TransitT_vert7)%>%mutate(Biomes=as.numeric(Biomes))
TransitTall=TransitT%>%mutate(Biomes =10)
TransitT=rbind(TransitT,TransitTall)
TransitT=TransitT%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                   layer=factor(layer,levels = c(paste0("TransitT_vert",1:7),"TransitT_vert_All")))

Variable="TransitT_vert7"
outline=quantile(TransitT3$value,0.99,na.rm=T)
TransitT%>%drop_na()%>%filter(layer==Variable)%>%mutate(value=log10(value))%>%
  select(Biomes,value)%>%mutate(Biomes=factor(Biomes,levels = 1:10))%>%
    dplyr::group_by(Biomes)%>%dplyr::summarise(Transit=mean(value,na.rm=T))
mean(TransitT[TransitT$layer==Variable,"value"],na.rm=T)
  

##HSD test
library(agricolae)
ABCDtest=function(layers="TransitT_vert1"){
  TH_com = HSD.test(aov(value~Biomes, data=TransitT3%>%filter(layer==layers&Biomes!=10)), trt = "Biomes")
  results=TH_com$groups
  results$Biomes=rownames(results)
  results$layer=layers
  return(results)
}

ABCD=lapply( c(paste0("TransitT_vert",1:7),"TransitT_vert_All"),ABCDtest)
ABCD=do.call(rbind,ABCD)
ABCD=ABCD%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                   layer=factor(layer,levels = c(paste0("TransitT_vert",1:7),"TransitT_vert_All")))


cols= rev(c(pal_d3("category10", alpha = 0.7)(9)[c(1,8,9,3,5,6,7,2,4)],"grey50"))
show_col(cols)

#violin and boxplot
PT=TransitT%>%drop_na()%>%filter(value<outline)%>%
  ggplot(aes(x=Biomes,y=value,fill=Biomes))+
  geom_violin(trim=T,col=NA)+
  geom_boxplot(width=0.15,position=position_dodge(0.9),
               outlier.colour = NA,color='black',show.legend = F
  ) +
  scale_color_manual(values =cols )+
  scale_fill_manual(values =cols )+
  # scale_y_continuous(limits = c(0,6))+
  stat_summary(aes(group=Biomes),position=position_dodge(0.78),geom = "point",
               fun = "mean", colour = "blue", size = 1,show.legend = F) +
  stat_summary(aes(group=Biomes,label=round(..y..,2)), position=position_dodge(0.78),fun=mean, geom="text", size=6,
               vjust = -0.5)+
  # scale_y_log10()+
  coord_flip()+
  facet_wrap(~ layer,ncol = 4)+
  scale_y_continuous(trans = "log10")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme(legend.position = "none")+
  geom_text(data=ABCD,aes(y=1e6,x=Biomes,label = groups, group = Biomes),
            position = position_dodge(0.9),
            size = 4, fontface = "bold")




##02 asoil Cage
Cage=Cage%>%select(NO,Biomes,XLat,Cage1:Cage7,Cage200all)%>%gather(layer,value,Cage1:Cage7,Cage200all)%>%mutate(Biomes=as.numeric(Biomes))
Cageall=Cage%>%mutate(Biomes =10)
Cage=rbind(Cage,Cageall)
Cage=Cage%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                             layer=factor(layer,levels = c(paste0("Cage",1:7),"Cage200all")))

##HSD test
library(agricolae)
ABCDtest=function(layers="Cage1"){
  TH_com = HSD.test(aov(value~Biomes, data=Cage3%>%filter(layer==layers&Biomes!=10)), trt = "Biomes")
  results=TH_com$groups
  results$Biomes=rownames(results)
  results$layer=layers
  return(results)
}

ABCD=lapply( c(paste0("Cage",1:7),"Cage200all"),ABCDtest)
ABCD=do.call(rbind,ABCD)
ABCD=ABCD%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                   layer=factor(layer,levels = c(paste0("Cage",1:7),"Cage200all")))


cols= rev(c(pal_d3("category10", alpha = 0.7)(9)[c(1,8,9,3,5,6,7,2,4)],"grey50"))

##violin and boxplot
PCage=Cage%>%ggplot(aes(x=Biomes,y=value,fill=Biomes))+
  geom_violin(trim=T,col=NA)+
  geom_boxplot(width=0.15,position=position_dodge(0.9),
               outlier.colour = NA,color='black',show.legend = F
  ) +
  scale_color_manual(values =cols )+
  scale_fill_manual(values =cols )+
  # scale_y_continuous(limits = c(0,6))+
  stat_summary(aes(group=Biomes),position=position_dodge(0.78),geom = "point",
               fun = "mean", colour = "blue", size = 1,show.legend = F) +
  coord_flip()+
  facet_wrap(~ layer,ncol = 4,scales = "free")+
  theme(legend.position = "none")+
  geom_text(data=ABCD,aes(y=6,x=Biomes,label = groups, group = Biomes),
            position = position_dodge(0.9),
            size = 4, fontface = "bold")






##03 ra=aefflux/asoil
Ratio=readRDS("./Ratio_200all.rds")
Ratio=Cage%>%select(NO,Biomes,XLat)
Ratio$Ratio1=TransitT$TransitT_vert1/Cage$Cage1
Ratio$Ratio2=TransitT$TransitT_vert2/Cage$Cage2
Ratio$Ratio3=TransitT$TransitT_vert3/Cage$Cage3
Ratio$Ratio4=TransitT$TransitT_vert4/Cage$Cage4
Ratio$Ratio5=TransitT$TransitT_vert5/Cage$Cage5
Ratio$Ratio6=TransitT$TransitT_vert6/Cage$Cage6
Ratio$Ratio7=TransitT$TransitT_vert7/Cage$Cage7
Ratio$Ratio200all=TransitT$TransitT_vert_All/Cage$Cage200all
saveRDS(Ratio,"./Ratio_Vert_200all.rds")


Ratio=Ratio%>%select(NO,Biomes,XLat,Ratio1:Ratio7,Ratio200all)%>%gather(layer,value,Ratio1:Ratio7,Ratio200all)%>%mutate(Biomes=as.numeric(Biomes))
Ratioall=Ratio%>%mutate(Biomes =10)
Ratio=rbind(Ratio,Ratioall)
Ratio=Ratio%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                     layer=factor(layer,levels = c(paste0("Ratio",1:7),"Ratio200all")))


#HSD test
library(agricolae)
ABCDtest=function(layers="Ratio1"){
  TH_com = HSD.test(aov(value~Biomes, data=Ratio3%>%filter(layer==layers&Biomes!=10)), trt = "Biomes")
  results=TH_com$groups
  results$Biomes=rownames(results)
  results$layer=layers
  return(results)
}

ABCD=lapply( c(paste0("Ratio",1:7),"Ratio200all"),ABCDtest)
ABCD=do.call(rbind,ABCD)
ABCD=ABCD%>%mutate(Biomes=factor(Biomes,levels = 10:1),
                   layer=factor(layer,levels = c(paste0("Ratio",1:7),"Ratio200all")))


cols= rev(c(pal_d3("category10", alpha = 0.7)(9)[c(1,8,9,3,5,6,7,2,4)],"grey50"))
show_col(cols)
# violin and boxplot
PRatio=Ratio%>%ggplot(aes(x=Biomes,y=value,fill=Biomes))+
  geom_violin(trim=T,col=NA)+
  geom_boxplot(width=0.15,position=position_dodge(0.9),
               outlier.colour = NA,color='black',show.legend = F
  ) +
  scale_color_manual(values =cols )+
  scale_fill_manual(values =cols )+
  # scale_y_continuous(limits = c(0,6))+
  stat_summary(aes(group=Biomes),position=position_dodge(0.78),geom = "point",
               fun = "mean", colour = "blue", size = 1,show.legend = F) +
  coord_flip()+
  facet_wrap(~ layer,ncol = 4)+
  scale_y_continuous(trans = "log10")+
  theme(legend.position = "none")+
  geom_text(data=ABCD,aes(y=4,x=Biomes,label = groups, group = Biomes),
            position = position_dodge(0.9),
            size = 4, fontface = "bold")


