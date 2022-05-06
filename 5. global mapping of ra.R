rm(list = ls())
library(raster)
library(dplyr)
library(tidyr)



##ra=aefflux/asoil
############################################
colslow=c("#A01F1C","#E72E31","#FF553E","#FFA98B")
colslow=colorRampPalette(colslow)(7)
colshigh=colorRampPalette(c("#E1EFF8","#41C3EA","#00ABDC"))(2)
color<-c(colslow,colshigh)
library(scales)
show_col(color)

RatioFun=function(Cratioras1){
  PR1=tm_shape(shp = Cratioras1,ylim = c(-70,90)) +
    tm_raster(
      title = "Ratio",
      palette = color,n=10,
      breaks = c(0,0.01,0.02,0.04,0.06,0.08,0.1,0.9,1,10),
      stretch.palette = F,
      style = "cont",
      legend.show = F,
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
      legend.bg.color = NA,
      legend.outside = F,
      legend.bg.alpha = 1
    )
  return(PR1)
}

Cratioras1=TransitRas[[1]]/Cagesras[[1]]
Cratioras2=TransitRas[[2]]/Cagesras[[2]]
Cratioras3=TransitRas[[3]]/Cagesras[[3]]
Cratioras4=TransitRas[[4]]/Cagesras[[4]]
Cratioras5=TransitRas[[5]]/Cagesras[[5]]
Cratioras6=TransitRas[[6]]/Cagesras[[6]]

PR1=RatioFun(Cratioras1)
PR2=RatioFun(Cratioras2)
PR3=RatioFun(Cratioras3)
PR4=RatioFun(Cratioras4)
PR5=RatioFun(Cratioras5)
PR6=RatioFun(Cratioras6)
Cratioras7=TransitRas[[7]]/Cagesras[[7]]
PR7=tm_shape(shp = Cratioras7,ylim = c(-70,90)) +
  tm_raster(
    title = "Ratio",
    palette = color,n=10,
    breaks = c(0,0.01,0.02,0.04,0.06,0.08,0.1,0.9,1,10),
    stretch.palette = F,
    style = "cont",
    legend.is.portrait = F,
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
    legend.bg.color = NA,
    legend.outside = T,
    legend.bg.alpha = 1
  )

windows(24,22)
print(PR1, vp = grid::viewport(0.26, 0.88, width = 0.5, height = 0.22))
print(PR2, vp = grid::viewport(0.7, 0.88, width = 0.5, height = 0.22))
print(PR3, vp = grid::viewport(0.26, 0.66, width = 0.45, height = 0.22))
print(PR4, vp = grid::viewport(0.7, 0.66, width = 0.45, height = 0.22))
print(PR5, vp = grid::viewport(0.26, 0.44, width = 0.45, height = 0.22))
print(PR6, vp = grid::viewport(0.7, 0.44, width = 0.45, height = 0.22))
print(PR7, vp = grid::viewport(0.358, 0.22, width = 0.65, height = 0.22))


##3.4 mapping all for 0-200cm
############################################
#3.4.1 asoil
Cageras=raster("./Products/Cage200.tif")
cols=c("#FDE725","#CFE11C","#A8DB34","#38BA76","#1F978A","#31688E","#404687")
color<-colorRampPalette(cols)(9)
P1=tm_shape(shp = Cageras,ylim = c(-70,90)) +
  tm_raster(
    title = "asoil",n=9,
    palette =  colorRampPalette(cols)(100),
    # breaks = c(0,500,1000,1500,2000,3000,5000,8000,10000),
    breaks = c(0,300,600,1000,1500,2000,3000,4000,5000),
    stretch.palette = F,
    style = "fixed",
    legend.hist = TRUE
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) +  
  
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
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )

#3.4.1 aefflux 
Transit=raster("./Products/Transit_All.tif")
P2=tm_shape(shp = Transit,ylim = c(-70,90)) +
  tm_raster(
    title = "aefflux",n=9,
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,10,20,30,50,80,100,200,500,4000),
    stretch.palette = F,
    style = "fixed",
    legend.hist = TRUE
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) +  
  
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

####ra=aefflux/asoil
Cratio=Transit/Cageras
colslow=c("#A01F1C","#E72E31","#FF553E","#FFA98B")
colslow=colorRampPalette(colslow)(7)
show_col(colslow)
colshigh=colorRampPalette(c("#E1EFF8","#41C3EA","#00ABDC"))(2)
color<-c(colslow,colshigh)
P3=tm_shape(shp = Cratio,ylim = c(-70,90)) +
  tm_raster(
    title = "Ratio",n=9,
    palette = color,
    breaks = c(0,0.01,0.02,0.04,0.06,0.08,0.1,0.99,1.01,10),
    stretch.palette = F,
    style = "fixed",
    legend.hist = TRUE
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) +  
  
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

library(grid)
windows(16,22)
grid.newpage()
print(P1, vp = grid::viewport(0.5, 0.8, width = 1, height = 0.3))
print(P2, vp = grid::viewport(0.485, 0.5, width = 1, height = 0.3))
print(P3, vp = grid::viewport(0.485, 0.2, width = 1, height = 0.3))

##3.5 mapping uncertainty for 0-200cm
############################################
cols=rev(c("#DC2121","#F89C53","#FFFFC1","#97D0A5","#3A8EBA"))
color<-colorRampPalette(cols)(9)
show_col(cols)

PUn1=tm_shape(shp = CageUn,ylim = c(-70,90)) +
  tm_raster(
    title = "CV",n=5,
    palette =  colorRampPalette(cols)(1000),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),
    stretch.palette = F,
    style = "cont",
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) + 
  
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
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


PUn2=tm_shape(shp = TransitUn,ylim = c(-70,90)) +
  tm_raster(
    title = "CV",
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),
    stretch.palette = F,
    style = "cont",
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) + 
  
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
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )


##Ura=Uasoil+Uaefflux
RatioUn=CageUn+TransitUn
PUn3=tm_shape(shp = RatioUn,ylim = c(-70,90)) +
  tm_raster(
    title = "CV",
    palette =  colorRampPalette(cols)(100),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),
    stretch.palette = F,
    style = "cont",
    legend.hist = F
  ) +
  tm_shape(shp = worldCombine) + 
  tm_borders(col = "grey", lwd = 0.7) +
  
  tm_graticules(alpha=1, lines = FALSE,lwd = 1, labels.size = 1,ticks = T,
                labels.inside.frame=F,
                x = seq(-150,150,50),
                y = seq(-60,80,20)) +  
  
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
    legend.outside = F,
    legend.bg.color = NA,
    legend.bg.alpha = 1
  )

library(grid)
windows(16,22)
print(PUn2, vp = grid::viewport(0.5, 0.8, width = 1, height = 0.3))
print(PUn1, vp = grid::viewport(0.5, 0.5, width = 1, height = 0.3))
print(PUn3, vp = grid::viewport(0.485, 0.2, width = 1, height = 0.3))

