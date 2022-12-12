#################################################
## Settings
#################################################

library(ggplot2)
library(cowplot)
library(viridis)
library(wesanderson)
library (tidyverse)
library(gridExtra)
library(wesanderson)

#################################################
## Figure 2 : presentation of the island & t°
#################################################



#################################################
## Figure 2 : presentation of the island & t°
#################################################

do_plot_fig2 = function(nr, nc, height_map, climate_maps, N_cycles, name = "fig2.png"){
  
  pal <- wes_palette("Zissou1", 10, type = "continuous")
  data <- expand.grid(X=nr, Y=nc)
  data$Z1 <- c(height_map)
  data$Z2 <- c(climate_maps[[1]])
  data$Z3 <- c(climate_maps[[N_cycles]])
  
  plt1 = ggplot() +
    geom_raster(data = data , aes(x = X, y=Y,fill = Z1)) + 
    scale_fill_stepsn(n.breaks = 12, colours = terrain.colors(12),na.value="white",guide = "none")+
    coord_fixed()+ theme(axis.text        = element_blank(),
                         axis.ticks       = element_blank(),
                         axis.title       = element_blank(),
                         panel.background = element_blank())
  
  dat.temp = data[,c("X","Y","Z2","Z3")]
  names (dat.temp) <- c('X','Y','Initital timestep','Final timestep')
  dat.temp$`Initital timestep`
  dat.temp$`Final timestep`
  library (tidyverse)
  kk = pivot_longer(dat.temp,cols=c(3,4),names_to = 'temperature',values_to = "Temp") 
  t1temp = kk %>% filter (temperature == 'Initial timestep') 
  kk$temperature = factor (kk$temperature,levels = c('Initial timestep','Final timestep'))
  
  plt2 = ggplot() +
    geom_raster(data = kk , aes(x = X, y=Y,fill = Temp)) + 
    facet_wrap(temperature~.)+
    scale_fill_gradientn(limits= c(0,16),colours = pal,na.value="white")+
    coord_fixed()+ theme(axis.text        = element_blank(),
                         axis.ticks       = element_blank(),
                         axis.title       = element_blank(),
                         panel.background = element_blank(),
                         strip.background = element_blank(),
                         strip.text  = element_blank())
  
  plt_temp = grid.arrange(plt1, plt2, ncol = 2, widths=c(1, 2.3))
  plot(plt_temp)
  ggsave(plot=plt_temp,name, width = 15, height = 5)
  return(plt_temp)
}

#################################################
## Figure 3 : presentation costs map
#################################################

do_plot_fig3 = function(nr,nc,cost1,cost2,cost3,name = "fig2.png"){
  data <- expand.grid(X=nr, Y=nc)
  data$C1 <- c(as.matrix(cost1))
  data$C2 <- c(as.matrix(cost2))
  data$C3 <- c(as.matrix(cost3))
  dat.temp = data[,c(2,1,3,4,5)]
  names (dat.temp) <- c('X','Y',"(a) uniform cost","(b) west-to-east cost","(c) altitude cost")
  kk2 = pivot_longer(dat.temp,cols=c(3,4,5),names_to = 'cost_type',values_to = "Cost") 
  kk2$temperature = factor (kk2$cost_type,levels = c("(a) uniform cost","(b) west-to-east cost","(c) altitude cost"))
  plt_cost = ggplot() + coord_fixed()+
    geom_raster(data = kk2 , aes(x = X, y=Y,fill = Cost)) + theme(panel.background = element_blank(),strip.background = element_blank() ,line = element_blank(),axis.title = element_blank(),axis.text = element_blank(), axis.ticks = element_blank())+
    facet_wrap(temperature~.)+
    scale_fill_gradientn(colours = viridis(10),limits=c(100,500))
  ggsave(name, width = 15, height = 5)
  return(plt_cost)
}

#################################################
## Figure 5 : results plot
#################################################

do_plot_fig5 = function(nr,nc,suitability_maps,choice,choix,){
  
  
  
  data <- expand.grid(X=nr, Y=nc)
  data$C1 <- c(as.matrix(cost1))
  data$C2 <- c(as.matrix(cost2))
  data$C3 <- c(as.matrix(cost3))
  dat.temp = data[,c(2,1,3,4,5)]
  names (dat.temp) <- c('X','Y',"(a) uniform cost","(b) west-to-east cost","(c) altitude cost")
  kk2 = pivot_longer(dat.temp,cols=c(3,4,5),names_to = 'cost_type',values_to = "Cost") 
  kk2$temperature = factor (kk2$cost_type,levels = c("(a) uniform cost","(b) west-to-east cost","(c) altitude cost"))
  plt_cost = ggplot() + coord_fixed()+
    geom_raster(data = kk2 , aes(x = X, y=Y,fill = Cost)) + theme(panel.background = element_blank(),strip.background = element_blank() ,line = element_blank(),axis.title = element_blank(),axis.text = element_blank(), axis.ticks = element_blank())+
    facet_wrap(temperature~.)+
    scale_fill_gradientn(colours = viridis(10),limits=c(100,500))
  ggsave(name, width = 15, height = 5)
  return(plt_cost)
}

cm0 = cm[[1]]
sp = apply(cm0,FUN = sum,2)



choix_ma = choix
choix_ma = choix_ma[(choix_ma[,1]!=-1),]

XY_ma = (choix[,c(1,2)]-1/2)/nrow
#points(XY_ma[,2], XY_ma[,1],col="white",pch=19)

possibilities = gsc
possibilities = possibilities +1


indices_pres = matrix(0,length(XY_pres[,1]))
#XY_pres = XY_pres[nrow:1,]
for (i in (1:length(XY_pres[,1]))){
  h = NA
  for (j in (1:length(possibilities[,1]))){
    if ((possibilities[j,1]==(nrow+1-XY_pres[i,1]))&(possibilities[j,2]==(XY_pres[i,2]))){
      h = possibilities[j,3]
      break
    }
  }
  if (!is.na(h)){
    indices_pres[i] = possibilities[h,3]
  }
}
print(indices_pres)

coloni = as.matrix(cm[[1]])[,indices_pres]
colo2 = 0.99999999 - coloni
colo3 = round(1-apply(colo2, FUN= prod,2),2)

#####
##### Afficher la situation SAALEE (si aucune action légitime n'est entreprise)
#####

rvb_tensor = array(0.5,c(nrow,ncol,3))
for (i in 1:length(c(colo3))){
  if (colo3[i]>=0.95){
    #0.5 + colo3[i]/2
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],1]=0.9#R
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],2]=0.0 #G
    rvb_tensor[possibilities[indices_pres[i],1],possibilities[indices_pres[i],2],3]=0.9#B    
  }
  
}
for (j in 1:length(cm[[1]][,1])){
  rvb_tensor[possibilities[j,1],possibilities[j,2],1]= 0.6
  rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0.6
  rvb_tensor[possibilities[j,1],possibilities[j,2],3]= 0.6
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_pres)){
  for (j in 1:length(cm[[1]][,1])){
    #○cm[[1]][j,indices_pres[i]]
    s = cm[[1]][j,indices_pres[i]]
    if (s!=0){
      rvb_tensor[possibilities[j,1],possibilities[j,2],1]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
      #rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0
      rvb_tensor[possibilities[j,1],possibilities[j,2],2]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
    }
  }
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_choix_ma2)){
  for (j in 1:length(cm[[1]][,1])){
    #○cm[[1]][j,indices_pres[i]]
    s = cm[[1]][j,indices_choix_ma2[i]]
    if (s!=0){
      rvb_tensor[possibilities[j,1],possibilities[j,2],1]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
      #rvb_tensor[possibilities[j,1],possibilities[j,2],2]= 0
      rvb_tensor[possibilities[j,1],possibilities[j,2],3]= s + rvb_tensor[possibilities[j,1],possibilities[j,2],1] -  rvb_tensor[possibilities[j,1],possibilities[j,2],1] * s
    }
  }
}
indices_choix_ma2 = indices_choix_ma[indices_choix_ma!=0] 
for (i in 1:length(indices_choix_ma2)){
  x = possibilities[indices_choix_ma2[i],1]
  y = possibilities[indices_choix_ma2[i],2]
  rvb_tensor[x,y,1]=0.3 #R
  rvb_tensor[x,y,2]=0.3 #G
  rvb_tensor[x,y,3]=0.9 #B
}
for (i in 1:length(XY_pres[,1])){
  if (colo3[i]>0.8){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0.3 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0.9 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0.3 #B
  }
  if (colo3[i]<=0.8){
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=0.9 #R
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],2]=0 #G
    rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],3]=0 #B
  }
}

rvb_tensor[is.na(map)]=1
rvb_tensor2 = round(rvb_tensor*255)
raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))

plotRGB(flip(raster_RGB))
#legend("topright",c("Surviving sites", "Naturally colonised sites", "Not surviving sites", "Planting sites AM", "Colonised sites due to AM"), cex=1.0, bty="y",
#       fill=c("green","yellow","red","blue","purple"),inset=.04)




















# ggplot() +
#   geom_raster(data = DSM_HARV_df , 
#               aes(x = x, y = y, 
#                   fill = HARV_dsmCrop)) + 
#   geom_raster(data = DSM_hill_HARV_df, 
#               aes(x = x, y = y, 
#                   alpha = HARV_DSMhill)) +  
#   scale_fill_viridis_c() +  
#   scale_alpha(range = c(0.15, 0.65), guide = "none") +  
#   ggtitle("Elevation with hillshade") +
#   coord_quickmap()
# 
# image.nan.better(map,col=col,zlim=range(map,na.rm=T),outside.below.color='brown',outside.above.color='brown',na.color='navy')
# levelplot(map,col.regions = c(terrain.colors(29)), main="Terrain simulated")
# 
# 
# cc0 = c_maps[[1]]
# data <- expand.grid(X=nr, Y=nc)
# colnames(data) = c("lontitude","latitude")
# 
# library(cowplot)
# 
# p1 <-levelplot(cc0,col.regions = c(viridis(29)), main="Temperatures simulated (2020)",at=seq(0, 16, length=30))
# #ggplot(data, aes(X, Y, fill= cc0[nrow(cc0):1,])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()
# cc1 = c_maps[[N_cycles]]
# #ggplot(data, aes(X, Y, fill= cc1[,ncol(cc1):1])) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=Tmarge) + theme_ipsum()
# p2 <-levelplot(cc1,col.regions = c(viridis(29)), main="Temperatures simulated (2070)",at=seq(0, 16, length=30))
# plot_grid(p1,p2)
# 
# # Plot suitabilities over time
# first_suit = 0*map
# for (t in 0:N_cycles){
#   back_suit = as.matrix(list_suit[[(N_cycles+1)-t]])
#   first_suit[back_suit>0.9]=((N_cycles+1)-(t))
# }
# 
# last_suit = 0*map
# for (t in 0:N_cycles){
#   forw_suit = as.matrix(list_suit[[t+1]])
#   last_suit[forw_suit>0.9]=t
# }
# 
# R = ((last_suit-first_suit)/N_cycles)[nrow(first_suit):1,]
# G = (first_suit/N_cycles)[nrow(first_suit):1,]
# B = (1 - last_suit/N_cycles)[nrow(first_suit):1,]
# rvb_tensor = array(0.5,c(nrow,ncol,3))
# rvb_tensor[,,1]=1*(first_suit==1 & last_suit==N_cycles) # red
# rvb_tensor[,,2]=1*((first_suit!=1 & last_suit==N_cycles)|(first_suit==1 & last_suit!=N_cycles)) # green
# rvb_tensor[,,3]=1*(first_suit==1 & last_suit!=N_cycles) # blue

# rvb_tensor[(rvb_tensor[,,1]+rvb_tensor[,,2]+rvb_tensor[,,3])==0]=0.5
# 
# #rvb_tensor[,,3]= rvb_tensor[,,3] - 1 * (R==0 & G==0 & B==1)
# rvb_tensor[rvb_tensor==-1 | FALSE]=0.5
# rvb_tensor[is.na(map)]=1
# 
# ggplot(data, aes(lontitude, latitude, fill= cost)) + geom_tile() + scale_fill_gradient(low="blue", high="red",limits=c(100,200)) + theme_ipsum()
# 






