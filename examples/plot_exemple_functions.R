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
library(dplyr)

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
  data$Z3 <- c(climate_maps[[N_cycles]])
  data$Z2 <- c(climate_maps[[1]])
  
  plt1 = ggplot() +
    geom_raster(data = data , aes(x = X, y=Y,fill = Z1)) + 
    scale_fill_stepsn(n.breaks = 12, colours = terrain.colors(12),na.value="white",guide = "none")+
    coord_fixed()+ theme(axis.text        = element_blank(),
                         axis.ticks       = element_blank(),
                         axis.title       = element_blank(),
                         panel.background = element_blank())
  
  dat.temp = data[,c("X","Y","Z2","Z3")]
  names (dat.temp) <- c('X','Y','Initial timestep','Final timestep')
  dat.temp$`Initital timestep`
  dat.temp$`Final timestep`
  library (tidyverse)
  kk = pivot_longer(dat.temp,cols=c(3,4),names_to = 'temperature',values_to = "Temp") 
  t1temp = kk %>% filter (temperature == 'Initial timestep') 
  kk$temperature = factor(kk$temperature,levels = c('Initial timestep','Final timestep'))
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

do_plot_fig2bis = function(nr,nc,cost1,cost2,cost3,name = "fig2.png"){
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
## Figure 3 : check la suitability
#################################################

do_plot_fig3 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        global_suitable_coordinates,
                        presence_map){
  # 0-preprocessing
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.95,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5
  rvb_tensor[,,1][array(suit_only_t0)] = 0.65
  rvb_tensor[,,2][array(suit_only_t0)] = 0.65
  rvb_tensor[,,3][array(suit_only_t0)] = 0.65
  rvb_tensor[,,1][array(suit_both)] = 0.75
  rvb_tensor[,,2][array(suit_both)] = 0.75
  rvb_tensor[,,3][array(suit_both)] = 0.75
  rvb_tensor[,,1][array(suit_only_tN)] = 0.85
  rvb_tensor[,,2][array(suit_only_tN)] = 0.85
  rvb_tensor[,,3][array(suit_only_tN)] = 0.85
  
  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  
  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.6
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  
  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 
  
  
  # 11 - plot tensor
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
  
}

#################################################
## Figure 4 : results plot
#################################################

do_plot_fig4 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        global_suitable_coordinates,
                        colonisation_matrices,
                        presence_map){
  # 0-preprocessing
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.95,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5
  rvb_tensor[,,1][array(suit_only_t0)] = 0.65
  rvb_tensor[,,2][array(suit_only_t0)] = 0.65
  rvb_tensor[,,3][array(suit_only_t0)] = 0.65
  rvb_tensor[,,1][array(suit_both)] = 0.75
  rvb_tensor[,,2][array(suit_both)] = 0.75
  rvb_tensor[,,3][array(suit_both)] = 0.75
  rvb_tensor[,,1][array(suit_only_tN)] = 0.85
  rvb_tensor[,,2][array(suit_only_tN)] = 0.85
  rvb_tensor[,,3][array(suit_only_tN)] = 0.85
  
  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  
  vectors_pres_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  colnames(vectors_pres_effects_on_tN) = c("index","i","effect")
  tmp = vectors_pres_effects_on_tN
  tmp$index = as.factor(tmp$index)
  tmp$effect = 1 - tmp$effect
  vectors_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  
  tmp3 = aggregate(effect~(index),tmp,FUN=function(x) prod(x))
  tmp3$index = as.numeric(as.character(tmp3$index))
  tmp3$effect = 1-tmp3$effect
  tmp3 = tmp3[tmp3$effect>0.1,]
  
  XY_effect = inner_join(tmp3,XY_to_index,by = c("index"))
  rvb_tensor[,,1][cbind(XY_effect$x,XY_effect$y)] = 0.5 - 0.3*XY_effect$effect
  rvb_tensor[,,2][cbind(XY_effect$x,XY_effect$y)] = 0.7 + 0.25*XY_effect$effect
  rvb_tensor[,,3][cbind(XY_effect$x,XY_effect$y)] = 0.5 - 0.3*XY_effect$effect
  
  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.9
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.3
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.3 
  
  surviving_sites = inner_join(XY_pres, XY_effect, by=c('x'='x', 'y'='y'))
  rvb_tensor[,,1][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.3
  rvb_tensor[,,2][cbind(surviving_sites[,1],surviving_sites[,2])]= 0.3
  rvb_tensor[,,3][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.9 
  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 
  
  
  # 11 - plot tensor
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
  
}

#################################################
## Figure 5 : results plot
#################################################

do_plot_fig5 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        choices,
                        global_suitable_coordinates,
                        colonisation_matrices,
                        presence_map,
                        cost){
  # 0-preprocessing
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.95,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5
  rvb_tensor[,,1][array(suit_only_t0)] = 0.65
  rvb_tensor[,,2][array(suit_only_t0)] = 0.65
  rvb_tensor[,,3][array(suit_only_t0)] = 0.65
  rvb_tensor[,,1][array(suit_both)] = 0.75
  rvb_tensor[,,2][array(suit_both)] = 0.75
  rvb_tensor[,,3][array(suit_both)] = 0.75
  rvb_tensor[,,1][array(suit_only_tN)] = 0.85
  rvb_tensor[,,2][array(suit_only_tN)] = 0.85
  rvb_tensor[,,3][array(suit_only_tN)] = 0.85
  
  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  
  vectors_pres_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  colnames(vectors_pres_effects_on_tN) = c("index","i","effect")
  tmp = vectors_pres_effects_on_tN
  tmp$index = as.factor(tmp$index)
  tmp$effect = 1 - tmp$effect
  vectors_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  choices = choices[choices[,1]!=0,] 
  
  vectors_intro_effects_on_tN = summary(Matrix(apply(choices, 1, function (x) cm[[x[3]]][,x[5]]),sparse=T))
  colnames(vectors_intro_effects_on_tN) = c("index","i","effect")
  tmp2 = vectors_intro_effects_on_tN
  tmp2$effect = 1 - tmp2$effect
  tmp2$index = as.factor(tmp2$index)
  
  tmp_tot = rbind(tmp,tmp2)
  tmp3 = aggregate(effect~(index),tmp_tot,FUN=function(x) prod(x))
  tmp3$index = as.numeric(as.character(tmp3$index))
  tmp3$effect = 1-tmp3$effect

  tmp3 = tmp3[tmp3$effect>0.1,]
  
  XY_effect = inner_join(tmp3,XY_to_index,by = c("index"))
  rvb_tensor[,,1][cbind(XY_effect$x,XY_effect$y)] = 0.5 - 0.3*XY_effect$effect
  rvb_tensor[,,2][cbind(XY_effect$x,XY_effect$y)] = 0.7 + 0.25*XY_effect$effect
  rvb_tensor[,,3][cbind(XY_effect$x,XY_effect$y)] = 0.5 - 0.3*XY_effect$effect
  
  # 3 - add finally conquered area (use : suitability_maps, choices)
  rvb_tensor[,,1][cbind(choices[,1]+1,choices[,2])] = 0.8-0.5*choices[,3]/N_cycles
  rvb_tensor[,,2][cbind(choices[,1]+1,choices[,2])] = 0.8-0.5*choices[,3]/N_cycles
  rvb_tensor[,,3][cbind(choices[,1]+1,choices[,2])] = 0.1

  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.9
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.3
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.3 
  
  surviving_sites = inner_join(XY_pres, XY_effect, by=c('x'='x', 'y'='y'))
  rvb_tensor[,,1][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.3
  rvb_tensor[,,2][cbind(surviving_sites[,1],surviving_sites[,2])]= 0.3
  rvb_tensor[,,3][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.9 
  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  
  # 11 - plot tensor
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  
  
  # 12 - do cost plot
  nr = seq(0,1,length=nrow)
  nc = seq(0,1,length=ncol)
  data <- expand.grid(X=nr, Y=nc)
  data$Cost <- c(as.matrix(cost1))
  ggplot() + 
    coord_fixed()+
    scale_fill_continuous(na.value = "white")+
    geom_raster(data = data , aes(x = X, y=Y,fill = Cost)) + 
    theme(panel.background = element_blank(),
          strip.background = element_blank() ,
          line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # 14 - all plot
  plotRGB(flip(raster_RGB,1))
  
  colnames(choices) = c("x","y","time_step","cost","index","impact")
  ggplot(as.data.frame(choices), aes(x=time_step)) + 
    xlim(0, 30)+
    ylim(0, 8)+
    geom_histogram(color="black", fill="white",bins = 30)+
    theme_minimal(base_size = 40)
}

#################################################
## Figure 6 : Density choice fonction
#################################################

do_plot_fig6 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        global_suitable_coordinates,
                        presence_map,
                        choicesX){
  # 0-preprocessing
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.95,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5
  rvb_tensor[,,1][array(suit_only_t0)] = 0.65
  rvb_tensor[,,2][array(suit_only_t0)] = 0.65
  rvb_tensor[,,3][array(suit_only_t0)] = 0.65
  rvb_tensor[,,1][array(suit_both)] = 0.75
  rvb_tensor[,,2][array(suit_both)] = 0.75
  rvb_tensor[,,3][array(suit_both)] = 0.75
  rvb_tensor[,,1][array(suit_only_tN)] = 0.85
  rvb_tensor[,,2][array(suit_only_tN)] = 0.85
  rvb_tensor[,,3][array(suit_only_tN)] = 0.85
  
  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  
  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.6
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  
  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 
  
  NN = length(choicesX)
  for (i in 1:NN){
    choice = choicesX[[i]]
    tmp_c= choice[choice[,1]!=0,]
    for (j in 1:length(tmp_c[,1])){
      rvb_tensor[,,1][cbind(tmp_c[j,1],tmp_c[j,2])] = 0.8
      rvb_tensor[,,2][cbind(tmp_c[j,1],tmp_c[j,2])] = rvb_tensor[,,2][cbind(tmp_c[j,1],tmp_c[j,2])] + rvb_tensor[,,2][cbind(tmp_c[j,1]+1,tmp_c[j,2]+1)]*0.8*0.01
      rvb_tensor[,,3][cbind(tmp_c[j,1],tmp_c[j,2])] = 0.2 
    }
  }
  # 11 - plot tensor
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
  
}




do_plot_fig7 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        global_suitable_coordinates,
                        presence_map,
                        countmap){
  # 0-preprocessing
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.95,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5
  rvb_tensor[,,1][array(suit_only_t0)] = 0.65
  rvb_tensor[,,2][array(suit_only_t0)] = 0.65
  rvb_tensor[,,3][array(suit_only_t0)] = 0.65
  rvb_tensor[,,1][array(suit_both)] = 0.75
  rvb_tensor[,,2][array(suit_both)] = 0.75
  rvb_tensor[,,3][array(suit_both)] = 0.75
  rvb_tensor[,,1][array(suit_only_tN)] = 0.85
  rvb_tensor[,,2][array(suit_only_tN)] = 0.85
  rvb_tensor[,,3][array(suit_only_tN)] = 0.85
  
  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  
  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.6
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.1
  
  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 
  
  xx = which(countmap33>0,arr.ind=T)[,1]
  yy = which(countmap33>0,arr.ind=T)[,2]
  for (i in 1:length(xx)){
    rvb_tensor[xx[i],yy[i],1] = 0.2
    rvb_tensor[xx[i],yy[i],2] = 0.2
    rvb_tensor[xx[i],yy[i],3] = 0.8
  }
  for (i in 1:length(xx)){
    rvb_tensor[xx[i],yy[i],1] = rvb_tensor[xx[i],yy[i],1] + 0.6 * countmap[xx[i],yy[i]]/100
    rvb_tensor[xx[i],yy[i],2] = rvb_tensor[xx[i],yy[i],2] - 0.0 * countmap[xx[i],yy[i]]/100
    rvb_tensor[xx[i],yy[i],3] = rvb_tensor[xx[i],yy[i],3] - 0.6 * countmap[xx[i],yy[i]]/100
  }

  
  # 11 - plot tensor
  rvb_tensor[is.na(rvb_tensor)] = NA
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
  
}







do_plot_fig8 = function(nrow,
                        ncol,
                        N_cycles,
                        height_map,
                        suitability_maps,
                        global_suitable_coordinates,
                        presence_map,
                        countmap33){
  # 0-preprocessing
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  rvb_tensor = array(0.9,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  # 2 - add suitability (use : suitability_maps)
  suit_t0 = suitability_maps[[1]]
  suit_tN = suitability_maps[[length(suitability_maps)]]
  suit_only_t0 = suit_t0>0.5 & suit_tN<0.5
  suit_both = suit_t0>0.5 & suit_tN>0.5
  suit_only_tN = suit_t0<0.5 & suit_tN>0.5

  # 3 - add naturally conquered area (use : suitability_maps, presence_map)
  XY_pres = data.frame(which(presence_map==1,arr.ind = T))
  colnames(XY_pres) = c("x","y")
  XY_to_index = data.frame(gsc+1)
  colnames(XY_to_index) = c("x","y","index")
  pres_index = inner_join(XY_pres,XY_to_index,by = c("x", "y"))
  

  
  # 9 - add initial presence site. By default red, then blue if it survived
  rvb_tensor[,,1][cbind(XY_pres[,1],XY_pres[,2])] = 0.9
  rvb_tensor[,,2][cbind(XY_pres[,1],XY_pres[,2])]= 0.5
  rvb_tensor[,,3][cbind(XY_pres[,1],XY_pres[,2])] = 0.3
  
  vectors_pres_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  colnames(vectors_pres_effects_on_tN) = c("index","i","effect")
  tmp = vectors_pres_effects_on_tN
  tmp$index = as.factor(tmp$index)
  tmp$effect = 1 - tmp$effect
  vectors_effects_on_tN = summary((cm[[1]][,pres_index[,3]]))
  
  tmp3 = aggregate(effect~(index),tmp,FUN=function(x) prod(x))
  tmp3$index = as.numeric(as.character(tmp3$index))
  tmp3$effect = 1-tmp3$effect
  tmp3 = tmp3[tmp3$effect>0.1,]
  
  XY_effect = inner_join(tmp3,XY_to_index,by = c("index"))
  
  surviving_sites = inner_join(XY_pres, XY_effect, by=c('x'='x', 'y'='y'))
  rvb_tensor[,,1][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.3
  rvb_tensor[,,2][cbind(surviving_sites[,1],surviving_sites[,2])]= 0.3
  rvb_tensor[,,3][cbind(surviving_sites[,1],surviving_sites[,2])] = 0.9 
  


  
  # 10 - add black belt
  external_coord = which(is.na(height_map[2:(nrow-1),2:(ncol-1)]),arr.ind=TRUE)
  is_on_border = apply(external_coord,1,FUN = function(x) 0!=sum(!is.na(height_map[(x[1]):(x[1]+2),(x[2]):(x[2]+2)])) )
  coord_border = external_coord[is_on_border,]
  rvb_tensor[,,1][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2
  rvb_tensor[,,2][cbind(coord_border[,1]+1,coord_border[,2]+1)]= 0.2
  rvb_tensor[,,3][cbind(coord_border[,1]+1,coord_border[,2]+1)] = 0.2 
  

  
  xx = which(countmap33>0,arr.ind=T)[,1]
  yy = which(countmap33>0,arr.ind=T)[,2]
  for (i in 1:length(xx)){
    rvb_tensor[xx[i],yy[i],1] = 0.8
    rvb_tensor[xx[i],yy[i],2] = 0.8
    rvb_tensor[xx[i],yy[i],3] = 0.8
  }
  xx = which(countmap33>20,arr.ind=T)[,1]
  yy = which(countmap33>20,arr.ind=T)[,2]
  for (i in 1:length(xx)){
    rvb_tensor[xx[i],yy[i],1] = 0.4
    rvb_tensor[xx[i],yy[i],2] = 0.4
    rvb_tensor[xx[i],yy[i],3] = 0.4
  }

  xx = which(countmap33>50,arr.ind=T)[,1]
  yy = which(countmap33>50,arr.ind=T)[,2]
  for (i in 1:length(xx)){
    rvb_tensor[xx[i],yy[i],1] = 0.2
    rvb_tensor[xx[i],yy[i],2] = 0.2
    rvb_tensor[xx[i],yy[i],3] = 0.2
  }
  
  # 11 - plot tensor
  rvb_tensor[is.na(rvb_tensor)] = NA
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
  
}























do_plot_fig001 = function(nrow,
                          ncol,
                          N_cycles,
                          height_map,
                          suitability_maps,
                          gsc,
                          colonisation_matrices,
                          presence_map,
                          cost){
  # 0-preprocessing
  cost = cost1
  val_tensor = array(0,c(nrow,ncol))
  ind_tensor = array(0,c(nrow,ncol))
  
  indices_pres2 = matrix(0,N_cycles,length(gsc[,1]))
  for (i in 1:N_cycles){
    cm_loc = cm[[i]]
    for (j in 1:length(gsc[,1])){
      indices_pres2[i,j] = sum(cm_loc[,j])
    }
  }
  opt_each_site = apply(indices_pres2,FUN=max,2)
  library(ramify)
  ind_each_site = (argmax(t(indices_pres2)))
  
  diag_cms = matrix(0,N_cycles,length(gsc[,1]))
  for (i in 1:N_cycles){
    cm_loc = cm[[i]]
    diag_cms[i,] = apply(cm_loc,2,sum)
  }
  
  maxi_diags = apply(diag_cms,FUN=max,2)
  
  plot(maxi_diags,opt_each_site)
  
  efficienty3 = opt_each_site / cost[gsc[,c(1,2)]+1]
  eff3 = efficienty3 / max(efficienty3)
  eff_tensor3 = array(0.5,c(nrow,ncol,3))
  
  for (i in 1:length(efficienty3)){
    val_tensor[gsc[i,1]+1,gsc[i,2]+1]= eff3[i]
    ind_tensor[gsc[i,1]+1,gsc[i,2]+1]= ind_each_site[i]
  }
  
  
  rvb_tensor = array(0.9,c(nrow,ncol,3))
  rvb_tensor[is.na(height_map)]=1
  
  rvb_tensor[,,1] = val_tensor
  
  
  # 1 - Structure of the RVB tensor (use : height maps, nrow, ncol)
  
  rvb_tensor[is.na(height_map)]=1
  rvb_tensor2 = round(rvb_tensor*255)
  raster_RGB = stack(raster(rvb_tensor2[,,1]),raster(rvb_tensor2[,,2]),raster(rvb_tensor2[,,3]))
  plt1 <- plotRGB(flip(raster_RGB,1))
}
do_plot_fig001(nrow,ncol,N_cycles,height_map,suitability_maps,gsc, cm,presence_map,cost1)


























opt_each_site = apply(indices_pres2,FUN=max,2)

diag_cms = matrix(0,N_cycles,length(gsc[,1]))
for (i in 1:N_cycles){
  cm_loc = cm[[i]]
  diag_cms[i,] = diag(cm_loc)
}

maxi_diags = apply(diag_cms,FUN=max,2)

plot(maxi_diags,opt_each_site)

### ### ###
efficienty3 = maxi_diags / cost[gsc[,c(1,2)]+1]
eff3 = efficienty3 / max(efficienty3)
#eff = (log(efficienty)-min(log(efficienty)))/(max(log(efficienty))-min(log(efficienty)))

eff_tensor3 = array(0.5,c(nrow,ncol,3))

for (i in 1:length(efficienty3)){
  #rvb_tensor[nrow+1-XY_pres[i,1],XY_pres[i,2],1]=1 #R
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,2]= eff3[i]  #G
  eff_tensor3[gsc[i,1]+1,gsc[i,2]+1,3]= 1-eff3[i]#B
}

eff_tensor3 [is.na(map)]=1
eff_tensor3 = round(eff_tensor3*255)
raster_RGB = stack(raster(eff_tensor3[,,1]),raster(eff_tensor3[,,2]),raster(eff_tensor3[,,3]))
#graphics.off()
plotRGB(flip(raster_RGB))
legend("topright",c("High effect ratio", "Low effect ratio"), cex=1.0, bty="y",
       fill=c("green","purple"),inset=.04)



































































library(ggplot2)
theme_set(theme_minimal())

a = data.frame(0,c(2,2))
a$X = c(0,1)
a$Colonization.probability = c(0,1)
sp <- ggplot(a, aes(X, Colonization.probability))+
  theme(legend.text = element_text(size=10))+
  geom_point(aes(color = Colonization.probability))
sp + scale_color_gradient(low = "darkgreen", high = "green3")


library(ggplot2)
theme_set(theme_minimal())

a = data.frame(0,c(2,2))
a$X = c(1,20)
a$Introduction.time = c(1,20)
sp <- ggplot(a, aes(X, Introduction.time))+
  theme(legend.text = element_text(size=10))+
  geom_point(aes(color = Introduction.time))
sp + scale_color_gradient(low = "yellow", high = "brown",breaks=c(1,10,20))









colfunc <- colorRampPalette(c("green", "darkgreen"))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

layout(matrix(1:2,ncol=2), width = c(1,1),height = c(1,1))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', 
     main = 'Colonization probability', cex.main=2.5)
text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex = 2)
rasterImage(legend_image, 0, 0, 1,1)




colfunc <- colorRampPalette(c("yellow", "brown"))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', 
     main = 'Colonization probability', cex.main=2)
text(x=1.3, y = seq(0,1,l=5), labels = seq(0,1,l=5),cex = 2)
rasterImage(legend_image, 0, 0, 1,1)



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1% - 10% choosen site", '11% - 50% choosen site', '50% - 100% choosen site',
                            'Surviving sites', 'Not surviving sites'), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('grey60', 'grey40', 'black', 'blue', 'orange'))
mtext("Frequency of choices", at=0.2, cex=2)



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend=c('Suitable conditions lost',
                             'Suitable conditions maintained', 
                             'New suitable conditions',
                             'Maintained range',
                             'Lost range'),
       pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('grey40', 'grey50', 'gray90', 'blue', 'red'))



plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Suitable conditions lost',
                            'Suitable conditions maintained', 
                            'New suitable conditions',
                            'Initially present sites'),
       pch=16, 
       pt.cex=3, 
       cex=1.5, 
       bty='n',
       col = c('grey40', 'grey50', 'gray60', 'green3'))

mtext("Species", at=0.2, cex=2)

nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
data <- expand.grid(X=nr, Y=nc)
cccc = 600-cost3
data$IntroductionTime <- c(as.matrix( (0 + 30*(cccc-min(cccc,na.rm=TRUE)) / (max(cccc,na.rm=TRUE)-min(cccc,na.rm=TRUE)))))
ggplot() + 
  coord_fixed()+
  scale_fill_continuous(na.value = "white",low="yellow", high="brown")+
  geom_raster(data = data , aes(x = X, y=Y,fill = IntroductionTime )) + 
  theme(panel.background = element_blank(),
        strip.background = element_blank() ,
        line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


countmap = array(0,c(nrow,ncol))


nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
data <- expand.grid(X=nr, Y=nc)
cccc = 600-cost3
data$IntroductionTime <- c(as.matrix( (0 + 30*(cccc-min(cccc,na.rm=TRUE)) / (max(cccc,na.rm=TRUE)-min(cccc,na.rm=TRUE)))))
ggplot() + 
  coord_fixed()+
  scale_fill_continuous(na.value = "white",low="yellow", high="brown")+
  geom_raster(data = data , aes(x = X, y=Y,fill = IntroductionTime )) + 
  theme(panel.background = element_blank(),
        strip.background = element_blank() ,
        line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

