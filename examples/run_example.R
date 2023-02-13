graphics.off()

#################################################
#################################################
### Definitive example script
#################################################
#################################################



### Imports
library("usethis")
library("roxygen2")
library("devtools")
library(Rcpp)
library(RcppEigen)
library(methods)
library(lattice)
library(ggplot2)
library(dplyr)
library(purrr)
library(Matrix)
library("viridis")
library("IESTR")
library(hrbrthemes)
library(grid)
library(raster)
source("examples/building_example_functions.R")
source("examples/plot_exemple_functions.R")
### Global Details
col = c('green3','green2','blue1','blue3')
col = terrain.colors(20)
hea = (viridis(10))

#################################################
## Initialisation & Variables
#################################################

set.seed(17)
rcpp_set_seed(17)

# Map - size
nrow = 50
ncol = 50
nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)

# Map - time caracteristics
N_cycles = 30
Tadd_cyc = 1.5

# Map - space caracteristics
Tmin_alt_0 = 10
Tmax_alt_0 = 15
Tadd_alt_1 = -8
Trange = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)
Tmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

# Algorithm - genetic algo caracteristics
npop = 500
nsur = 100
ngen = 30

# Algorithm - optimum condition values
threshold = 50
confidence = 0.9

# Species - caracteristics
Trange_spe = c(8.8,9.1,11.6,11.9)
migr_spe = array(0.01, c(3, 3))
migr_spe[2,2] = 1

# Species - presence - area 1
nb_cell_occuped1 = 10
xmin = 0.15
xmax = 0.34
ymin = 0.2
ymax = 0.34
lim1 = c(xmin,xmax,ymin,ymax)

# Species - presence - area 2
nb_cell_occuped2 = 10
xmin = 0.72
xmax = 0.9
ymin = 0.72
ymax = 0.9
lim2 = c(xmin,xmax,ymin,ymax)

# Species - presence - area 3
nb_cell_occuped3 = 10
xmin = 0.05
xmax = 0.21
ymin = 0.25
ymax = 0.45
lim3 = c(xmin,xmax,ymin,ymax)

#################################################
## Maps construction
#################################################

# Common for all the figures
height_map = height_map(nrow,ncol)
climate_maps = climat_maps_maker(height_map,Tadd_cyc,N_cycles)
suitability_maps = suitability_maps(climate_maps,Trange_spe)
presence_1st_area = presence_map(nrow,ncol,suitability_maps,height_map,lim2,nb_cell_occuped1)
#presence_2nd_area = presence_map(nrow,ncol,suitability_maps,height_map,lim2,nb_cell_occuped2)
#presence_3nd_area = presence_map(nrow,ncol,suitability_maps,height_map,lim3,nb_cell_occuped3)
presence_4nd_area = presence_map_2(nrow,ncol,suitability_maps,height_map,N_cycles,40)

presence_map = presence_1st_area + presence_4nd_area



presence_map[is.na(presence_map)] = 0
presence_map[0!=(presence_map)] = 1
presence_map = Matrix(presence_map, sparse = T)
presence_map = as(presence_map,"dgCMatrix")[nrow:1,]

cost1 = cost_map(nrow,ncol,height_map,cost_mapping1)
cost2 = cost_map(nrow,ncol,height_map,cost_mapping2)
cost3 = cost_map(nrow,ncol,height_map,cost_mapping3)

#################################################
## Maps plot
#################################################

do_plot_fig2(nr, nc, height_map, climate_maps, N_cycles)
#do_plot_fig3(nr, nc, cost1,cost2, cost3)

#################################################
## Results Standards
#################################################

gss = rcpp_global_suitable_sites(suitability_maps)
gsc = rcpp_global_suitable_coordinates(gss)
ltm = rcpp_local_transition_matrix(gss,gsc,migr_spe)
tm = rcpp_transition_matrices(suitability_maps,ltm,gsc)
cm = rcpp_colonisation_matrices(tm)




vs = rcpp_viable_sites(cm)
vt = rcpp_viable_triplets(vs,cm,gsc,gss,cost1)
vv = rcpp_viable_values(vt,vs,gss,cm)
#cv = rcpp_get_current_vector(pres,cm,gss)

do_plot_fig3(nrow,ncol,N_cycles,height_map,suitability_maps,gsc,presence_map)


ph = rcpp_pheromons(vt)
ecp = rcpp_eval_current_prob(threshold,presence_map,cm,gss)
ntp = threshold - which(cumsum(ecp)>0.95)[1] 

po = rcpp_generate_population(ph,gss,npop,ntp)
resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map,cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
choix1 = rcpp_result_to_choice(resultat1,vt)

do_plot_fig4(nrow,ncol,N_cycles,height_map,suitability_maps,
             global_suitable_coordinates,colonisation_matrices,presence_map)

do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix1,
             global_suitable_coordinates,colonisation_matrices,presence_map)


#################################################
## Resultant plot
#################################################

vt3 = rcpp_viable_triplets(vs,cm,gsc,gss,cost3)
vv3 = rcpp_viable_values(vt3,vs,gss,cm)
ph3 = rcpp_pheromons(vt3)
po3 = rcpp_generate_population(ph3,gss,npop,ntp)
resultat3 = rcpp_algorithm_opt(ph3,vt3,po3,cost3,presence_map,cm,gss,vv3,threshold,confidence,npop,nsur,ngen,ntp)
choix3 = rcpp_result_to_choice(resultat3,vt3)
do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix3,
             global_suitable_coordinates,colonisation_matrices,presence_map)

#################################################
## Run the program over and over
#################################################

choix = list()

NN = 10

for(i in 1:NN){
  message(i)
  resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map,cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
  choix1 = rcpp_result_to_choice(resultat1,vt)
  choix[[i]] = choix1
}

countmap3 = array(0,c(nrow,ncol))
for (i in 1:NN){
  choice = choix[[i]]
  tmp_c= choice[choice[,1]!=0,]
  
  for (j in 1:length(tmp_c[,1])){
    countmap3[tmp_c[j,1],tmp_c[j,2]] = countmap3[tmp_c[j,1],tmp_c[j,2]] +1
  }
}


do_plot_fig6(nrow,ncol,N_cycles,height_map,suitability_maps,gsc,presence_map,choix)

choix10 = list()

NN = 200

for(i in 1:NN){
  message(i)
  resultat1 = rcpp_algorithm_opt(ph,vt,po,cost3,presence_map,cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
  choix1 = rcpp_result_to_choice(resultat1,vt)
  choix10[[i]] = choix1
}

countmap33 = array(0,c(nrow,ncol))
for (i in 1:NN){
  choice = choix10[[i]]
  tmp_c= choice[choice[,1]!=0,]
  
  for (j in 1:length(tmp_c[,1])){
    countmap33[tmp_c[j,1],tmp_c[j,2]] = countmap33[tmp_c[j,1],tmp_c[j,2]] +1
  }
}


do_plot_fig6(nrow,ncol,N_cycles,height_map,suitability_maps,gsc,presence_map,choix10)
























nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)
data <- expand.grid(X=nr, Y=nc)
cccc = countmap3 + (cost3*0)
data$IntroductionTime <- c(as.matrix( cccc ))
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









npops = c(50,100,200,400,600,800)
nsurs = c(0.02,0.04,0.08,0.16,0.24,0.48)

resultat_comparaison_cost = Matrix(0,6,6)
resultat_comparaison_time = Matrix(0,6,6)

for (i in 1:6){
  for (j in 1:6){
    
    message("Simulation; i=",i," j=",j)
    
    start.time <- Sys.time()
    
    npop = npops[i]
    nsur = nsurs[j] * npop
    
    po = rcpp_generate_population(ph,gss,npop,ntp)
    resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map,cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
    choix1 = rcpp_result_to_choice(resultat1,vt)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    
    resultat_comparaison_cost[i,j]=sum(choix1[,4])
    resultat_comparaison_time[i,j]=as.numeric(time.taken)
  }
}



choix = list()

# Algorithm - genetic algo caracteristics
npop = 300
nsur = 100
ngen = 30

ph = rcpp_pheromons(vt)
ecp = rcpp_eval_current_prob(threshold,presence_map,cm,gss)
ntp = threshold - which(cumsum(ecp)>0.95)[1] 

po = rcpp_generate_population(ph,gss,npop,ntp)

for(i in 1:1000){
  resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map,cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
  choix1 = rcpp_result_to_choice(resultat1,vt)
  choix[[i]] = choix1
}









