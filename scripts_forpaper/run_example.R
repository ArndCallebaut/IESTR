graphics.off()

#################################################
#################################################
### Definitive example script
#################################################
#################################################

### Imports
library(IESTR)
library(dplyr)
library(purrr)
library(Matrix)
library(grid)
library(raster)
library(lattice)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(methods)

### call for support scripts.
source("./scripts_forpaper/building_example_functions.R")
source("./scripts_forpaper/plot_exemple_functions.R")

### Global Details
col = terrain.colors(20)
hea = (viridis(10))

#################################################
## Initialisation & Variable
#################################################

# Setting seed for R and C++
set.seed(17)
rcpp_set_seed(17)

# Algorithm - genetic algo caracteristics
npop = 300 #size of the initial population
nsur = 100 #size of the surviving population
ngen = 20 #number of generations

# Algorithm - optimum condition values
threshold = 50
confidence = 0.8

#################################################
## Maps construction
#################################################

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
ltmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

# Species - caracteristics
Trange_spe = c(8.8,9.1,11.6,11.9)
migr_spe = array(0.01, c(3, 3))
migr_spe[2,2] = 1

# Species - presence - area 2
nb_cell_occuped1 = 10
xmin = 0.72
xmax = 0.9
ymin = 0.72
ymax = 0.9
lim2 = c(xmin,xmax,ymin,ymax)

# Construction of the maps used for the case study
height_map = height_map(nrow,ncol)
climate_maps = climat_maps_maker(height_map,Tadd_cyc,N_cycles)
suitability_maps = suitability_maps(climate_maps,Trange_spe)
presence_1st_area = presence_map_2(nrow,ncol,suitability_maps,height_map,N_cycles, 40) # find 40 sites in suitable sites that will not be suitable in the futur
presence_2nd_area = presence_map(nrow,ncol,suitability_maps,height_map,lim2,nb_cell_occuped1) 
presence_map = presence_concatenate(presence_1st_area ,presence_2nd_area,nrow)
cost1 = cost_map(nrow,ncol,height_map,cost_mapping1)
cost3 = cost_map(nrow,ncol,height_map,cost_mapping3)

# A plot of the height map & temperatures at the beginning/at the end
do_plot_fig2(nr, nc, height_map, climate_maps, N_cycles)

#################################################
## IESTR - case n°1 : Uniform costs 
#################################################

# All of this part describes how IESTR is used, step by step by step.
# There are to main result : 
# tm = the transition matrix, that allow to know the impact on the final state of the system of any introduction
# choix1 = the calculated optimal introductions in space and time.
# These can be directly calculated from the "IESTR_process" function for a quickier use of IESTR functionnalities.

#Utility indices of suitable sites used in the package; 
gss = rcpp_global_suitable_sites(suitability_maps)
gsc = rcpp_global_suitable_coordinates(gss)

do_plot_fig3(nrow,ncol,N_cycles,height_map,suitability_maps,gsc,presence_map) # inital presence of the species

#Creation of the transition matrices
#If you have more informations about spread dynamics than just a migration Kernel, this can be included here. (Contact the author for more infos)
sm = rcpp_spread_matrix(gss,gsc,migr_spe)
ltm = rcpp_local_transition_matrix(suitability_maps,sm,gsc)
tm = rcpp_transition_matrix(ltm)

#Utility reindicing of values of the transition matrices
vs = rcpp_viable_sites(tm)
vt = rcpp_viable_triplets(vs,tm,gsc,gss,cost1)
vv = rcpp_viable_values(vt,vs,gss,tm)

#Pre-processing for the optimisation part
ph = rcpp_pheromons(vt) #Weight for each viable pair of time&site of introduction (more weight = more chance to be choosen)
ecp = rcpp_eval_current_prob(threshold,presence_map,tm,gss) #probability of number of site with presence at the end of the period.
ntp = threshold - which(cumsum(ecp)>0.95)[1] #evaluated maximum number of sites of introduction necessary to reach the threshold.
po = rcpp_generate_population(ph,gss,npop,ntp) #generate the initial population for the genetic algorithm.

#Genetic algorithm use to optimise introduction choices
resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map,tm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)

#Results filtered
choix1 = rcpp_result_to_choice(resultat1,vt)

#Show how would evolve the system without introduction
do_plot_fig4(nrow,ncol,N_cycles,height_map,suitability_maps,
             gsc,tm,presence_map)

#Show how would evolve the system with the introduction
do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix1,
             gsc,tm,presence_map)


#################################################
## IESTR - case n°2 : Non-uniform costs 
#################################################

vt3 = rcpp_viable_triplets(vs,tm,gsc,gss,cost3)
vv3 = rcpp_viable_values(vt3,vs,gss,tm)
ph3 = rcpp_pheromons(vt3)
po3 = rcpp_generate_population(ph3,gss,npop,ntp)
resultat3 = rcpp_algorithm_opt(ph3,vt3,po3,cost3,presence_map,tm,gss,vv3,threshold,confidence,npop,nsur,ngen,ntp)
choix3 = rcpp_result_to_choice(resultat3,vt3)

#Show how would evolve the system with the introduction
do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix3,
             global_suitable_coordinates,tm,presence_map)

