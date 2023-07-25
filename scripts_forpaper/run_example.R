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

### call for support scripts related to run examples and plot results 
source("./scripts_forpaper/building_example_functions.R")
source("./scripts_forpaper/plot_exemple_functions.R")

### Global color defaults
col = terrain.colors(20)
hea = (viridis(10))

#################################################
## Initialisation & Variable
#################################################

# Setting seed for R and C++
set.seed(17)
rcpp_set_seed(17)

### Algorithm - optimum condition values
#Threshold = Goal number of sites of species presence at the end of the period
threshold = 50 
#Confidence = Minimum probability that your goal is achieved during the optimization process 
#(the lower the confidnece the bigger the risk, related to flexibility in the solutions)
confidence = 0.8 

#################################################
## Building suitability maps 
#################################################

### Map - size
nrow = 50
ncol = 50
nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)

#### Map - time characteristics
#N time steps
N_cycles = 30
#Temp increase at the end of the time period
Tadd_cyc = 1.5

#### Map - space characteristics
Tmin_alt_0 = 10
Tmax_alt_0 = 15
Tadd_alt_1 = -8
Trange = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)
ltmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

### Species - characteristics
#Temperature ranges
Trange_spe = c(8.8,9.1,11.6,11.9)

# Species - presence 
nb_cell_occuped1 = 10
xmin = 0.72
xmax = 0.9
ymin = 0.72
ymax = 0.9
lim2 = c(xmin,xmax,ymin,ymax)

# Building maps used for the case study
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


#Dispersal and establishment probability kernel;
migr_spe = array(0.01, c(3, 3))
migr_spe[2,2] = 1


#################################################
## IESTR - case n°1 : Uniform costs 
#################################################

# All of this section describes how IESTR is used, step by step by step.
# There are two main results : 
# tm = the transition matrix, that allows to know the impact on the final state of the system of any introduction
# choix1 = the calculated optimal introductions in space and time.
# These can be directly calculated from the "IESTR_process" function, a wrapper for the different steps shown here.

#STEP 1 Internal indexing of sites used by the package; ====
#These build support matrices to speed up information; 
gss = rcpp_global_suitable_sites(suitability_maps)
gsc = rcpp_global_suitable_coordinates(gss)
# inital presence of the species
do_plot_fig3(nrow,ncol,N_cycles,height_map,suitability_maps,gsc,presence_map) 

#STEP 2 Creation of the transition matrices ====
#get the spread matrix [to get at the migration potential of sites]
sm = rcpp_spread_matrix(gss,gsc,migr_spe)
#obtain local transition matrix
ltm = rcpp_local_transition_matrix(suitability_maps,sm,gsc)
#obtain transition matrix
tm = rcpp_transition_matrix(ltm)

#STEP 3 Reindexing of values of the transition matrices ====
#simplifying to potential solutions
vs = rcpp_viable_sites(tm)
vt = rcpp_viable_triplets(vs,tm,gsc,gss,cost1)
vv = rcpp_viable_values(vt,vs,gss,tm)

#STEP 4 Pre-processing for the optimization part ====

#Weight for each viable time-site combinations for species introduction 
#  (more weight = more chances to be chosen)
ph = rcpp_pheromons(vt) 
#occurrence probability without introduction (suitability and migration only) to meet the species range goal (threshold) 
ecp = rcpp_eval_current_prob(threshold,presence_map,tm,gss) 
#Maximum number of sites of introduction necessary to reach the threshold.
ntp = threshold - which(cumsum(ecp)>0.95)[1] 

#Algorithm - genetic algo characteristics
npop = 300 #size of the initial population
nsur = 100 #size of the surviving population
ngen = 20 #number of generations

#generate the initial population for the genetic algorithm.
po = rcpp_generate_population(ph,gss,npop,ntp) 

#STEP 5 Genetic algorithm use to optimise introduction choices ====
resultat1 = rcpp_algorithm_opt(pheromons = ph,
                               viablesTriplets = vt,
                               population0 = po,
                               costMatrix = cost1,
                               currentPresenceMatrix = presence_map,
                               transitionmatrices = tm,
                               globalSuitableSites = gss,
                               viablesValues = vv,
                               threshold = threshold,
                               confidence = confidence,
                               npop = npop,
                               nsur = nsur,
                               ngen = ngen,
                               nbtoplant = ntp)

#STEP 6 Results filtered ====
choix1 = rcpp_result_to_choice(lastPopulation = resultat1,
                               viablesTriplets = vt)

#STEP 7 Show how would evolve the system without introduction ====
do_plot_fig4(nrow,ncol,N_cycles,height_map,suitability_maps,
             gsc,tm,presence_map)

#STEP 8 Show how would evolve the system with the introduction ====
do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix1,
             gsc,tm,presence_map)

#ARNAUD
#add IESTR_process example

#################################################
## IESTR - case n°2 : Non-uniform costs 
#################################################

vt3 = rcpp_viable_triplets(vs,tm,gsc,gss,cost3)
vv3 = rcpp_viable_values(vt3,vs,gss,tm)
ph3 = rcpp_pheromons(vt3)
po3 = rcpp_generate_population(ph3,gss,npop,ntp)
resultat3 = rcpp_algorithm_opt(ph3,vt3,po3,cost3,presence_map,tm,
                               gss,vv3,threshold,
                               confidence,
                               npop,nsur,ngen,ntp)
choix3 = rcpp_result_to_choice(resultat3,vt3)

#Show how would evolve the system with the introduction
do_plot_fig5(nrow,ncol,N_cycles,height_map,suitability_maps,choix3,
             global_suitable_coordinates,tm,presence_map)

