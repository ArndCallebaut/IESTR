

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
hea = viridis(10)

#################################################
## Initialisation & Variable
#################################################

# Setting seed for R and C++
set.seed(98)
rcpp_set_seed(98)

### Algorithm - optimum condition values
#Threshold = Goal number3 of sites of species presence at the end of the period
threshold = 50 
#Confidence = Minimum probability that your goal is achieved during the optimization process 
#(the lower the confidnece the bigger the risk, related to flexibility in the solutions)
confidence = 0.9

#################################################
## Building suitability maps of the case study
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

vt3 = rcpp_viable_triplets(vs,tm,gsc,gss,cost3)
vv3 = rcpp_viable_values(vt3,vs,gss,tm)
ph3 = rcpp_pheromons(vt3)

#STEP 4 Pre-processing for the optimization part ====

#Weight for each viable time-site combinations for species introduction 
#  (more weight = more chances to be chosen)
ph = rcpp_pheromons(vt) 
#occurrence probability without introduction (suitability and migration only) to meet the species range goal (threshold) 
ecp = rcpp_eval_current_prob(threshold,presence_map,tm,gss) 
#Maximum number3 of sites of introduction necessary to reach the threshold.
ntp = 20

#Algorithm - genetic algo characteristics

npops = c(100,500,1000,2500,4000,10000,15000)
nsurs = c(10,50,100,250,400,1000,1500) #taux de survivants posé à 10%

ngen = 30
n = length(npops)
K = 20

number3 = matrix(0,K,n)
number32 = matrix(0,K,n)
choixx = list()
choixx2 = list()
m=1
for (j in 1:K){
  message("CALL FOR THE DAY - iteration ", j)
    for (i in 1:n){
      message(i)
      npop = npops[i]
      nsur = nsurs[i]
      
      po = rcpp_generate_population(ph,gss,npop,ntp)
      resultat = rcpp_algorithm_opt2(ph,vt,po,cost1,presence_map,tm,
                                     gss,vv,threshold,
                                     confidence,
                                     npop,nsur,ngen,ntp)
      
      po3 = rcpp_generate_population(ph3,gss,npop,ntp)
      resultat3 = rcpp_algorithm_opt2(ph3,vt3,po3,cost3,presence_map,tm,
                                      gss,vv3,threshold,
                                      confidence,
                                      npop,nsur,ngen,ntp)
      
      choix1 = rcpp_result_to_choice(resultat,vt)
      choix3 = rcpp_result_to_choice(resultat3,vt3)
      
      res = choix3[choix3[,4]!=0,]
      number3[j,i] = sum(res[,4])
      choixx[[m]]=choix3
      
      res = choix1[choix1[,4]!=0,]
      number32[j,i] = sum(res[,4])
      choixx2[[m]]=choix1
      
      m=m+1
    }
}
print(number3)



optimised_values = c(number3)
population_size= c(rep("  100",K),rep("  500",K),rep(" 1000",K),rep(" 2500",K),rep(" 4000",K),rep("10000",K),rep("15000",K))
tab = as.data.frame(cbind(population_size,optimised_values))
tab$category = as.factor(population_size)
tab$optimised_values = as.numeric(tab$optimised_values)
p <- ggplot(tab, aes(x=(population_size), y=optimised_values)) + 
  geom_boxplot()+theme(axis.text=element_text(size=25),
                       axis.title=element_text(size=25, face="bold"))
p

