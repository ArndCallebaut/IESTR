#' IESTR_process
#' 
#' Whole IESTR process
#'
#' @param suitability_maps list of matrix of suitability for each timestep.
#' @param presence_map sparse matrix of presence of the species.
#' @param migr_spe matrix of probability of spread in neighbour cells at each time step.
#' @param cost Cost of introduction of the species in each cell.
#' @param threshold Number of sites with presence is wanted at the end of the study period
#' @param confidence Probability that the threshold is reached is wanted (must be <1)
#' @param npop Genetic algorithm : size of the initial population
#' @param nsur Genetic algorithm : size of the survival population at each generation
#' @param ngen Genetic algorithm : number of generation
#' @return the transition matrices & the choices of introductions
#' @export
IESTR_process = function(suitability_maps,
                         presence_map,
                         migr_spe,
                         cost,
                         threshold,
                         confidence,
                         npop,
                         nsur,
                         ngen){
  
  nrow = length(presence_map[,1])
  ncol = length(presence_map[1,])
  N_cycles = length(suitability_maps)
  
  main_result = list()
  
  #Utility indices of suitable sites used in the package; 
  gss = rcpp_global_suitable_sites(suitability_maps)
  gsc = rcpp_global_suitable_coordinates(gss)
  
  #Creation of the transition matrices
  #If you have more informations about spread dynamics than just a migration Kernel, this can be included here. (Contact the author for more infos)
  sm = rcpp_spread_matrix(gss,gsc,migr_spe)
  ltm = rcpp_local_transition_matrix(suitability_maps,sm,gsc)
  tm = rcpp_transition_matrix(ltm)
  
  main_result[[1]] = tm
  
  #Utility reindicing of values of the transition matrices
  vs = rcpp_viable_sites(tm)
  vt = rcpp_viable_triplets(vs,tm,gsc,gss,cost)
  vv = rcpp_viable_values(vt,vs,gss,tm)
  
  #Pre-processing for the optimisation part
  ph = rcpp_pheromons(vt) #Weight for each viable pair of time&site of introduction (more weight = more chance to be choosen)
  ecp = rcpp_eval_current_prob(threshold,presence_map,tm,gss) #probability of number of site with presence at the end of the period.
  ntp = threshold - which(cumsum(ecp)>0.95)[1] #evaluated maximum number of sites of introduction necessary to reach the threshold.
  po = rcpp_generate_population(ph,gss,npop,ntp) #generate the initial population for the genetic algorithm.
  
  #Genetic algorithm use to optimise introduction choices
  resultat1 = rcpp_algorithm_opt2(ph,vt,po,cost,presence_map,tm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
  
  #Results filtered
  choix1 = rcpp_result_to_choice(resultat1,vt)
  
  main_result[[2]] = choix1
  
  return(main_result)
}
