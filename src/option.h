#include <RcppEigen.h>
//#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
//#include <Rcpp.h>
#include <numeric>    
#include <algorithm>    
#include <bits/stdc++.h> 
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <cstdlib>
#include <chrono>
#include <thread>
#include "maths_tools.h";
#include "local_optimising_planting_choice.h";
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

Eigen::SparseMatrix<double> rcpp_global_suitable_sites(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix);

Rcpp::NumericMatrix rcpp_global_suitable_coordinates(Eigen::SparseMatrix<double>& globalSuitableSites);

Eigen::SparseMatrix<double> rcpp_local_transition_matrix(Eigen::SparseMatrix<double> globalSuitableSites,Rcpp::NumericMatrix globalSuitableCoordinates,Rcpp::NumericMatrix migrationKernel);

std::list<Eigen::SparseMatrix<double>> rcpp_transition_matrices(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix,Eigen::SparseMatrix<double> localTransitionMatrix,Rcpp::NumericMatrix globalSuitableCoordinates);

std::vector<Eigen::SparseMatrix<double>> rcpp_colonisation_matrices(std::list<Eigen::SparseMatrix<double>> transitionMatrices);

Eigen::SparseMatrix<double> rcpp_viable_sites(std::vector<Eigen::SparseMatrix<double>> colonisationMatrices);

/*
Eigen::SparseVector<double> rcpp_get_current_vector(Eigen::SparseMatrix<double> currentPresenceMatrix,
                                                    std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                                    Eigen::SparseMatrix<double> globalSuitableSites);*/

Rcpp::NumericMatrix rcpp_viable_triplets(Eigen::SparseMatrix<double> viableSites,
                                         std::list<Eigen::SparseMatrix<double>> colonisationMatrices,
                                         Rcpp::NumericMatrix globalSuitableCoordinates,
                                         Eigen::SparseMatrix<double> globalSuitableSites,
                                         Eigen::SparseMatrix<double> costMatrix);

Eigen::SparseMatrix<double> rcpp_viable_values(Rcpp::NumericMatrix viablesTriplets,
                                               Eigen::SparseMatrix<double> viableSites,
                                               Eigen::SparseMatrix<double> globalSuitableSites,
                                               std::vector<Eigen::SparseMatrix<double>> colonisationMatrices);

NumericVector rcpp_eval_current_prob(int threshold,
                                     Eigen::SparseMatrix<double> currentPresenceMatrix,
                                     std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                     Eigen::SparseMatrix<double> globalSuitableSites);

Rcpp::NumericVector rcpp_pheromons(Rcpp::NumericMatrix viablesTriplets);
  

  
Rcpp::NumericVector rcpp_generate_population(Rcpp::NumericVector pheromons,Eigen::SparseMatrix<double> globalSuitableSites,int npop, int nbtoplant);

Rcpp::NumericVector rcpp_algorithm_opt(Rcpp::NumericVector pheromons,
                                       Rcpp::NumericMatrix viablesTriplets,
                                       Rcpp::NumericMatrix population0,
                                       Eigen::SparseMatrix<double> costMatrix,
                                       Eigen::SparseMatrix<double> currentPresenceMatrix,
                                       std::vector<Eigen::SparseMatrix<double>> colonisationMatrices,
                                       Eigen::SparseMatrix<double> globalSuitableSites,
                                       Eigen::SparseMatrix<double> viablesValues,
                                       int threshold,
                                       double confidence,
                                       int npop,
                                       int nsur,
                                       int ngen,
                                       int nbtoplant);





Rcpp::NumericMatrix rcpp_result_to_choice(Rcpp::NumericMatrix lastPopulation,Rcpp::NumericMatrix viablesTriplets);






//arma::sp_mat global_suitable_sites(std::list<arma::sp_mat> consecutiveSuitabilityMatrix);
