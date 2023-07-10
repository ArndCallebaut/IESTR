


#include <RcppEigen.h>
#include <iostream>
#include <vector>
#include <Rcpp.h>
#include <numeric>      
#include <algorithm>    
#include <bits/stdc++.h> 
#include <iostream>
#include <stdlib.h> 
#include <string>
#include <cstdlib>
#include "maths_tools.h";

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(Matrix)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

//' (obsolete) test if removing each site can allow to still reach the threshold with a given confidence
//' then choose randomly to remove it from the individual.
//'
//' @param viablesValues2 matrix of effects of introduction of each choice.
//' @param threshold number of sites colonised desired.
//' @param confidence matrix of effects of introduction of each choice.
//' @param viablesTriplets link between sites&time and viablesValues2.
//' @param population total population of the current generation in the evolutionnary algorithm.
//' @param pop index of the studied individual of the population
//' @param current current probabilities of presence of the species at the end of the cycles, with the unchanged individual. 
int half_random_add_site(Eigen::SparseMatrix<double> viablesValues2, 
                             int threshold, 
                             double confidence, 
                             Rcpp::NumericMatrix viablesTriplets, 
                             Rcpp::NumericMatrix population, 
                             int pop,  
                             Eigen::SparseVector<double> current){
  
  int n = population.cols();
  if (n==0){
    return(-1);
  }
  int elimination_index = -1;
  double gain = 0;
  for (int s=0 ; s < n ; ++s){
    if(population(pop,s)!=-1){
      
      Eigen::SparseVector<double> cur = current;
      for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,s)); it; ++it) {
        cur.coeffRef(it.index()-1) = (cur.coeffRef(it.index()-1)-it.value())/(1-it.value());
      }
      Rcpp::NumericVector vec = eval_probabilityVector(cur,threshold);
      if (vec(threshold)>=confidence){
        if(gain<viablesTriplets(population(pop,s),3)){
          elimination_index = s;
          gain = viablesTriplets(population(pop,s),3);
          break;
        }
      }
    }
  }
  if (elimination_index ==-1){
    return(-1);
  }
  if (elimination_index >=0){
    for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,elimination_index)); it; ++it) {
      current.coeffRef(it.index()-1) = (current.coeffRef(it.index()-1)-it.value())/(1-it.value());
    }
    population(pop,elimination_index) = -1;
    return(elimination_index);
  }
  return(-1);
}

//' test if removing each site can allow to still reach the threshold with a given confidence
//' then choose the more expensive choice to remove it from the individual.
//'
//' @param viablesValues2 matrix of effects of introduction of each choice.
//' @param threshold number of sites colonised desired.
//' @param confidence matrix of effects of introduction of each choice.
//' @param viablesTriplets link between sites&time and viablesValues2.
//' @param population total population of the current generation in the evolutionnary algorithm.
//' @param pop index of the studied individual of the population
//' @param current current probabilities of presence of the species at the end of the cycles, with the unchanged individual. 
int optimise_remove_site(Eigen::SparseMatrix<double> viablesValues2, 
                         int threshold, 
                         double confidence, 
                         Rcpp::NumericMatrix viablesTriplets, 
                         Rcpp::NumericMatrix population, 
                         int pop, 
                         Eigen::SparseVector<double> current){
  // We want to remove a site
  int n = population.cols();
  if (n==0){
    return(-1);
  }
  int elimination_index = -1;
  double gain = 0;
  for (int s=0 ; s < n ; ++s){
    if(population(pop,s)!=-1){
      Eigen::SparseVector<double> cur = current;
      for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,s)); it; ++it) {
        cur.coeffRef(it.index()-1) = (cur.coeffRef(it.index()-1)-it.value())/(1-it.value());
      }
      Rcpp::NumericVector vec = eval_probabilityVector(cur,threshold);
      if (vec(threshold)>=confidence){
        if(gain<viablesTriplets(population(pop,s),3)){
          elimination_index = s;
          gain = viablesTriplets(population(pop,s),3);
          break;
        }
      }
    }
  }
  
  if (elimination_index ==-1){
    return(-1);
  }
  
  if (elimination_index >=0){
    for (Eigen::SparseMatrix<double>::InnerIterator it(viablesValues2,population(pop,elimination_index)); it; ++it) {
      current.coeffRef(it.index()-1) = (current.coeffRef(it.index()-1)-it.value())/(1-it.value());
    }
    population(pop,elimination_index) = -1;
    return(elimination_index);
  }
  return(-1);
}

//' test if adding a site can allow to reach the threshold with a given confidence
//' then choose the less expensive choice.
//'
//' @param viablesValues2 matrix of effects of introduction of each choice.
//' @param threshold number of sites colonised desired.
//' @param confidence desired confidence of reaching a certain number of sites
//' @param viablesTriplets link between sites&time and viablesValues2.
//' @param population total population of the current generation in the evolutionnary algorithm.
//' @param pop index of the studied individual of the population
//' @param current current probabilities of presence of the species at the end of the cycles, with the unchanged individual. 
void optimise_add_site(Eigen::SparseMatrix<double> viablesValues2, 
                       int threshold, 
                       double confidence, 
                       Rcpp::NumericMatrix viablesTriplets, 
                       Rcpp::NumericMatrix population, 
                       int pop, 
                       Eigen::SparseVector<double> current){
  // We want to add a site
  int n1 = viablesValues2.rows();
  if (n1==0){
    return;
  }
  int n2 = population.cols();
  int added_index = -1;
  double cost = 1000000;
  int s;
  for (s=0 ; s < n2 ; ++s){
    if(population(pop,s)==-1){
      for(int m=0; m<n1;++m){
        Rcpp::NumericVector vec = eval_probabilityVector_adding(current,viablesValues2.col(m),threshold);
        if (vec(threshold)>=confidence){
          if(cost>viablesTriplets(m,4)){
            added_index = m;
            cost = viablesTriplets(m,4);
            break;
          }
        }
      }
      break;
    }
  }
  if (added_index ==-1){
    return;
  }

  if (added_index >=0){
    population(pop,s) = added_index;
    return;
  }
  return;
}