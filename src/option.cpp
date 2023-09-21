#include <RcppEigen.h>
#include <iostream>
#include <vector>
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


using namespace std;
using namespace Rcpp;
using namespace Eigen;
using namespace std::chrono;

//' rcpp_global_suitable_sites.
//' 
//' Establishes from the consecutiveSuitabilityMatrix what are the
//' study sites in the map and attribute them an index used in 
//' other functions.
//' Warning : the indices start at 1.
//'
//' @param consecutiveSuitabilityMatrix (std::list<Eigen::SparseMatrix<double>>) list of matrix of suitability for each timestep.
//' @return suitable sites as a sparse matrix
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_global_suitable_sites(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix){
  int nbr = consecutiveSuitabilityMatrix.back().rows();
  int nbc = consecutiveSuitabilityMatrix.back().cols();
  Eigen::SparseMatrix<double> globalSuitableSites(nbr,nbc);
  Rcpp::NumericVector vec_nbpp(consecutiveSuitabilityMatrix.size());
  int k = 0;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    vec_nbpp(k) = thisSuitabilityMatrix.nonZeros();
    ++k;
  }
  int tnbpp = sum(vec_nbpp);
  Rcpp::NumericMatrix possibilities(tnbpp,5);
  int j=1;
  for (auto const& thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    for (k=0; k<thisSuitabilityMatrix.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(thisSuitabilityMatrix,k); it; ++it)
      {
        if(globalSuitableSites.coeffRef(it.row(),it.col())==0){
          globalSuitableSites.coeffRef(it.row(),it.col()) = double(j);
          j++;
        }
      }
    }
  globalSuitableSites.makeCompressed();
  return(globalSuitableSites);
}

//' rcpp_global_suitable_coordinates.
//' 
//' Calculate the coordinates of the suitable sites.
//'
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) support matrix containing index of each sites
//' @return coordinates of the suitable sites as a matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_global_suitable_coordinates(Eigen::SparseMatrix<double> globalSuitableSites){
  int nbsites = globalSuitableSites.nonZeros();
  Rcpp::NumericMatrix globalSuitableCoordinates(nbsites,3);
  for (int k=0; k<globalSuitableSites.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(globalSuitableSites,k); it; ++it)
    {
      globalSuitableCoordinates(it.value()-1,0)=it.row();
      globalSuitableCoordinates(it.value()-1,1)=it.col();
      globalSuitableCoordinates(it.value()-1,2)=int(it.value())-1;
    }
  return(globalSuitableCoordinates);
}

//' rcpp_spread_matrix
//' 
//' Calculate the spread_matrix.
//'
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) support matrix containing index of each sites
//' @param globalSuitableCoordinates (Rcpp::NumericMatrix) support matrix coordinates of each sites
//' @param migrationKernel (Rcpp::NumericMatrix) matrix modelling spread properties of the species.
//' @return Spread matrix used to model spread behavior of the species without suitability
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_spread_matrix(Eigen::SparseMatrix<double> globalSuitableSites,
                                               Rcpp::NumericMatrix globalSuitableCoordinates,
                                               Rcpp::NumericMatrix migrationKernel){
  int nbr = globalSuitableSites.rows();
  int nbc = globalSuitableSites.cols();
  int i;
  int j;  
  int migrationRange = (migrationKernel.rows()-1)/2;
  int nbsites = globalSuitableSites.nonZeros();
  Eigen::SparseMatrix<double> spreadmatrix(nbsites,nbsites);
  for (int index=0;index<nbsites;++index){
    i = globalSuitableCoordinates(index,0);
    j = globalSuitableCoordinates(index,1);
    for (int h=max(0,i-migrationRange); h<min(nbr,i+migrationRange+1); ++h){
      for (int k=max(0,j-migrationRange); k<min(nbc,j+migrationRange+1); ++k){
        if (globalSuitableSites.coeffRef(h,k)!=0 && migrationKernel(h+migrationRange-i,k+migrationRange-j)!=0){
          spreadmatrix.insert(globalSuitableSites.coeffRef(h,k)-1,index)=migrationKernel(h+migrationRange-i,k+migrationRange-j);
        }
      }
    }
  }
  spreadmatrix.pruned(0.01);
  spreadmatrix.makeCompressed();
  return (spreadmatrix);
}

//' rcpp_local_transition_matrix
//' 
//' Calculate the local transition matrices from the spread matrix and suitability matrices
//'
//' @param consecutiveSuitabilityMatrix (std::list<Eigen::SparseMatrix<double>>) list of suitability matrix
//' @param spreadmatrix (Eigen::SparseMatrix<double>) matrix used to model spread behavior of the species without suitability
//' @param globalSuitableCoordinates (Rcpp::NumericMatrix) support matrix coordinates of each sites
//' @return list of matrix used to model spread behavior of the species with suitability, between two adjacent timestep
//' @export
// [[Rcpp::export]]

std::list<Eigen::SparseMatrix<double>> rcpp_local_transition_matrix(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix,
                                                                Eigen::SparseMatrix<double> spreadmatrix,
                                                                Rcpp::NumericMatrix globalSuitableCoordinates){
  std::list<Eigen::SparseMatrix<double>> localtransitionmatrices;
  Eigen::SparseMatrix<double> localtransitionmatrix;
  int i=0;
  for (auto thisSuitabilityMatrix : consecutiveSuitabilityMatrix){
    localtransitionmatrix = spreadmatrix;
    for (int k=0; k<spreadmatrix.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(spreadmatrix,k); it; ++it)
      {
        localtransitionmatrix.coeffRef(it.row(),it.col()) *= thisSuitabilityMatrix.coeffRef(globalSuitableCoordinates(it.row(),0),globalSuitableCoordinates(it.row(),1));
      }
    }
    localtransitionmatrix.prune(0.01);
    ++i;
    localtransitionmatrix.makeCompressed();
    localtransitionmatrices.push_back(localtransitionmatrix);
  }
  localtransitionmatrices.reverse();
  return(localtransitionmatrices);
}

//' rcpp_transition_matrix
//' 
//' Calculate the transition matrices from the local transition matrices
//'
//' @param localtransitionmatrices (std::list<Eigen::SparseMatrix<double>>) list of matrix used to model spread behavior of the species with suitability, between two adjacent timestep
//' @return list of transition matrices from any timestep to the last timestep.
//' @export
// [[Rcpp::export]]
std::vector<Eigen::SparseMatrix<double>> rcpp_transition_matrix(std::list<Eigen::SparseMatrix<double>> localtransitionmatrices){
  std::vector<Eigen::SparseMatrix<double>> transitionmatrices;
  transitionmatrices.push_back(localtransitionmatrices.front());
  int k=1;
  bool first = true;
  for (auto const& elt : localtransitionmatrices) {
    if (first){first = false;}
    else{
      transitionmatrices.push_back(proba_matrix_mult3(transitionmatrices.back(),elt));
      ++k;
    }
  }
  //transitionmatrices.pruned(0.01);
  reverse(transitionmatrices.begin(),transitionmatrices.end());
  return(transitionmatrices);
}

//' rcpp_viable_sites
//' 
//' Calculate sites&times introduction considered "viable" 
//' eg that are not totally bested by other timings of introduction for the same site, 
//' and that are efficient enough. (colonises more than one site on average).
//'
//' @param transitionmatrices (Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @return boolean sparse matrix, each row correspond to a time, each column to a site. (1=viable; 0=not viable)
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_viable_sites(std::vector<Eigen::SparseMatrix<double>> transitionmatrices){
  int nbsites = transitionmatrices.front().rows();
  int nbperiod = transitionmatrices.size();
  Eigen::SparseMatrix<double> viableSites(nbperiod,nbsites);
  bool av1;
  bool av2;
  int h;
  int k;
  double var;
  Eigen::SparseVector<int> suit1;
  Eigen::SparseVector<int> suit2;
  bool nonzero1;
  bool nonzero2;
    
  for (int i=0;i<nbsites;++i){
    Rcpp::NumericVector viability (nbperiod);
    for (int j=0;j<nbperiod;++j){
      Eigen::SparseVector<int> suit2 = transitionmatrices[j].row(i);
      if(suit2.nonZeros()!=0){
        viability(j)=1;
      }
    }
    h=0;
    for (int ss=0; ss<(nbperiod-1); ++ss){
      h=ss;
      suit1 = transitionmatrices[ss].row(i);
      if ( viability(h)==0 ){
        break;
      }
      for (int sss=ss+1; sss<nbperiod; ++sss){
        k=sss;
        suit2 = transitionmatrices[sss].row(i);
        nonzero1 = false;
        nonzero2 = false;
        if (viability(k)==0 || viability(h)==0 ){
          break;
        }
        av1 = false;
        av2 = false;
        if (k>=h){
          ++k;
          break;
        }
        for (int j=0;j<nbsites;++j){
          if(suit1.coeffRef(j,i)!=0){
            nonzero1 = true;
          }
          if(suit2.coeffRef(j,i)!=0){
            nonzero2 = true;
          }
          var = suit1.coeffRef(j)-suit2.coeffRef(j);
          if(var>0){
            av1=true;
          }
          if(var<0){
            av2=true;
          }
          if(av1 && av2){
            break;
          }
        }
        if(!nonzero1){
          // si le premier site ne donne aucun resultat
          viability(h)=0;
        }
        if(!nonzero2){
          // si le second site ne donne aucun resultat
          viability(k)=0;
        }
        if (av1 && !av2){
          viability(k)=0;
        }
        if (av2 && !av1){
          viability(h)=0;
        }
        if (!av2 && !av1){
          viability(k)=0;
        }
        ++k;
      }
      ++h;
    }
    
    for (int j=0;j<nbperiod;++j){
      if (viability[j]==1){
        if((transitionmatrices.at(j).col(i)).sum()>=0.99){
          viableSites.coeffRef(int(j),int(i))=1;
          //viableSites2.coeffRef(int(i),int(j))=1;
        }
      }
    }
  }
  viableSites.makeCompressed();
  return(viableSites);
}

//' rcpp_viable_triplet
//' 
//' Returns matrix summing up information of each viable site&time pair.
//' col1&2 = X & Y coordinates
//' col3 = time of introduction
//' col4 = cost of introduction
//' col5 = index of introduction
//' col6 = mean number of final colonised sites due to introduction
//'
//' @param viableSites (Eigen::SparseMatrix<double>) boolean sparse matrix, each row correspond to a time, each column to a site. (1=viable; 0=not viable)
//' @param transitionmatrices (std::list<Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @param globalSuitableCoordinates (Rcpp::NumericMatrix) support matrix coordinates of each sites
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) support matrix containing index of each sites
//' @param costMatrix (Eigen::SparseMatrix<double>) matrix of cost of introduction
//' @return For each viable site : index, position on map, time, and cost.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_viable_triplets(Eigen::SparseMatrix<double> viableSites,
                                         std::list<Eigen::SparseMatrix<double>> transitionmatrices,
                                         Rcpp::NumericMatrix globalSuitableCoordinates,
                                         Eigen::SparseMatrix<double> globalSuitableSites,
                                         Eigen::SparseMatrix<double> costMatrix){
  int nboc = viableSites.nonZeros();
  Rcpp::NumericMatrix viablesTriplets(nboc,6);
  
  
  int t = 0;
  int index = 0;
  for (auto colma : transitionmatrices){
    Eigen::SparseVector<double> loc = viableSites.row(t);
    for (Eigen::SparseVector<double>::InnerIterator it(loc); it; ++it){
      viablesTriplets(index,0) = globalSuitableCoordinates(it.index(),0);
      viablesTriplets(index,1) = globalSuitableCoordinates(it.index(),1);
      viablesTriplets(index,2) = t;
      viablesTriplets(index,3) = costMatrix.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      viablesTriplets(index,4) = globalSuitableSites.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1));
      viablesTriplets(index,5) = colma.col(globalSuitableSites.coeffRef(viablesTriplets(index,0),viablesTriplets(index,1))).sum();
      ++index;
    }
    ++t;
  }
  return (viablesTriplets);
}

//' rcpp_viable_values
//' 
//' Summing up the effect on the final state of the system of each viable site&time pair.
//'
//' @param viablesTriplets (Rcpp::NumericMatrix) For each viable site : index, position on map, time, and cost.
//' @param viableSites (Eigen::SparseMatrix<double>) boolean sparse matrix, each row correspond to a time, each column to a site. (1=viable; 0=not viable)
//' @param globalSuitableSites (Eigen::SparseMatrix<double> ) support matrix containing index of each sites
//' @param transitionmatrices (std::vector<Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @return for each viable site : index, position on map, time, and cost.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> rcpp_viable_values(Rcpp::NumericMatrix viablesTriplets,
                                                Eigen::SparseMatrix<double> viableSites,
                                                Eigen::SparseMatrix<double> globalSuitableSites,
                                                std::vector<Eigen::SparseMatrix<double>> transitionmatrices){
  int nbsites = globalSuitableSites.nonZeros();
  int nboc = viableSites.nonZeros();
  int t = 0;
  int index = 0;
  Eigen::SparseMatrix<double> viableSites2 = viableSites.transpose();
  Eigen::SparseMatrix<double> viablesValues(nbsites,nboc);
  Eigen::SparseVector<double> colo2;
  
  for (int k=0; k<viableSites2.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(viableSites2,k); it; ++it)
    {
      if(it.value()==0){
        continue;
      }
      int xx = viablesTriplets(index,0);
      int yy = viablesTriplets(index,1);
      double val = globalSuitableSites.coeffRef(xx,yy) ;
      colo2 = transitionmatrices.at(it.col()).col(it.row());
      
      for (Eigen::SparseVector<double>::InnerIterator it2(colo2); it2; ++it2){
        viablesValues.insert(it2.index(),index) = it2.value();
      }
      ++index;
    }
    viablesValues.makeCompressed();
    return(viablesValues);
}

   


//' rcpp_eval_current_prob
//' 
//' Evaluate how would evolve the current present species without introduction.
//'
//' @param threshold (int) number of site with presence we want to obtain.
//' @param currentPresenceMatrix (Eigen::SparseMatrix<double>) matrix of presence of the species.
//' @param transitionmatrices (std::vector<Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) support matrix containing index of each sites
//' @return for each viable site : index, position on map, time, and cost.
//' @export
// [[Rcpp::export]]
NumericVector rcpp_eval_current_prob(int threshold,
                                     Eigen::SparseMatrix<double> currentPresenceMatrix,
                                     std::vector<Eigen::SparseMatrix<double>> transitionmatrices,
                                     Eigen::SparseMatrix<double> globalSuitableSites){
  
  int nbsites = (globalSuitableSites).nonZeros();
  int x;
  int y;
  int s;
  double v;
  Eigen::SparseVector<double> currentVector(nbsites);
  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      for (Eigen::SparseMatrix<double>::InnerIterator it2(transitionmatrices.back(),s-1); it2; ++it2) {
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
      }
    }
  NumericVector res = eval_probabilityVector(currentVector,threshold);
  return(res);
}

//' rcpp_pheromons
//' 
//' Evaluate the potential of each site&time pair, and give them weights for the optimization algorithm.
//'
//' @param viablesTriplets (Rcpp::NumericMatrix) information about each of the viable site&time pair.
//' @return a pre-ranking of each site&time pair for the algorithm.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pheromons(Rcpp::NumericMatrix viablesTriplets){
  int nboc = viablesTriplets.nrow();
  Rcpp::NumericVector pheromons(nboc);
  for (int i=0 ; i<nboc ; i++){
    pheromons[i] = viablesTriplets(i,5)/viablesTriplets(i,3);
  } 
  //pheromons = 1 / (1+pheromons-min(pheromons));
  pheromons = pheromons / sum(pheromons);
  return(pheromons);
}

//' rcpp_generate_population
//' 
//' Generate initial population for the genetic algorithm.
//'
//' @param pheromons (Rcpp::NumericVector) weights of relative importance of each site&time pair.
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) support matrix containing index of each sites
//' @param npop (int) number of individuals in the initial population
//' @param nbtoplant (int) pre-evaluate number of introduction necessary (overestimated, usually)
//' @return an initial population of set of site&time choices of introduction
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_generate_population(Rcpp::NumericVector pheromons,
                                             Eigen::SparseMatrix<double> globalSuitableSites,int npop, int nbtoplant){
  Rcpp::NumericMatrix population(npop,nbtoplant);
  Rcpp::NumericVector permutation(nbtoplant);
  for (int j=0; j<npop; ++j) {
    permutation=generate_permutation4(nbtoplant,pheromons);
    population(j,_) = permutation ;
  }
  return(population);
}


//' rcpp_algorithm_opt
//' 
//' (please use rcpp_algorithm_opt2, this version is slower) Genetic algorithm used to optimised choices of introduction under the constraint of reaching a dertain number of 
//' presence at the end of the study period (threshold) with a minimal probability (confidence) while miimizing the cost
//' of all the introductions
//'
//' @param pheromons (Rcpp::NumericVector) weights of relative importance of each site&time pair.
//' @param viablesTriplets (Rcpp::NumericMatrix) For each viable site : index, position on map, time, and cost.
//' @param population0 (Rcpp::NumericMatrix) An initial population of set of site&time choices of introduction
//' @param costMatrix (Eigen::SparseMatrix<double>) matrix of costs of introduction, sparse.
//' @param currentPresenceMatrix (Eigen::SparseMatrix<double>) matrix of current presence of the species, sparse.
//' @param transitionmatrices (std::vector<Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) matrix of viable sites to use for the optimiation part.
//' @param viablesValues (Eigen::SparseMatrix<double>) viables values from transitions matrices from viable sites.
//' @param threshold (int) number of presence site to reach at the end of the studied period, with the introductions
//' @param confidence (double) minimal probability to reach the threshold accepted by the optimisation algorithm 
//' @param npop (int) genetic algorithm : size of the population
//' @param nsur (int) genetic algorithm : size of the surviving population
//' @param ngen (int) genetic algorithm : number of generation
//' @param nbtoplant evaluated maximum number of introduction necessary to satisfy the constraints of final presence.
//' @return an optimised population matrix of set of site&time choices of introduction.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_algorithm_opt(Rcpp::NumericVector pheromons,
                                       Rcpp::NumericMatrix viablesTriplets,
                                       Rcpp::NumericMatrix population0,
                                       Eigen::SparseMatrix<double> costMatrix,
                                       Eigen::SparseMatrix<double> currentPresenceMatrix,
                                       std::vector<Eigen::SparseMatrix<double>> transitionmatrices,
                                       Eigen::SparseMatrix<double> globalSuitableSites,
                                       Eigen::SparseMatrix<double> viablesValues,
                                       int threshold,
                                       double confidence,
                                       int npop,
                                       int nsur,
                                       int ngen,
                                       int nbtoplant){
  
  // **************************************************************************
  // ** 0 INITIALISATION ******************************************************
  // **************************************************************************
  
  int index;    
  int g;
  int x;
  int y;
  int s;
  double v;
  double cost;
  int indexmini;
  double cout_mini;
  int papa;
  int mama;
  int f;
  double v1;
  double v2;
  int progress_bar_var;
  
  Rcpp::NumericMatrix population=population0;
  Rcpp::NumericMatrix survivors(npop,nbtoplant);
  Rcpp::NumericVector ranks(npop);
  
  // **************************************************************************
  // ** 1 EVALUATION FINAL STATE **********************************************
  // **************************************************************************
  
  int nboc = viablesTriplets.nrow();
  int nbsites = (globalSuitableSites).nonZeros();
  Eigen::SparseVector<double> currentVector(nbsites);
  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      
      for (Eigen::SparseMatrix<double>::InnerIterator it2(transitionmatrices.back(),s-1); it2; ++it2) {
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
      }
    }
    
  // **************************************************************************
  // ** 2 SELECTION LOOP PROCESS **********************************************
  // **************************************************************************
  
  if(g%10!=0){
    std::cout << "Error, ngen must be multiple of 10. ngen setted to 10." << std::endl;
    ngen = 10;
  }
  Rcpp::NumericMatrix evaluation(npop,2);
  int ref_10 = (int) ((double) ngen)/(10);
  
  for (g=0; g<ngen; ++g) {
    
    // **************************************************************************
    // ** 2.0 LOOP PRINT PROCESS ************************************************
    // **************************************************************************
    
    if(g == ngen - 1) {
      std::cout << "100.0%" << std::endl;
    } else if( g%ref_10==0) {
      std::cout << ((double)(g))/((double)(ngen))*100.0 << "%-" << std::endl;
    }
    
    // **************************************************************************
    // ** 2.1 EVALUATION EACH SET OF CHOICES ************************************
    // **************************************************************************
    
    for (int j=0; j<npop;++j) {
      Eigen::SparseVector<double> npm = currentVector;
      Eigen::SparseVector<double> lea(nbsites);
      cost = 0;
      
      for (int h=0;h<nbtoplant;++h){
        index = population(j,h);
        if(index != -1){
          x = viablesTriplets(index,0);
          y = viablesTriplets(index,1);

          cost = cost + costMatrix.coeffRef(x,y);
          lea = viablesValues.col(index);

          for (SparseVector<double>::InnerIterator it(lea); it; ++it)
          {
            it.value();
            it.index();
            v1 = npm.coeffRef(it.index()-1);
            v2 = it.value();
            npm.coeffRef(it.index()-1) = v1 + v2 - v1 * v2;
          }
        } 
      }
      
      Rcpp::NumericVector evaluate = eval_probabilityVector(npm, threshold+1);
      evaluation(j,1)=cost;
      
      if((evaluate(threshold)+evaluate(threshold+1))>=confidence){ 
        evaluation(j,0)=1;
        while(evaluate(threshold+1)>=confidence){
          int index_eliminated = optimise_remove_site(viablesValues,threshold,confidence,viablesTriplets,population,j,npm);
          if (index_eliminated==-1){
            break;
          }
          
          if (index_eliminated!=-1){
            evaluate = eval_probabilityVector(npm, threshold+1);
          }
        }
      }
      else{
        if(evaluate(threshold)>=confidence){ 
          optimise_add_site(viablesValues,threshold,confidence,viablesTriplets,population,j,npm);
          evaluation(j,0)=1;
        }
        else{
          evaluation(j,0)=0;
          evaluation(j,1)=10000000;
        }
      }
    }
    
    // **************************************************************************
    // ** 2.1 EVALUATION EACH SET OF CHOICES ************************************
    // **************************************************************************
    
    
    int saving;
    for (int h=0;h<npop;++h){
      ranks[h]=h;
    }
    for (int h=0;h<nsur;++h){
      indexmini = (h);
      cout_mini = evaluation(ranks(h),1);
      for (int k=h+1;k<npop;++k){
        if(cout_mini>evaluation(ranks(k),1)){
          indexmini = (k);
          cout_mini = evaluation(ranks(k),1);
        }
      }
      saving = ranks(h);
      ranks(h)=ranks(indexmini);
      ranks(indexmini)=saving;
    }
    
    for (int h=0;h<nsur;++h){
      survivors(h,_) = population(int(ranks(h)),_);
    }
    
    
    if (1==1 || g!=ngen){
      int mid = randomfunc3(nbtoplant);
      pheromons = 1 / (1+evaluation(_,1)-min(evaluation(_,1)));
      pheromons = pheromons / sum(pheromons);
      for (int h=nsur;h<npop;++h){
        papa = index_random_choice_non_uniform(pheromons);
        mama = index_random_choice_non_uniform(pheromons);
        for (int j=0;j<nbtoplant;++j){
          f = randomfunc3(2);
          if (f==1){
            survivors(h,j) = population(papa,j);
          }
          else{
            survivors(h,j) = population(mama,j);
          }
          
          f = randomfunc3(5);
          
          if (f==1){
            survivors(h,j) =  randomfunc3(nboc);
          }
        }
      }
      }
    if (true || g!=ngen){
      population = survivors;
      }
  
    }
    
  return(population);
}

//' rcpp_algorithm_opt2
//' 
//' Genetic algorithm used to optimised choices of introduction under the constraint of reaching a dertain number of 
//' presence at the end of the study period (threshold) with a minimal probability (confidence) while miimizing the cost
//' of all the introductions
//'
//' @param pheromons (Rcpp::NumericVector) weights of relative importance of each site&time pair.
//' @param viablesTriplets (Rcpp::NumericMatrix) For each viable site : index, position on map, time, and cost.
//' @param population0 (Rcpp::NumericMatrix) An initial population of set of site&time choices of introduction
//' @param costMatrix (Eigen::SparseMatrix<double>) matrix of costs of introduction, sparse.
//' @param currentPresenceMatrix (Eigen::SparseMatrix<double>) matrix of current presence of the species, sparse.
//' @param transitionmatrices (std::vector<Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.
//' @param globalSuitableSites (Eigen::SparseMatrix<double>) matrix of viable sites to use for the optimiation part.
//' @param viablesValues (Eigen::SparseMatrix<double>) viables values from transitions matrices from viable sites.
//' @param threshold (int) number of presence site to reach at the end of the studied period, with the introductions
//' @param confidence (double) minimal probability to reach the threshold accepted by the optimisation algorithm 
//' @param npop (int) genetic algorithm : size of the population
//' @param nsur (int) genetic algorithm : size of the surviving population
//' @param ngen (int) genetic algorithm : number of generations
//' @param nbtoplant evaluated maximum number of introduction necessary to satisfy the constraints of final presence.
//' @return an optimised population matrix of set of site&time choices of introduction.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_algorithm_opt2(Rcpp::NumericVector pheromons,
                                       Rcpp::NumericMatrix viablesTriplets,
                                       Rcpp::NumericMatrix population0,
                                       Eigen::SparseMatrix<double> costMatrix,
                                       Eigen::SparseMatrix<double> currentPresenceMatrix,
                                       std::vector<Eigen::SparseMatrix<double>> transitionmatrices,
                                       Eigen::SparseMatrix<double> globalSuitableSites,
                                       Eigen::SparseMatrix<double> viablesValues,
                                       int threshold,
                                       double confidence,
                                       int npop,
                                       int nsur,
                                       int ngen,
                                       int nbtoplant){
  
  // **************************************************************************
  // ** 0 INITIALISATION ******************************************************
  // **************************************************************************
  
  int index;    
  int g;
  int x;
  int y;
  int s;
  double v;
  double cost;
  int indexmini;
  double cout_mini;
  int papa;
  int mama;
  int f;
  double v1;
  double v2;
  int progress_bar_var;
  
  double best = 10000;
  
  Rcpp::NumericMatrix population=population0;
  Rcpp::NumericMatrix survivors(npop,nbtoplant);
  Rcpp::NumericVector ranks(npop);
  Rcpp::NumericVector pheromons2(npop);
  
  // **************************************************************************
  // ** 1 EVALUATION FINAL STATE **********************************************
  // **************************************************************************
  
  int nboc = viablesTriplets.nrow();
  int nbsites = (globalSuitableSites).nonZeros();
  Eigen::SparseVector<double> currentVector(nbsites);
  for (int k=0; k<currentPresenceMatrix.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(currentPresenceMatrix,k); it; ++it)
    {
      x = it.row();
      y = it.row();
      s = globalSuitableSites.coeffRef(it.row(),it.col());
      
      for (Eigen::SparseMatrix<double>::InnerIterator it2(transitionmatrices.back(),s-1); it2; ++it2) {
        v = it2.value();
        currentVector.coeffRef(it2.row()) = currentVector.coeffRef(it2.row()) +v-v*currentVector.coeffRef(it2.row());
      }
    }
    
    // **************************************************************************
    // ** 2 SELECTION LOOP PROCESS **********************************************
    // **************************************************************************
    
  Rcpp::NumericMatrix evaluation(npop,2);
  int ref_10 = (int) ((double) ngen)/(10) ;
  
  if(g%10!=0){
    std::cout << "Error, ngen must be multiple of 10. ngen setted to 10." << std::endl;
    ngen = 10;
  }
  
  for (g=0; g<ngen; ++g) {
    
    // **************************************************************************
    // ** 2.0 LOOP PRINT PROCESS ************************************************
    // **************************************************************************
    
    
    if(g == ngen - 1) {
      std::cout << "100.0%" << std::endl;
    } else if( g%ref_10==0) {
      std::cout << ((double)(g))/((double)(ngen))*100.0 << "%-" << std::endl;
    }
    
    // **************************************************************************
    // ** 2.1 EVALUATION EACH SET OF CHOICES ************************************
    // **************************************************************************
    
    //a loop over each individuals of the population of set of choices.
    for (int j=0; j<npop;++j) {
      //vector of probability of presence, will be modified by the choices of introduction.
      Eigen::SparseVector<double> npm = currentVector;
      //vector that will contain the introduction effect on each site.
      Eigen::SparseVector<double> lea(nbsites);
      //cost of the individual
      cost = 0;
      //evaluation=0 means the inidividual isn't validated; if it allow to reach the threshold of number of presence, the value change to 1.
      evaluation(j,0)=0;
      //then is the total cost of introductions, by default set to a very very high value.
      evaluation(j,1)=9999999999999;

      
      // Add one by one the introduction and the resulting probabability of presence on the final probability of presence vector.
      // if it's enought, the loop stop and any further introduction is cancelled and remove from the individual.
      
      // completed : check if a individual has reached the threshold, even without taking into account all of the introductions.
      bool completed = false;
      
    
      for (int h=0;h<nbtoplant;++h){
        
        if (completed){
          population(j,h)= -1;
        }
        
        // we obtain the index of the introduction
        index = population(j,h);
        
        if((index!=-1) & !completed){
          x = viablesTriplets(index,0);
          y = viablesTriplets(index,1);
          
          cost = cost + costMatrix.coeffRef(x,y);
          lea = viablesValues.col(index); // get the vector of resulting presence of the introduction indexed
          
          // then we add it to npm, .
          for (SparseVector<double>::InnerIterator it(lea); it; ++it)
          {
            it.value();
            it.index();
            v1 = npm.coeffRef(it.index()-1);
            v2 = it.value();
            npm.coeffRef(it.index()-1) = v1 + v2 - v1 * v2;
          }
          
          //We evaluate if, with these probabilities of presence at the end, we reach a sum of presence > threshold, with a certain confidence
          Rcpp::NumericVector evaluate = eval_probabilityVector(npm, threshold);
          if(evaluate(threshold)>=confidence){
            completed = true;
            evaluation(j,0)=1;
            evaluation(j,1)=cost;
          }
        } 
      }
    }
    
    // **************************************************************************
    // ** 2.2 RANKING EACH SET OF CHOICES ***************************************
    // **************************************************************************
    
    // ranking each individual
    int saving;
    for (int h=0;h<npop;++h){
      ranks[h]=h;
    }
    for (int h=0;h<nsur;++h){
      indexmini = (h);
      cout_mini = evaluation(ranks(h),1);
      for (int k=h+1;k<npop;++k){
        if(cout_mini>evaluation(ranks(k),1)){
          indexmini = (k);
          cout_mini = evaluation(ranks(k),1);
        }
      }
      saving = ranks(h);
      ranks(h)=ranks(indexmini);
      ranks(indexmini)=saving;
    }
    
    // building the first surviving individuals
    for (int h=0;h<nsur;++h){
      if(evaluation(int(ranks(h)),0)==1){
        survivors(h,_) = population(int(ranks(h)),_);
      }
      else{
        // in case a potential survivor don't reach the threshold, he's replaced by the best one.
        // this situation shouldn't happen in normal cases.
        survivors(h,_) = population(int(ranks(1)),_);
      }
        
    }
    
    // building childs individuals
    if (1==1 || g!=ngen){
      int mid = randomfunc3(nbtoplant);
      for (int h=nsur;h<npop;++h){
        // get parents from the survivors 
        papa = int(ranks(randomfunc3(nsur)));
        mama = int(ranks(randomfunc3(nsur)));
        for (int j=0;j<nbtoplant;++j){
          f = randomfunc3(2);//for each sites of an individual, randomly choose the site from the "papa" or "mama" individual
          if (f==1){
            survivors(h,j) = population(papa,j);
          }
          else{
            survivors(h,j) = population(mama,j);
          }
          f = randomfunc3(100);//there is a 1% chance of randomly change a site.  
          if (f==1){
            survivors(h,j) = index_random_choice_non_uniform(pheromons);
          }
        }
      }
    }
    if (true || g!=ngen){
      population = survivors;
    }
    
  
  }
  return(population);
}


//' Rewriting the final population given by the genetic algorithm with coordinates and times of introduction
//'
//' @param lastPopulation (Rcpp::NumericMatrix) a matrix of an optimised population obtain with the genetic algorithm
//' @param viablesTriplets (Rcpp::NumericMatrix) a matrix of information about each viable site.
//' @return the choices of introduction in coordinates and timesteps. 1st&2nd columns are the XY coordinates, 3rd is the timesteps, 4st is the cost, 5st is the index on the grid.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_result_to_choice(Rcpp::NumericMatrix lastPopulation,
                                          Rcpp::NumericMatrix viablesTriplets){
  int nbtoplant =(lastPopulation).ncol();
  Rcpp::NumericMatrix  choices_matrix(nbtoplant,6);  
  if(true){
    for (int j=0;j<nbtoplant;++j){
      //survivors(h,j) = population(int(ranks(h%nsur)),j);
      if( lastPopulation(1,j)!=-1){
        choices_matrix(j,0) = viablesTriplets(lastPopulation(1,j),0)+1;
        choices_matrix(j,1) = viablesTriplets(lastPopulation(1,j),1)+1;
        choices_matrix(j,2) = viablesTriplets(lastPopulation(1,j),2)+1;
        choices_matrix(j,3) = viablesTriplets(lastPopulation(1,j),3);
        choices_matrix(j,4) = viablesTriplets(lastPopulation(1,j),4)+1;
        //lastPopulation(j,5) = (viablesValues2.col(lastPopulation(1,j))).sum();
      }
      else {
        lastPopulation(j,0) = -1;
        lastPopulation(j,1) = -1;
        lastPopulation(j,2) = -1;
        lastPopulation(j,3) = -1;
      }
    }
  }
  return (choices_matrix);
}

