


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
#include <random>
#include <chrono>
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Rcpp;
using namespace Eigen;

#include <Rcpp.h>

//' Set the random seed for C++ 
//'
//' @param seed seed value
// [[Rcpp::export]]
void rcpp_set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

//' return random value from ununiform probability set. 
//'
//' @param ununiform_probabilities0 Rcpp vector of probabilities for each indices value
//[[Rcpp::export]]
int index_random_choice_non_uniform (Rcpp::NumericVector ununiform_probabilities0){
  Rcpp::NumericVector ununiform_probabilities = ununiform_probabilities0/Rcpp::sum(ununiform_probabilities0);
  Rcpp::NumericVector cumulative_probabilities(ununiform_probabilities.length());
  cumulative_probabilities(0) = ununiform_probabilities(0);
  for (int h=1; h<ununiform_probabilities.length() ; h++){
    cumulative_probabilities(h) = ununiform_probabilities(h)+cumulative_probabilities(h-1);  
  }
  double randy = ((double) rand() / (RAND_MAX)) ;
  for (int h = 0 ; h < ununiform_probabilities.length() ; h++){
    if (randy<cumulative_probabilities(h)){
      return(int(h));
    }
  }
  return(ununiform_probabilities.length()-1);
}

//' return the probability that the number of success of Bernouillis event with probabilities given is between
//' zero and the given threshold 
//'
//' @param probabilityVector vector of probabilities for each bernouilli event considered.
//' @param threshold1 maximum value of success tested, inclued cases with more success.
NumericVector eval_probabilityVector (Eigen::SparseVector<double> probabilityVector, int threshold1){
  
  int const threshold = threshold1;
  Rcpp::NumericVector nScores(threshold+1);
  int h;
  for (h=0; h<=threshold; ++h){
    nScores(h) = 0;
  }
  nScores(0)=1;
  for (Eigen::SparseVector<double>::InnerIterator it(probabilityVector); it; ++it)
  {
    nScores(threshold) =  it.value() * nScores(threshold-1) +  nScores(threshold);
    for (h=(threshold-1); h>=1; --h){
      nScores(h) =  it.value() * nScores(h-1) +  (1-it.value()) * nScores(h);
    }
    nScores[0] = ((1-it.value()) * nScores[0]);
  }
  return(nScores) ;
}

//' return the evaluation of number of success like eval_probabilitVector, with two vector summed.
//'
//' @param probabilityVector first vector
//' @param probabilityVector2 second vector
//' @param threshold1 maximum value of success tested, inclued cases with more success.
NumericVector eval_probabilityVector_adding (Eigen::SparseVector<double> probabilityVector, Eigen::SparseVector<double> probabilityVector2,int threshold1){
  
  Eigen::SparseVector<double> cury = probabilityVector;
  for (InIterVec i_(probabilityVector2); i_; ++i_){
    //Rcpp::Rcout << " i=" << i_.index() << " value=" << i_.value() << std::endl;
    cury.coeffRef(i_.index()) = cury.coeffRef(i_.index()) + i_.value() - cury.coeffRef(i_.index()) * i_.value()  ;
  }
  return(eval_probabilityVector(cury,threshold1)) ;
}

//' get a random value between 0 & (j-1)
//'
//' @param j limit-1 of random value
int randomfunc3(int j){
  return rand() % j;
}

//' generate a permutation.
//'
//' @param  vector of probabilities for each bernouilli event considered.
//' @param threshold1 maximum value of success tested, inclued cases with more success.
NumericVector generate_permutation3(int permutation_size, int total_size){
  
  Rcpp::NumericVector v(total_size);
  for (int i = 0; i < total_size; i++){
    v(i) = i ;
  }
  
  Rcpp::NumericVector result(permutation_size);
  int indice_max = total_size;
  int alea;
  int saver;
  for (int i = 0; i < permutation_size; i++){
    alea = randomfunc3(indice_max);
    //std::cout << "OK here it's done..."<<v(alea) << std::endl;
    result(i) = v(alea);
    saver = v(indice_max-1);
    v(indice_max-1) = v(alea);
    v(alea) = saver;
    indice_max --;
  }
  return(result);
}

//' return a permutation of numbers with ununiform probabilities of being chosen.
//'
//' @param permutation_size size of the permutation
//' @param ununiform_probabilities0 relative probabilities of being chosen
//[[Rcpp::export]]
NumericVector generate_permutation4(int permutation_size, Rcpp::NumericVector ununiform_probabilities0){

  Rcpp::NumericVector ununiform_probabilities = ununiform_probabilities0;
  int total_size = ununiform_probabilities.length();
  Rcpp::NumericVector v(total_size);
  
  for (int i = 0; i < total_size; i++){
    v(i) = i ;
  }

  Rcpp::NumericVector result(permutation_size);
  int indice_max = total_size;
  int alea;
  int saver;
  
  for (int i = 0; i < permutation_size; i++){

    alea = index_random_choice_non_uniform(ununiform_probabilities);
    result(i) = v(alea);
    saver = v(indice_max-1);
    v(indice_max-1) = v(alea);
    v(alea) = saver;
    
    ununiform_probabilities(v(alea)) = ununiform_probabilities(indice_max-1);
    ununiform_probabilities(indice_max-1) = 0;
    indice_max --;
    //ununiform_probabilities.erase((indice_max-1-i));
  }
  return(result);
}

//' return the special matrix product between two probability matrix. (obsolete version)
//'
//' @param A first probability matrix
//' @param B second probability matrix
Eigen::SparseMatrix<double> proba_matrix_mult(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  // obsolete version 
  int n = A.rows();
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      val1 = it.value();
      ind1 = it.row();
      if(val1>0.01){
        for (int h=0; h<B.outerSize(); ++h)
          for (Eigen::SparseMatrix<double>::InnerIterator it2(B,h); it2; ++it2)
          {
            if (it.col()==it2.row()){
              val2 = it2.value()*val1;
              ind2 = it2.col();
              if (val2>0.001){
                result_matrix.coeffRef(ind1,ind2) = result_matrix.coeffRef(ind1,ind2) + val2 - val2* result_matrix.coeffRef(ind1,ind2);
              }
            }
          }
      }
    }
    return (result_matrix);
}

//' (obsolete version) return the special matrix product between two probability matrix.
//'
//' @param A first probability matrix
//' @param B second probability matrix
Eigen::SparseMatrix<double> proba_matrix_mult2(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  // probabilistic matrix multiplication of A*B
  int n = A.rows();
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind;
  double val2;
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      val1 = it.value();
      ind = it.col();
      for (int h=0; h<n;++h){
        val2 = B.coeffRef(ind,h);
        if(val2!=0){
          result_matrix.coeffRef(it.row(),h)+= val1 + val2 - val1*val2;
        }
      }
    }
    return (result_matrix);
}

//' return the special matrix product between two probability matrix.
//'
//' @param A first probability matrix
//' @param B second probability matrix
//[[Rcpp::export]]
Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> A,Eigen::SparseMatrix<double> B){
  int n = A.rows();
  Eigen::SparseVector<double> local2(n);
  Eigen::SparseMatrix<double>result_matrix(n,n);
  double val1;
  int ind1;
  int ind2;
  double val2;
  double val3;
  int i;
  std::cout << "A..."<<A.nonZeros()<<" ------- B..."<<B.nonZeros()<<std::endl;
  double prec = 0.0005;
  
  for (int k=0; k<B.outerSize(); ++k){
    Eigen::SparseVector<double> local(n);
        for (Eigen::SparseMatrix<double>::InnerIterator it(B,k); it; ++it)
          {
            val1 = it.value();
            ind1 = it.row();
            if(val1>prec){
              local2 = val1 * A.col(ind1);
              for (Eigen::SparseVector<double>::InnerIterator it(local2); it; ++it)
                {
                local.coeffRef(it.index()) = local.coeffRef(it.index()) +it.value() - local.coeffRef(it.index()) *it.value() ;
                }
              }
            }
        for (InIterVec i_(local); i_; ++i_){
          if (i_.value()>prec){
            result_matrix.insert(i_.index(),k) = i_.value();
          }
        }
    }
  result_matrix.pruned(prec);
  result_matrix.makeCompressed();
  return (result_matrix);
}




