#ifndef LOPC
#define LOPC

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

int half_random_remove_site(Eigen::SparseMatrix<double> viablesValues2, 
                    int threshold, double confidence, 
                    Rcpp::NumericMatrix viablesTriplets, 
                    Rcpp::NumericMatrix population, 
                    int pop,
                    Eigen::SparseVector<double> current);

int optimise_remove_site(Eigen::SparseMatrix<double> viablesValues2, 
                         int threshold, 
                         double confidence, 
                         Rcpp::NumericMatrix viablesTriplets, 
                         Rcpp::NumericMatrix population, 
                         int pop,
                         Eigen::SparseVector<double> current);

void optimise_add_site(Eigen::SparseMatrix<double> viablesValues2, 
                       int threshold, double confidence, 
                       Rcpp::NumericMatrix viablesTriplets, 
                       Rcpp::NumericMatrix population,
                       int pop,
                       Eigen::SparseVector<double> current);

int optimise_remove_site2(Eigen::SparseMatrix<double> viablesValues2, 
                         int threshold, 
                         double confidence, 
                         Rcpp::NumericMatrix viablesTriplets, 
                         Rcpp::NumericMatrix population, 
                         int pop,
                         Eigen::SparseVector<double> current);

void optimise_add_site2(Eigen::SparseMatrix<double> viablesValues2, 
                       int threshold, double confidence, 
                       Rcpp::NumericMatrix viablesTriplets, 
                       Rcpp::NumericMatrix population,
                       int pop,
                       Eigen::SparseVector<double> current);

#endif