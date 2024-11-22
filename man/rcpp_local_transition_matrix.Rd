% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_local_transition_matrix}
\alias{rcpp_local_transition_matrix}
\title{rcpp_local_transition_matrix}
\usage{
rcpp_local_transition_matrix(
  consecutiveSuitabilityMatrix,
  spreadmatrix,
  globalSuitableCoordinates
)
}
\arguments{
\item{consecutiveSuitabilityMatrix}{(std::list<Eigen::SparseMatrix<double>>) list of suitability matrix}

\item{spreadmatrix}{(Eigen::SparseMatrix<double>) matrix used to model spread behavior of the species without suitability}

\item{globalSuitableCoordinates}{(Rcpp::NumericMatrix) support matrix coordinates of each sites}
}
\value{
list of matrix used to model spread behavior of the species with suitability, between two adjacent timestep
}
\description{
Calculate the local transition matrices from the spread matrix and suitability matrices
}
