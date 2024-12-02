% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_viable_sites}
\alias{rcpp_viable_sites}
\title{rcpp_viable_sites}
\usage{
rcpp_viable_sites(transitionmatrices)
}
\arguments{
\item{transitionmatrices}{(Eigen::SparseMatrix<double>>) list of transition matrices from any timestep to the last timestep.}
}
\value{
boolean sparse matrix, each row correspond to a time, each column to a site. (1=viable; 0=not viable)
}
\description{
Calculate sites&times introduction considered "viable" 
eg that are not totally bested by other timings of introduction for the same site, 
and that are efficient enough. (colonises more than one site on average).
}
