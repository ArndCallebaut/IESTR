// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_set_seed
void rcpp_set_seed(unsigned int seed);
RcppExport SEXP _IESTR_rcpp_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_set_seed(seed);
    return R_NilValue;
END_RCPP
}
// index_random_choice_non_uniform
int index_random_choice_non_uniform(Rcpp::NumericVector ununiform_probabilities0);
RcppExport SEXP _IESTR_index_random_choice_non_uniform(SEXP ununiform_probabilities0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ununiform_probabilities0(ununiform_probabilities0SEXP);
    rcpp_result_gen = Rcpp::wrap(index_random_choice_non_uniform(ununiform_probabilities0));
    return rcpp_result_gen;
END_RCPP
}
// generate_permutation4
NumericVector generate_permutation4(int permutation_size, Rcpp::NumericVector ununiform_probabilities0);
RcppExport SEXP _IESTR_generate_permutation4(SEXP permutation_sizeSEXP, SEXP ununiform_probabilities0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type permutation_size(permutation_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ununiform_probabilities0(ununiform_probabilities0SEXP);
    rcpp_result_gen = Rcpp::wrap(generate_permutation4(permutation_size, ununiform_probabilities0));
    return rcpp_result_gen;
END_RCPP
}
// proba_matrix_mult3
Eigen::SparseMatrix<double> proba_matrix_mult3(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B);
RcppExport SEXP _IESTR_proba_matrix_mult3(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(proba_matrix_mult3(A, B));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_global_suitable_sites
Eigen::SparseMatrix<double> rcpp_global_suitable_sites(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix);
RcppExport SEXP _IESTR_rcpp_global_suitable_sites(SEXP consecutiveSuitabilityMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::list<Eigen::SparseMatrix<double>> >::type consecutiveSuitabilityMatrix(consecutiveSuitabilityMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_global_suitable_sites(consecutiveSuitabilityMatrix));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_global_suitable_coordinates
Rcpp::NumericMatrix rcpp_global_suitable_coordinates(Eigen::SparseMatrix<double> globalSuitableSites);
RcppExport SEXP _IESTR_rcpp_global_suitable_coordinates(SEXP globalSuitableSitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_global_suitable_coordinates(globalSuitableSites));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_spread_matrix
Eigen::SparseMatrix<double> rcpp_spread_matrix(Eigen::SparseMatrix<double> globalSuitableSites, Rcpp::NumericMatrix globalSuitableCoordinates, Rcpp::NumericMatrix migrationKernel);
RcppExport SEXP _IESTR_rcpp_spread_matrix(SEXP globalSuitableSitesSEXP, SEXP globalSuitableCoordinatesSEXP, SEXP migrationKernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type globalSuitableCoordinates(globalSuitableCoordinatesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type migrationKernel(migrationKernelSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_spread_matrix(globalSuitableSites, globalSuitableCoordinates, migrationKernel));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_local_transition_matrix
std::list<Eigen::SparseMatrix<double>> rcpp_local_transition_matrix(std::list<Eigen::SparseMatrix<double>> consecutiveSuitabilityMatrix, Eigen::SparseMatrix<double> localTransitionMatrix, Rcpp::NumericMatrix globalSuitableCoordinates);
RcppExport SEXP _IESTR_rcpp_local_transition_matrix(SEXP consecutiveSuitabilityMatrixSEXP, SEXP localTransitionMatrixSEXP, SEXP globalSuitableCoordinatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::list<Eigen::SparseMatrix<double>> >::type consecutiveSuitabilityMatrix(consecutiveSuitabilityMatrixSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type localTransitionMatrix(localTransitionMatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type globalSuitableCoordinates(globalSuitableCoordinatesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_local_transition_matrix(consecutiveSuitabilityMatrix, localTransitionMatrix, globalSuitableCoordinates));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_transition_matrix
std::vector<Eigen::SparseMatrix<double>> rcpp_transition_matrix(std::list<Eigen::SparseMatrix<double>> transitionMatrices);
RcppExport SEXP _IESTR_rcpp_transition_matrix(SEXP transitionMatricesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::list<Eigen::SparseMatrix<double>> >::type transitionMatrices(transitionMatricesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_transition_matrix(transitionMatrices));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_viable_sites
Eigen::SparseMatrix<double> rcpp_viable_sites(std::vector<Eigen::SparseMatrix<double>> colonisationMatrices);
RcppExport SEXP _IESTR_rcpp_viable_sites(SEXP colonisationMatricesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<Eigen::SparseMatrix<double>> >::type colonisationMatrices(colonisationMatricesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_viable_sites(colonisationMatrices));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_viable_triplets
Rcpp::NumericMatrix rcpp_viable_triplets(Eigen::SparseMatrix<double> viableSites, std::list<Eigen::SparseMatrix<double>> colonisationMatrices, Rcpp::NumericMatrix globalSuitableCoordinates, Eigen::SparseMatrix<double> globalSuitableSites, Eigen::SparseMatrix<double> costMatrix);
RcppExport SEXP _IESTR_rcpp_viable_triplets(SEXP viableSitesSEXP, SEXP colonisationMatricesSEXP, SEXP globalSuitableCoordinatesSEXP, SEXP globalSuitableSitesSEXP, SEXP costMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type viableSites(viableSitesSEXP);
    Rcpp::traits::input_parameter< std::list<Eigen::SparseMatrix<double>> >::type colonisationMatrices(colonisationMatricesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type globalSuitableCoordinates(globalSuitableCoordinatesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type costMatrix(costMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_viable_triplets(viableSites, colonisationMatrices, globalSuitableCoordinates, globalSuitableSites, costMatrix));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_viable_values
Eigen::SparseMatrix<double> rcpp_viable_values(Rcpp::NumericMatrix viablesTriplets, Eigen::SparseMatrix<double> viableSites, Eigen::SparseMatrix<double> globalSuitableSites, std::vector<Eigen::SparseMatrix<double>> colonisationMatrices);
RcppExport SEXP _IESTR_rcpp_viable_values(SEXP viablesTripletsSEXP, SEXP viableSitesSEXP, SEXP globalSuitableSitesSEXP, SEXP colonisationMatricesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type viablesTriplets(viablesTripletsSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type viableSites(viableSitesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    Rcpp::traits::input_parameter< std::vector<Eigen::SparseMatrix<double>> >::type colonisationMatrices(colonisationMatricesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_viable_values(viablesTriplets, viableSites, globalSuitableSites, colonisationMatrices));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_eval_current_prob
NumericVector rcpp_eval_current_prob(int threshold, Eigen::SparseMatrix<double> currentPresenceMatrix, std::vector<Eigen::SparseMatrix<double>> colonisationMatrices, Eigen::SparseMatrix<double> globalSuitableSites);
RcppExport SEXP _IESTR_rcpp_eval_current_prob(SEXP thresholdSEXP, SEXP currentPresenceMatrixSEXP, SEXP colonisationMatricesSEXP, SEXP globalSuitableSitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type currentPresenceMatrix(currentPresenceMatrixSEXP);
    Rcpp::traits::input_parameter< std::vector<Eigen::SparseMatrix<double>> >::type colonisationMatrices(colonisationMatricesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_eval_current_prob(threshold, currentPresenceMatrix, colonisationMatrices, globalSuitableSites));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_pheromons
Rcpp::NumericVector rcpp_pheromons(Rcpp::NumericMatrix viablesTriplets);
RcppExport SEXP _IESTR_rcpp_pheromons(SEXP viablesTripletsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type viablesTriplets(viablesTripletsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_pheromons(viablesTriplets));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_generate_population
Rcpp::NumericVector rcpp_generate_population(Rcpp::NumericVector pheromons, Eigen::SparseMatrix<double> globalSuitableSites, int npop, int nbtoplant);
RcppExport SEXP _IESTR_rcpp_generate_population(SEXP pheromonsSEXP, SEXP globalSuitableSitesSEXP, SEXP npopSEXP, SEXP nbtoplantSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pheromons(pheromonsSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type nbtoplant(nbtoplantSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_generate_population(pheromons, globalSuitableSites, npop, nbtoplant));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_algorithm_opt
Rcpp::NumericVector rcpp_algorithm_opt(Rcpp::NumericVector pheromons, Rcpp::NumericMatrix viablesTriplets, Rcpp::NumericMatrix population0, Eigen::SparseMatrix<double> costMatrix, Eigen::SparseMatrix<double> currentPresenceMatrix, std::vector<Eigen::SparseMatrix<double>> colonisationMatrices, Eigen::SparseMatrix<double> globalSuitableSites, Eigen::SparseMatrix<double> viablesValues, int threshold, double confidence, int npop, int nsur, int ngen, int nbtoplant);
RcppExport SEXP _IESTR_rcpp_algorithm_opt(SEXP pheromonsSEXP, SEXP viablesTripletsSEXP, SEXP population0SEXP, SEXP costMatrixSEXP, SEXP currentPresenceMatrixSEXP, SEXP colonisationMatricesSEXP, SEXP globalSuitableSitesSEXP, SEXP viablesValuesSEXP, SEXP thresholdSEXP, SEXP confidenceSEXP, SEXP npopSEXP, SEXP nsurSEXP, SEXP ngenSEXP, SEXP nbtoplantSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pheromons(pheromonsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type viablesTriplets(viablesTripletsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type population0(population0SEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type currentPresenceMatrix(currentPresenceMatrixSEXP);
    Rcpp::traits::input_parameter< std::vector<Eigen::SparseMatrix<double>> >::type colonisationMatrices(colonisationMatricesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type globalSuitableSites(globalSuitableSitesSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type viablesValues(viablesValuesSEXP);
    Rcpp::traits::input_parameter< int >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type confidence(confidenceSEXP);
    Rcpp::traits::input_parameter< int >::type npop(npopSEXP);
    Rcpp::traits::input_parameter< int >::type nsur(nsurSEXP);
    Rcpp::traits::input_parameter< int >::type ngen(ngenSEXP);
    Rcpp::traits::input_parameter< int >::type nbtoplant(nbtoplantSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_algorithm_opt(pheromons, viablesTriplets, population0, costMatrix, currentPresenceMatrix, colonisationMatrices, globalSuitableSites, viablesValues, threshold, confidence, npop, nsur, ngen, nbtoplant));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_result_to_choice
Rcpp::NumericMatrix rcpp_result_to_choice(Rcpp::NumericMatrix lastPopulation, Rcpp::NumericMatrix viablesTriplets);
RcppExport SEXP _IESTR_rcpp_result_to_choice(SEXP lastPopulationSEXP, SEXP viablesTripletsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type lastPopulation(lastPopulationSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type viablesTriplets(viablesTripletsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_result_to_choice(lastPopulation, viablesTriplets));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IESTR_rcpp_set_seed", (DL_FUNC) &_IESTR_rcpp_set_seed, 1},
    {"_IESTR_index_random_choice_non_uniform", (DL_FUNC) &_IESTR_index_random_choice_non_uniform, 1},
    {"_IESTR_generate_permutation4", (DL_FUNC) &_IESTR_generate_permutation4, 2},
    {"_IESTR_proba_matrix_mult3", (DL_FUNC) &_IESTR_proba_matrix_mult3, 2},
    {"_IESTR_rcpp_global_suitable_sites", (DL_FUNC) &_IESTR_rcpp_global_suitable_sites, 1},
    {"_IESTR_rcpp_global_suitable_coordinates", (DL_FUNC) &_IESTR_rcpp_global_suitable_coordinates, 1},
    {"_IESTR_rcpp_spread_matrix", (DL_FUNC) &_IESTR_rcpp_spread_matrix, 3},
    {"_IESTR_rcpp_local_transition_matrix", (DL_FUNC) &_IESTR_rcpp_local_transition_matrix, 3},
    {"_IESTR_rcpp_transition_matrix", (DL_FUNC) &_IESTR_rcpp_transition_matrix, 1},
    {"_IESTR_rcpp_viable_sites", (DL_FUNC) &_IESTR_rcpp_viable_sites, 1},
    {"_IESTR_rcpp_viable_triplets", (DL_FUNC) &_IESTR_rcpp_viable_triplets, 5},
    {"_IESTR_rcpp_viable_values", (DL_FUNC) &_IESTR_rcpp_viable_values, 4},
    {"_IESTR_rcpp_eval_current_prob", (DL_FUNC) &_IESTR_rcpp_eval_current_prob, 4},
    {"_IESTR_rcpp_pheromons", (DL_FUNC) &_IESTR_rcpp_pheromons, 1},
    {"_IESTR_rcpp_generate_population", (DL_FUNC) &_IESTR_rcpp_generate_population, 4},
    {"_IESTR_rcpp_algorithm_opt", (DL_FUNC) &_IESTR_rcpp_algorithm_opt, 14},
    {"_IESTR_rcpp_result_to_choice", (DL_FUNC) &_IESTR_rcpp_result_to_choice, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_IESTR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
