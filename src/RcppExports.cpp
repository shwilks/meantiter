// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calc_mean_titer_negll
double calc_mean_titer_negll(const double& predicted_mean, const arma::vec& max_titers, const arma::vec& min_titers, const double& titer_sd);
RcppExport SEXP _meantiter_calc_mean_titer_negll(SEXP predicted_meanSEXP, SEXP max_titersSEXP, SEXP min_titersSEXP, SEXP titer_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type predicted_mean(predicted_meanSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type max_titers(max_titersSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type min_titers(min_titersSEXP);
    Rcpp::traits::input_parameter< const double& >::type titer_sd(titer_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_mean_titer_negll(predicted_mean, max_titers, min_titers, titer_sd));
    return rcpp_result_gen;
END_RCPP
}
// calc_mean_titer_negll_by_par
double calc_mean_titer_negll_by_par(const arma::vec& pars, const arma::vec& max_titers, const arma::vec& min_titers, double titer_sd);
RcppExport SEXP _meantiter_calc_mean_titer_negll_by_par(SEXP parsSEXP, SEXP max_titersSEXP, SEXP min_titersSEXP, SEXP titer_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type max_titers(max_titersSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type min_titers(min_titersSEXP);
    Rcpp::traits::input_parameter< double >::type titer_sd(titer_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_mean_titer_negll_by_par(pars, max_titers, min_titers, titer_sd));
    return rcpp_result_gen;
END_RCPP
}
// calc_mean_titer_ci_by_par
double calc_mean_titer_ci_by_par(const arma::vec& pars, const arma::vec& max_titers, const arma::vec& min_titers, double titer_sd, double target_negll);
RcppExport SEXP _meantiter_calc_mean_titer_ci_by_par(SEXP parsSEXP, SEXP max_titersSEXP, SEXP min_titersSEXP, SEXP titer_sdSEXP, SEXP target_negllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type max_titers(max_titersSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type min_titers(min_titersSEXP);
    Rcpp::traits::input_parameter< double >::type titer_sd(titer_sdSEXP);
    Rcpp::traits::input_parameter< double >::type target_negll(target_negllSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_mean_titer_ci_by_par(pars, max_titers, min_titers, titer_sd, target_negll));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_meantiter_calc_mean_titer_negll", (DL_FUNC) &_meantiter_calc_mean_titer_negll, 4},
    {"_meantiter_calc_mean_titer_negll_by_par", (DL_FUNC) &_meantiter_calc_mean_titer_negll_by_par, 4},
    {"_meantiter_calc_mean_titer_ci_by_par", (DL_FUNC) &_meantiter_calc_mean_titer_ci_by_par, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_meantiter(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}