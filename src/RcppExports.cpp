// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppSMC.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcppsmc_logNormWeightsCpp
Rcpp::List rcppsmc_logNormWeightsCpp(arma::vec logWeights);
RcppExport SEXP _RcppSMCkalman_rcppsmc_logNormWeightsCpp(SEXP logWeightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type logWeights(logWeightsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppsmc_logNormWeightsCpp(logWeights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppSMCkalman_rcppsmc_logNormWeightsCpp", (DL_FUNC) &_RcppSMCkalman_rcppsmc_logNormWeightsCpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppSMCkalman(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}