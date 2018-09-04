// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _FASTR_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// YPR
NumericMatrix YPR(NumericVector age, NumericVector sel, NumericVector m, NumericVector mat, NumericVector wgt, double alpha, double beta, double sigma, double inc);
RcppExport SEXP _FASTR_YPR(SEXP ageSEXP, SEXP selSEXP, SEXP mSEXP, SEXP matSEXP, SEXP wgtSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP sigmaSEXP, SEXP incSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type age(ageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sel(selSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mat(matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wgt(wgtSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type inc(incSEXP);
    rcpp_result_gen = Rcpp::wrap(YPR(age, sel, m, mat, wgt, alpha, beta, sigma, inc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FASTR_rcpp_hello_world", (DL_FUNC) &_FASTR_rcpp_hello_world, 0},
    {"_FASTR_YPR", (DL_FUNC) &_FASTR_YPR, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_FASTR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
