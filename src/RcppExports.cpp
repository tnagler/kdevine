// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kern_gauss
NumericVector kern_gauss(const NumericVector x);
RcppExport SEXP _kdevine_kern_gauss(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(kern_gauss(x));
    return rcpp_result_gen;
END_RCPP
}
// ikern_gauss
NumericVector ikern_gauss(const NumericVector x);
RcppExport SEXP _kdevine_ikern_gauss(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ikern_gauss(x));
    return rcpp_result_gen;
END_RCPP
}
// eval_kde1d
NumericVector eval_kde1d(const NumericVector xsort, const NumericVector xev, const double xmin, const double xmax, const double bw);
RcppExport SEXP _kdevine_eval_kde1d(SEXP xsortSEXP, SEXP xevSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type xsort(xsortSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xev(xevSEXP);
    Rcpp::traits::input_parameter< const double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const double >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_kde1d(xsort, xev, xmin, xmax, bw));
    return rcpp_result_gen;
END_RCPP
}
// eval_pkde1d
NumericVector eval_pkde1d(const NumericVector x, const NumericVector xev, const double xmin, const double xmax, const double bw);
RcppExport SEXP _kdevine_eval_pkde1d(SEXP xSEXP, SEXP xevSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type xev(xevSEXP);
    Rcpp::traits::input_parameter< const double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const double >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_pkde1d(x, xev, xmin, xmax, bw));
    return rcpp_result_gen;
END_RCPP
}
// eval_qkde1d
NumericVector eval_qkde1d(const NumericVector x, const NumericVector qev, const double xmin, const double xmax, const double bw);
RcppExport SEXP _kdevine_eval_qkde1d(SEXP xSEXP, SEXP qevSEXP, SEXP xminSEXP, SEXP xmaxSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type qev(qevSEXP);
    Rcpp::traits::input_parameter< const double >::type xmin(xminSEXP);
    Rcpp::traits::input_parameter< const double >::type xmax(xmaxSEXP);
    Rcpp::traits::input_parameter< const double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_qkde1d(x, qev, xmin, xmax, bw));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kdevine_kern_gauss", (DL_FUNC) &_kdevine_kern_gauss, 1},
    {"_kdevine_ikern_gauss", (DL_FUNC) &_kdevine_ikern_gauss, 1},
    {"_kdevine_eval_kde1d", (DL_FUNC) &_kdevine_eval_kde1d, 5},
    {"_kdevine_eval_pkde1d", (DL_FUNC) &_kdevine_eval_pkde1d, 5},
    {"_kdevine_eval_qkde1d", (DL_FUNC) &_kdevine_eval_qkde1d, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_kdevine(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
