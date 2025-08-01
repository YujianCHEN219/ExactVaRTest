// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fb_lrind_fastcpp
List fb_lrind_fastcpp(int n, double alpha, double prune_threshold);
RcppExport SEXP _ExactVaRTest_fb_lrind_fastcpp(SEXP nSEXP, SEXP alphaSEXP, SEXP prune_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type prune_threshold(prune_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(fb_lrind_fastcpp(n, alpha, prune_threshold));
    return rcpp_result_gen;
END_RCPP
}
// fb_lrcc_fastcpp
List fb_lrcc_fastcpp(int n, double alpha, double prune_threshold);
RcppExport SEXP _ExactVaRTest_fb_lrcc_fastcpp(SEXP nSEXP, SEXP alphaSEXP, SEXP prune_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type prune_threshold(prune_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(fb_lrcc_fastcpp(n, alpha, prune_threshold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ExactVaRTest_fb_lrind_fastcpp", (DL_FUNC) &_ExactVaRTest_fb_lrind_fastcpp, 3},
    {"_ExactVaRTest_fb_lrcc_fastcpp", (DL_FUNC) &_ExactVaRTest_fb_lrcc_fastcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ExactVaRTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
