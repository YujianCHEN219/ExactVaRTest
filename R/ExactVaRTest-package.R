#' ExactVaRTest – Exact Finite-Sample VaR Back-Testing
#'
#' Provides fast dynamic-programming algorithms (C++/Rcpp) – with pure-R
#' fall-backs – for the exact finite-sample distributions of Christoffersen’s
#' (1998) Unconditional Coverage (UC), Independence (IND) and Conditional
#' Coverage (CC) tests.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib ExactVaRTest, .registration = TRUE
"_PACKAGE"
