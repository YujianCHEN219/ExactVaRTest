#' ExactVaRTest – Exact Finite-Sample VaR Back-Testing
#'
#' Provides fast dynamic‑programming algorithms (C++/Rcpp) – with pure‑R
#' fall‑backs – for the exact finite‑sample distributions and p‑values of
#' Christoffersen’s (1998) VaR back‑tests: Independence (IND) and Conditional
#' Coverage (CC) tests, and the Unconditional Coverage (UC) test via closed‑form
#' binomial enumeration.

#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib ExactVaRTest, .registration = TRUE
"_PACKAGE"
