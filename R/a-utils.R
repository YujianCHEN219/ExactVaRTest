# --------------------------------------------------------------------
#  Internal utility helpers  ——  NOT exported
# --------------------------------------------------------------------

#' Machine-safe lower bound used across the package
#'
#' A single place to store the numeric constant `1e-15`, so we do not
#' rewrite it in every algorithm.  **Internal use only.**
#'
#' @keywords internal
#' @noRd
EPS <- 1e-15

#' Numerically safe logarithm
#'
#' Replaces values smaller than `EPS` with `EPS` before taking `log`,
#' preventing `-Inf` when probabilities underflow.  
#' **Internal use only, therefore _no export tag_.**
#'
#' @param x Numeric vector.
#' @return `log(pmax(x, EPS))`
#' @keywords internal
#' @noRd
safe_log <- function(x) {
  log(pmax(x, EPS))
}
