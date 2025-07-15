#  -----------------------------------------------------------------------
#  High-level wrappers that auto-select the fastest available engine
#  -----------------------------------------------------------------------

#' Exact LR_ind distribution (auto-select engine)
#'
#' Returns the exact finite-sample distribution of Christoffersen’s
#' independence statistic \eqn{LR_{\mathrm{ind}}}.
#' The function first tries the Rcpp implementation
#' (`fb_lrind_fastcpp`) and silently falls back to the pure-R reference
#' version (`fb_lrind_R`) if the shared library is not available.
#'
#' @param n Integer. Sample size (\eqn{n \ge 1}{n >= 1}).
#' @param alpha Numeric in (0,1). Exception probability.
#' @param prune_threshold Numeric. States whose probability drops below
#'   this value are pruned in the dynamic-programming recursion.
#'
#' @return A data frame with two columns
#'   \describe{
#'     \item{LR}{unique LR values (ascending)}
#'     \item{prob}{their exact probabilities (sum to 1)}
#'   }
#' @export
lr_ind_dist <- function(n, alpha = 0.05, prune_threshold = 1e-15) {
  if (exists("fb_lrind_fastcpp", mode = "function")) {
    tryCatch(
      fb_lrind_fastcpp(n, alpha, prune_threshold),
      error = function(e) fb_lrind_R(n, alpha, prune_threshold)
    )
  } else {
    fb_lrind_R(n, alpha, prune_threshold)
  }
}

#' Exact LR_cc distribution (auto-select engine)
#'
#' Returns the exact finite-sample distribution of Christoffersen’s
#' conditional-coverage statistic \eqn{LR_{\mathrm{cc}}}.
#' The fast Rcpp backend (`fb_lrcc_fastcpp`) is used when available; if
#' compilation fails the function falls back to the pure-R algorithm
#' (`fb_lrcc_R`).
#'
#' @inheritParams lr_ind_dist
#'
#' @return A data frame with two columns `LR` and `prob`.
#' @export
lr_cc_dist <- function(n, alpha = 0.05, prune_threshold = 1e-15) {
  if (exists("fb_lrcc_fastcpp", mode = "function")) {
    tryCatch(
      fb_lrcc_fastcpp(n, alpha, prune_threshold),
      error = function(e) fb_lrcc_R(n, alpha, prune_threshold)
    )
  } else {
    fb_lrcc_R(n, alpha, prune_threshold)
  }
}

