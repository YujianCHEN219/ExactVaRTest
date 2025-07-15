# -----------------------------------------------------------------------
#  Wrappers that auto-select the fastest available engine
# -----------------------------------------------------------------------

#' Exact LR_ind distribution (auto-select engine)
#'
#' Returns the finite-sample distribution of Christoffersen’s independence
#' statistic \eqn{LR_{\mathrm{ind}}}.
#'
#' @param n Integer sample size (\eqn{n \ge 1}).
#' @param alpha Exception probability \eqn{\alpha \in (0,1)}.
#' @param prune_threshold Probability below which states are pruned by the
#'   dynamic-programming recursion.
#'
#' @return A data frame with two columns  
#'   \describe{  
#'     \item{LR}{unique LR values (ascending)}  
#'     \item{prob}{exact probabilities (sum to 1)}  
#'   }
#'
#' @examples
#' ## exact LR_ind distribution for n = 8, alpha = 0.05
#' lr_ind_dist(8, 0.05)
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
#' Returns the finite-sample distribution of Christoffersen’s
#' conditional-coverage statistic \eqn{LR_{\mathrm{cc}}}.
#'
#' @inheritParams lr_ind_dist
#'
#' @return A data frame with columns \code{LR} and \code{prob}.
#'
#' @examples
#' ## exact LR_cc distribution for n = 8, alpha = 0.05
#' lr_cc_dist(8, 0.05)
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

# -----------------------------------------------------------------------
#  Low-level helpers: LR statistics (pure R)
# -----------------------------------------------------------------------

EPS_ <- 1e-15                       # local safe-log constant

#' Christoffersen LR_ind statistic
#'
#' @param x     0/1 exception series.
#' @param alpha Exception probability.
#'
#' @return Numeric LR_ind statistic.
#'
#' @examples
#' set.seed(1)
#' x <- rbinom(50, 1, 0.05)
#' lr_ind_stat(x, 0.05)
#' @export
lr_ind_stat <- function(x, alpha = 0.05) {
  n <- length(x)
  if (n < 2) return(0)
  T00 <- T01 <- T10 <- T11 <- 0L
  for (t in 2:n) {
    if (x[t-1L] == 0L && x[t] == 0L)      T00 <- T00 + 1L else
      if (x[t-1L] == 0L && x[t] == 1L)      T01 <- T01 + 1L else
        if (x[t-1L] == 1L && x[t] == 0L)      T10 <- T10 + 1L else
          T11 <- T11 + 1L
  }
  T0   <- T00 + T10
  T1   <- T01 + T11
  pHat <- if (n > 1) T1 / (n - 1) else 0
  num  <- T0 * log(pmax(1 - pHat, EPS_)) + T1 * log(pmax(pHat,  EPS_))
  pi01 <- if ((T00 + T01) > 0) T01 / (T00 + T01) else 1
  pi11 <- if ((T10 + T11) > 0) T11 / (T10 + T11) else 1
  den  <- T00 * log(pmax(1 - pi01, EPS_)) + T01 * log(pmax(pi01, EPS_)) +
    T10 * log(pmax(1 - pi11, EPS_)) + T11 * log(pmax(pi11, EPS_))
  -2 * (num - den)
}

#' Christoffersen LR_cc statistic
#'
#' Computes \eqn{LR_{\mathrm{cc}} = LR_{\mathrm{uc}} + LR_{\mathrm{ind}}}.
#'
#' @inheritParams lr_ind_stat
#'
#' @return Numeric LR_cc statistic.
#'
#' @examples
#' set.seed(1)
#' x <- rbinom(50, 1, 0.05)
#' lr_cc_stat(x, 0.05)
#' @export
lr_cc_stat <- function(x, alpha = 0.05) {
  n  <- length(x)
  c1 <- sum(x)
  
  ## LR_uc
  p_   <- max(min(alpha, 1 - EPS_), EPS_)
  phat <- if (c1 == 0) 0 else if (c1 == n) 1 else c1 / n
  ph_  <- max(min(phat, 1 - EPS_), EPS_)
  
  lr_uc <- -2 * ( c1 * log(p_)      + (n - c1) * log(1 - p_) -
                    c1 * log(ph_)     - (n - c1) * log(1 - ph_) )
  
  lr_uc + lr_ind_stat(x, alpha)
}

# -----------------------------------------------------------------------
#  Convenience helpers: exact p-values
# -----------------------------------------------------------------------

#' Exact p-value for LR_ind
#'
#' @param lr_obs Observed \eqn{LR_{\mathrm{ind}}} value.
#' @inheritParams lr_ind_dist
#'
#' @return Numeric p-value.
#'
#' @examples
#' set.seed(1)
#' x  <- rbinom(250, 1, 0.05)
#' lr <- lr_ind_stat(x, 0.05)
#' pval_lr_ind(lr, length(x), 0.05)
#' @export
pval_lr_ind <- function(lr_obs, n, alpha = 0.05, prune_threshold = 1e-15) {
  dist <- lr_ind_dist(n, alpha, prune_threshold)
  if (!length(dist$LR)) return(NA_real_)
  sum(dist$prob[dist$LR >= lr_obs])
}

#' Exact p-value for LR_cc
#'
#' @param lr_obs Observed \eqn{LR_{\mathrm{cc}}} value.
#' @inheritParams lr_ind_dist
#'
#' @return Numeric p-value.
#'
#' @examples
#' set.seed(1)
#' x  <- rbinom(250, 1, 0.05)
#' lr <- lr_cc_stat(x, 0.05)
#' pval_lr_cc(lr, length(x), 0.05)
#' @export
pval_lr_cc <- function(lr_obs, n, alpha = 0.05, prune_threshold = 1e-15) {
  dist <- lr_cc_dist(n, alpha, prune_threshold)
  if (!length(dist$LR)) return(NA_real_)
  sum(dist$prob[dist$LR >= lr_obs])
}

# -----------------------------------------------------------------------
#  One-shot back-test helper
# -----------------------------------------------------------------------

#' Exact finite-sample back-test for a VaR exception series
#'
#' Computes the LR statistic, its exact finite-sample p-value and a
#' reject/accept decision **in one call**.  
#' The returned object has class `"ExactVaRBacktest"` with a
#' tailor-made `print()` method.
#'
#' @param x 0/1 exception series.
#' @param alpha Exception probability used by the VaR model.
#' @param type `"ind"` (independence) or `"cc"` (conditional coverage).
#'   Partial matching allowed.
#' @param sig Significance level for the decision rule. Default `0.05`.
#' @param prune_threshold Passed to the dynamic-programming engine.
#'
#' @return An object of class `"ExactVaRBacktest"` (a named list).
#'
#' @examples
#' set.seed(42)
#' x <- rbinom(300, 1, 0.05)
#' backtest_lr(x, alpha = 0.05, type = "cc")
#' @export
backtest_lr <- function(x,
                        alpha = 0.05,
                        type  = c("ind", "cc"),
                        sig   = 0.05,
                        prune_threshold = 1e-15) {
  
  type <- match.arg(type)
  n    <- length(x)
  if (n < 1) stop("Series 'x' must have positive length.")
  
  if (type == "ind") {
    stat <- lr_ind_stat(x, alpha)
    pval <- pval_lr_ind(stat, n, alpha, prune_threshold)
  } else {
    stat <- lr_cc_stat(x, alpha)
    pval <- pval_lr_cc(stat, n, alpha, prune_threshold)
  }
  
  obj <- list(
    stat   = stat,
    pval   = pval,
    reject = (pval < sig),
    type   = type,
    alpha  = alpha,
    sig    = sig,
    n      = n
  )
  class(obj) <- "ExactVaRBacktest"
  obj
}

# -----------------------------------------------------------------------
# Printer
# -----------------------------------------------------------------------

#' @export
print.ExactVaRBacktest <- function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  
  digits <- max(1L, as.integer(digits))          # keep at least 1 digit
  
  test_name <- if (x$type == "ind")
    "Independence (LR_ind)"
  else
    "Conditional-coverage (LR_cc)"
  
  stat_fmt <- formatC(x$stat, format = "f", digits = digits)
  pval_fmt <- formatC(x$pval, format = "f", digits = digits)
  lvl_fmt  <- formatC(x$sig * 100, format = "f", digits = 2)
  
  decision <- if (x$reject)
    paste0("REJECT null at ", lvl_fmt, "% level")
  else
    paste0("fail to reject at ", lvl_fmt, "% level")
  
  cat("Exact finite-sample back-test\n",
      "--------------------------------\n",
      sprintf("%-15s: %s\n",  "Test",           test_name),
      sprintf("%-15s: %d\n",  "Sample size",    x$n),
      sprintf("%-15s: %.4f\n","Model alpha",    x$alpha),
      sprintf("%-15s: %.4f\n","Signif. level",  x$sig),
      sprintf("%-15s: %s\n",  "LR statistic",   stat_fmt),
      sprintf("%-15s: %s\n",  "Exact p-value",  pval_fmt),
      sprintf("%-15s: %s\n",  "Decision",       decision),
      sep = "")
  invisible(x)
}

