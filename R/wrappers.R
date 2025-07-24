# -----------------------------------------------------------------------
#  Wrappers that auto‑select the fastest available engine
# -----------------------------------------------------------------------

#' @importFrom stats dbinom

#  Numeric constant used by lr_*_stat()
EPS <- 1e-15

#' Exact LR_ind distribution (auto‑select engine)
#'
#' Returns the finite‑sample distribution of Christoffersen’s independence
#' statistic \eqn{LR_{\mathrm{ind}}}.
#'
#' @inheritParams lr_uc_dist
#' @param prune_threshold Probability below which states are pruned by the
#'   dynamic‑programming recursion.
#' @return A list with elements \code{LR} and \code{prob}.
#'
#' @examples
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

#' Exact LR_cc (and LR_uc) distribution (auto‑select engine)
#'
#' Returns the finite‑sample distribution of Christoffersen’s conditional‑coverage
#' statistic \eqn{LR_{\mathrm{cc}}}.  The returned list also includes the matching
#' unconditional‑coverage distribution \eqn{LR_{\mathrm{uc}}}, produced by the same
#' dynamic‑programming run.
#'
#' @inheritParams lr_ind_dist
#' @return A list with elements \code{LR_cc}, \code{prob_cc}, \code{LR_uc},
#'   \code{prob_uc}.
#'
#' @examples
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

#' Exact LR_uc distribution (closed‑form binomial)
#'
#' @param n Integer sample size (\eqn{n \ge 1}).
#' @param alpha Exception probability \eqn{\alpha \in (0,1)}.
#' @return A list with elements \code{LR} and \code{prob}.
#'
#' @examples
#' lr_uc_dist(8, 0.01)
#' @export
lr_uc_dist <- function(n, alpha = 0.05) {
  c1   <- 0:n
  prob <- dbinom(c1, n, alpha)
  p_   <- max(min(alpha, 1 - EPS), EPS)
  phat <- pmax(pmin(c1 / n, 1 - EPS), EPS)
  LR   <- -2 * ( c1 * log(p_) + (n - c1) * log(1 - p_) -
                   c1 * log(phat) - (n - c1) * log(1 - phat) )
  list(LR = LR, prob = prob)
}

# -----------------------------------------------------------------------
#  Low‑level helpers: LR statistics (pure R)
# -----------------------------------------------------------------------

#' Christoffersen LR_ind statistic
#'
#' @param x     0/1 exception series.
#' @param alpha Exception probability.
#' @return Numeric LR_ind statistic.
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
  num  <- T0 * log(pmax(1 - pHat, EPS)) + T1 * log(pmax(pHat, EPS))
  pi01 <- if ((T00 + T01) > 0) T01 / (T00 + T01) else 1
  pi11 <- if ((T10 + T11) > 0) T11 / (T10 + T11) else 1
  den  <- T00 * log(pmax(1 - pi01, EPS)) + T01 * log(pmax(pi01, EPS)) +
    T10 * log(pmax(1 - pi11, EPS)) + T11 * log(pmax(pi11, EPS))
  -2 * (num - den)
}

#' Christoffersen LR_cc statistic
#'
#' @inheritParams lr_ind_stat
#' @return Numeric LR_cc statistic.
#' @export
lr_cc_stat <- function(x, alpha = 0.05) {
  n  <- length(x)
  c1 <- sum(x)
  p_   <- max(min(alpha, 1 - EPS), EPS)
  phat <- if (c1 == 0) 0 else if (c1 == n) 1 else c1 / n
  ph_  <- max(min(phat, 1 - EPS), EPS)
  lr_uc <- -2 * ( c1 * log(p_) + (n - c1) * log(1 - p_) -
                    c1 * log(ph_) - (n - c1) * log(1 - ph_) )
  lr_uc + lr_ind_stat(x, alpha)
}


#' Christoffersen LR_uc statistic
#' 
#' @param x     0/1 exception series.
#' @param alpha Exception probability.
#' @return Numeric LR_uc statistic.
#' @export              
lr_uc_stat <- function(x, alpha = 0.05) {
  n  <- length(x)
  c1 <- sum(x)
  p_   <- max(min(alpha, 1 - EPS), EPS)
  phat <- if (c1 == 0) 0 else if (c1 == n) 1 else c1 / n
  ph_  <- max(min(phat, 1 - EPS), EPS)
  -2 * ( c1 * log(p_) + (n - c1) * log(1 - p_) -
           c1 * log(ph_) - (n - c1) * log(1 - ph_) )
}

# -----------------------------------------------------------------------
#  Convenience helpers: exact p-values
# -----------------------------------------------------------------------

#' Exact p-value for LR_uc
#'
#' @param lr_obs Observed LR_uc statistic.
#' @param n      Sample size.
#' @param alpha  Exception probability.
#' @return Numeric p-value.
#' @export
pval_lr_uc <- function(lr_obs, n, alpha = 0.05) {
  dist <- lr_uc_dist(n, alpha)
  if (!length(dist$LR)) return(NA_real_)
  sum(dist$prob[dist$LR >= lr_obs])
}

#' Exact p-value for LR_ind
#'
#' @param lr_obs           Observed LR_ind statistic.
#' @param n                Sample size.
#' @param alpha            Exception probability.
#' @param prune_threshold  State-pruning threshold for DP engine.
#' @return Numeric p-value.
#' @export
pval_lr_ind <- function(lr_obs, n, alpha = 0.05, prune_threshold = 1e-15) {
  dist <- lr_ind_dist(n, alpha, prune_threshold)
  if (!length(dist$LR)) return(NA_real_)
  sum(dist$prob[dist$LR >= lr_obs])
}

#' Exact p-value for LR_cc
#'
#' @param lr_obs           Observed LR_cc statistic.
#' @param n                Sample size.
#' @param alpha            Exception probability.
#' @param prune_threshold  State-pruning threshold for DP engine.
#' @return Numeric p-value.
#' @export
pval_lr_cc <- function(lr_obs, n, alpha = 0.05, prune_threshold = 1e-15) {
  dist <- lr_cc_dist(n, alpha, prune_threshold)
  if (!length(dist$LR_cc)) return(NA_real_)
  sum(dist$prob_cc[dist$LR_cc >= lr_obs])
}


# -----------------------------------------------------------------------
#  Print method for back‑test results
# -----------------------------------------------------------------------

#' Print method for ExactVaRBacktest
#'
#' @param x      An \code{ExactVaRBacktest} object.
#' @param digits Number of digits to print.
#' @param ...    Unused.
#' @method print ExactVaRBacktest
#' @export
print.ExactVaRBacktest <- function(x,
                                   digits = max(3L, getOption("digits") - 3L),
                                   ...) {
  digits <- max(1L, as.integer(digits))
  test_name <- switch(
    x$type,
    uc  = "Unconditional coverage (LR_uc)",        
    ind = "Independence (LR_ind)",
    cc  = "Conditional coverage (LR_cc)"        
  )
  
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
# -----------------------------------------------------------------------
#  Single‑test back‑test helper
# -----------------------------------------------------------------------

#' Exact finite‑sample back‑test for a VaR exception series
#'
#' @inheritParams lr_ind_stat
#' @param type `"uc"`, `"ind"` or `"cc"`.
#' @param sig  Significance level (default `0.05`).
#' @param prune_threshold Passed to the dynamic‑programming engine.
#' @return An object of class `"ExactVaRBacktest"`.
#'
#' @examples
#' set.seed(123)
#' x <- rbinom(250, 1, 0.01)
#' backtest_lr(x, alpha = 0.01, type = "uc")
#' @export
backtest_lr <- function(x,
                        alpha = 0.05,
                        type  = c("uc", "ind", "cc"),
                        sig   = 0.05,
                        prune_threshold = 1e-15) {
  
  type <- match.arg(type)
  n    <- length(x)
  if (n < 1) stop("Series 'x' must have positive length.")
  
  if (type == "uc") {
    stat <- lr_uc_stat(x, alpha)
    pval <- pval_lr_uc(stat, n, alpha)
  } else if (type == "ind") {
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
#  All in one three‑test back‑tester
# -----------------------------------------------------------------------

#' Exact UC/IND/CC back‑tests in one call
#'
#' @inheritParams lr_ind_stat
#' @param sig  Significance level (default `0.05`).
#' @param prune_threshold Passed to the dynamic programming engine.
#' @return An object of class `"ExactVaRBacktestAll"`.
#'
#' @examples
#' set.seed(1)
#' x <- rbinom(300, 1, 0.02)
#' backtest_all(x, alpha = 0.02)
#' @export
backtest_all <- function(x,
                         alpha = 0.05,
                         sig   = 0.05,
                         prune_threshold = 1e-15) {
  
  n <- length(x)
  if (n < 1) stop("Series 'x' must have positive length.")
  
  dist_cc  <- lr_cc_dist(n, alpha, prune_threshold)
  dist_ind <- lr_ind_dist(n, alpha, prune_threshold)
  
  stat_uc  <- lr_uc_stat(x, alpha)
  stat_ind <- lr_ind_stat(x, alpha)
  stat_cc  <- lr_cc_stat(x, alpha)
  
  p_uc  <- sum(dist_cc$prob_uc[dist_cc$LR_uc >= stat_uc])
  p_cc  <- sum(dist_cc$prob_cc[dist_cc$LR_cc >= stat_cc])
  p_ind <- sum(dist_ind$prob   [dist_ind$LR    >= stat_ind])
  
  obj <- list(
    uc  = list(stat = stat_uc,  pval = p_uc,  reject = p_uc  < sig),
    ind = list(stat = stat_ind, pval = p_ind, reject = p_ind < sig),
    cc  = list(stat = stat_cc,  pval = p_cc,  reject = p_cc  < sig),
    sig   = sig,
    alpha = alpha,
    n     = n
  )
  class(obj) <- "ExactVaRBacktestAll"
  obj
}

# -----------------------------------------------------------------------
#  Print method for ExactVaRBacktestAll
# -----------------------------------------------------------------------

#' @export
print.ExactVaRBacktestAll <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  digits <- max(1L, as.integer(digits))
  
  # header printed once
  cat("Exact finite sample backtest\n\n",
      sprintf("%-15s: %d\n",  "Sample size",   x$n),
      sprintf("%-15s: %.4f\n","Model alpha",   x$alpha),
      sprintf("%-15s: %.4f\n","Signif. level", x$sig),
      sep = "")
  
  print_one <- function(res, title) {
    stat_fmt <- formatC(res$stat, format = "f", digits = digits)
    pval_fmt <- formatC(res$pval, format = "f", digits = digits)
    lvl_fmt  <- formatC(x$sig * 100, format = "f", digits = 2)
    decision <- if (res$reject)
      paste0("REJECT null at ", lvl_fmt, "% level")
    else
      paste0("fail to reject at ", lvl_fmt, "% level")
    
    cat("--------------------------------\n",
        sprintf("%-15s: %s\n",  "Test",          title),
        sprintf("%-15s: %s\n",  "LR statistic",  stat_fmt),
        sprintf("%-15s: %s\n",  "Exact p-value", pval_fmt),
        sprintf("%-15s: %s\n",  "Decision",      decision),
        sep = "")
  }
  
  print_one(x$uc,  "Unconditional coverage (LR_uc)")
  print_one(x$ind, "Independence (LR_ind)")
  print_one(x$cc,  "Conditional coverage (LR_cc)")
  
  invisible(x)
}
