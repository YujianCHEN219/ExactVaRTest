## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
agg_probs <- function(res) {
  tapply(res$prob, res$LR, sum, simplify = TRUE)
}

collapse_close <- function(res, tol = 1e-12) {
  o  <- order(res$LR)
  LR <- res$LR[o]; P <- res$prob[o]
  keep <- c(TRUE, abs(diff(LR)) > tol)
  grp  <- cumsum(keep)
  LRc  <- tapply(LR, grp, `[`, 1L)
  Pc   <- tapply(P,  grp, sum)
  list(LR = as.numeric(LRc), prob = as.numeric(Pc) / sum(Pc))
}

strip_tiny <- function(res, prob_tol = 1e-12) {
  keep <- res$prob > prob_tol
  if (all(keep)) return(res)
  res <- list(LR = res$LR[keep], prob = res$prob[keep])
  res$prob <- res$prob / sum(res$prob)
  res
}

check_lr_dist <- function(gen_cpp, gen_R, n_vec = c(5, 10),
                          alpha = 0.05, tol = 1e-8) {
  for (n in n_vec) {
    dist_cpp <- gen_cpp(n, alpha)
    dist_R   <- gen_R(n, alpha)
    
    ## basic shape
    expect_gt(length(dist_cpp$LR), 0)
    
    ## probabilities must be finite, non-negative and sum to 1
    expect_true(all(is.finite(dist_cpp$prob)))
    expect_true(all(dist_cpp$prob >= 0))
    expect_equal(sum(dist_cpp$prob), 1, tolerance = tol)
    
    ## compare C++ vs. R after aggregating duplicate LR values
    a_cpp <- agg_probs(dist_cpp)
    a_R   <- agg_probs(dist_R)
    
    expect_equal(sort(as.numeric(names(a_cpp))),
                 sort(as.numeric(names(a_R))),
                 tolerance = tol)
    
    a_R <- a_R[names(a_cpp)]               
    expect_equal(unname(a_cpp), unname(a_R), tolerance = tol)
  }
}

## tolerance settings
TOL        <- 1e-8      # numeric comparison
COLL_TOL   <- 1e-12     # collapse nearly-equal LR
PROB_TOL   <- 1e-12     # drop negligible probability

## ------------------------------------------------------------------
## LR_ind : C++ vs. R  (n = 40)
## ------------------------------------------------------------------
test_that("lr_ind_dist – C++ and R engines numerically identical", {
  n     <- 40
  alpha <- 0.05
  
  res_cpp <- lr_ind_dist(n, alpha)
  res_R   <- fb_lrind_R(n, alpha)
  
  res_cpp <- strip_tiny(collapse_close(res_cpp, COLL_TOL), PROB_TOL)
  res_R   <- strip_tiny(collapse_close(res_R,   COLL_TOL), PROB_TOL)
  
  expect_equal(res_cpp$LR,   res_R$LR,   tolerance = TOL)
  expect_equal(res_cpp$prob, res_R$prob, tolerance = TOL)
  
  expect_true(all(is.finite(res_cpp$prob)))
  expect_true(all(res_cpp$prob >= 0))
  expect_equal(sum(res_cpp$prob), 1, tolerance = TOL)
})

## ------------------------------------------------------------------
## LR_cc : C++ vs. R  (n = 40)
## ------------------------------------------------------------------
test_that("lr_cc_dist – C++ and R engines numerically identical", {
  n     <- 40
  alpha <- 0.05
  
  res_cpp <- lr_cc_dist(n, alpha)
  res_R   <- fb_lrcc_R(n, alpha)
  
  res_cpp <- strip_tiny(collapse_close(res_cpp, COLL_TOL), PROB_TOL)
  res_R   <- strip_tiny(collapse_close(res_R,   COLL_TOL), PROB_TOL)
  
  expect_equal(res_cpp$LR,   res_R$LR,   tolerance = TOL)
  expect_equal(res_cpp$prob, res_R$prob, tolerance = TOL)
  
  expect_true(all(is.finite(res_cpp$prob)))
  expect_true(all(res_cpp$prob >= 0))
  expect_equal(sum(res_cpp$prob), 1, tolerance = TOL)
})

## ------------------------------------------------------------------
## Lightweight finite-sample sanity checks (shared for ind & cc)
## ------------------------------------------------------------------
test_that("lr_ind_dist / lr_cc_dist give valid finite-sample distributions", {
  skip_on_cran()
  skip_if_not_installed("ExactVaRTest")
  
  check_lr_dist(lr_ind_dist, fb_lrind_R, n_vec = c(5, 10), alpha = 0.05, tol = TOL)
  check_lr_dist(lr_cc_dist,  fb_lrcc_R,  n_vec = c(5, 10), alpha = 0.05, tol = TOL)
})
