## ------------------------------------------------------------------
## Helpers
## ------------------------------------------------------------------
agg_probs <- function(res) {
  tapply(res$prob, res$LR, sum, simplify = TRUE)
}

collapse_close <- function(res, digits = 10) {
  res$LR <- round(res$LR, digits)
  agg <- tapply(res$prob, res$LR, sum)
  list(LR = as.numeric(names(agg)),
       prob = as.numeric(agg) / sum(agg))
}

strip_tiny <- function(res, prob_tol = 1e-12) {
  keep <- res$prob > prob_tol
  if (all(keep)) return(res)
  res <- list(LR = res$LR[keep], prob = res$prob[keep])
  res$prob <- res$prob / sum(res$prob)
  res
}

## tolerance settings
TOL      <- 1e-8      # numeric comparison
PROB_TOL <- 1e-12     # drop negligible probability

check_lr_dist <- function(gen_cpp, gen_R,
                          n_vec  = c(5, 10),
                          alpha  = 0.05,
                          tol    = 1e-8) {
  for (n in n_vec) {
    dist_cpp <- strip_tiny(collapse_close(gen_cpp(n, alpha)), PROB_TOL)
    dist_R   <- strip_tiny(collapse_close(gen_R  (n, alpha)), PROB_TOL)
    
    allLR <- sort(unique(c(dist_cpp$LR, dist_R$LR)))
    p_cpp <- dist_cpp$prob[match(allLR, dist_cpp$LR)]
    p_R   <- dist_R$prob  [match(allLR, dist_R$LR)]
    p_cpp[is.na(p_cpp)] <- 0
    p_R  [is.na(p_R)]   <- 0
    
    expect_true(all(is.finite(p_cpp)))
    expect_true(all(p_cpp >= 0))
    expect_equal(sum(p_cpp), 1, tolerance = tol)
    
    expect_equal(cumsum(p_cpp), cumsum(p_R), tolerance = tol)
  }
}

## ------------------------------------------------------------------
## LR_ind : C++ vs. R  (n = 40)
## ------------------------------------------------------------------
test_that("lr_ind_dist – C++ and R engines numerically identical", {
  n     <- 40
  alpha <- 0.05
  
  res_cpp <- strip_tiny(collapse_close(lr_ind_dist(n, alpha)), PROB_TOL)
  res_R   <- strip_tiny(collapse_close(fb_lrind_R (n, alpha)), PROB_TOL)
  
  expect_equal(cumsum(res_cpp$prob), cumsum(res_R$prob), tolerance = TOL)
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
  
  res_cpp <- strip_tiny(collapse_close(lr_cc_dist(n, alpha)), PROB_TOL)
  res_R   <- strip_tiny(collapse_close(fb_lrcc_R (n, alpha)), PROB_TOL)
  
  expect_equal(cumsum(res_cpp$prob), cumsum(res_R$prob), tolerance = TOL)
  expect_true(all(is.finite(res_cpp$prob)))
  expect_true(all(res_cpp$prob >= 0))
  expect_equal(sum(res_cpp$prob), 1, tolerance = TOL)
})

## ------------------------------------------------------------------
## Lightweight finite-sample sanity checks
## ------------------------------------------------------------------
test_that("lr_ind_dist / lr_cc_dist give valid finite-sample distributions", {
  skip_on_cran()
  skip_if_not_installed("ExactVaRTest")
  
  check_lr_dist(lr_ind_dist, fb_lrind_R, n_vec = c(5, 10),
                alpha = 0.05, tol = TOL)
  check_lr_dist(lr_cc_dist,  fb_lrcc_R,  n_vec = c(5, 10),
                alpha = 0.05, tol = TOL)
})
