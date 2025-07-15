## ------------------------------------------------------------------
##  Pure-R exact distribution of Christoffersen's LR_cc
## ------------------------------------------------------------------
fb_lrcc_R <- function(n,
                      alpha           = 0.05,
                      prune_threshold = 1e-15) {
  
  if (n < 2)
    return(list(LR = 0, prob = 1))
  
  LRuc_count <- function(c1, n, p) {
    p_   <- max(min(p, 1 - EPS), EPS)
    phat <- if (c1 == 0) 0 else if (c1 == n) 1 else c1 / n
    ph_  <- max(min(phat, 1 - EPS), EPS)
    -2 * (c1 * safe_log(p_) + (n - c1) * safe_log(1 - p_) -
            c1 * safe_log(ph_) - (n - c1) * safe_log(1 - ph_))
  }
  
  combine_states <- function(mat) {
    if (!nrow(mat)) return(mat)
    o   <- do.call(order, lapply(seq_len(6), function(i) mat[, i]))
    mat <- mat[o, , drop = FALSE]
    idx <- c(which(rowSums(abs(diff(mat[, 1:6]))) != 0), nrow(mat))
    out <- mat[idx, , drop = FALSE]
    start <- 1L
    for (i in seq_along(idx)) {
      s <- start; e <- idx[i]
      if (s < e) out[i, 7] <- sum(mat[s:e, 7])
      start <- e + 1L
    }
    out
  }
  
  S <- matrix(c(0, 0, 0, 0, 0, 0, 1 - alpha,
                1, 1, 0, 0, 0, 0, alpha),
              nrow = 2, byrow = TRUE)
  
  for (k in seq_len(n - 1)) {
    S <- S[S[, 7] >= prune_threshold, , drop = FALSE]
    if (!nrow(S)) break
    
    tmp0 <- S
    tmp1 <- S
    
    tmp0[, 3] <- tmp0[, 3] + (tmp0[, 1] == 0)
    tmp0[, 4] <- tmp0[, 4] + (tmp0[, 1] == 1)
    tmp0[, 1] <- 0
    tmp0[, 7] <- tmp0[, 7] * (1 - alpha)
    tmp0 <- tmp0[tmp0[, 7] >= prune_threshold, , drop = FALSE]
    
    tmp1[, 5] <- tmp1[, 5] + (tmp1[, 1] == 0)
    tmp1[, 6] <- tmp1[, 6] + (tmp1[, 1] == 1)
    tmp1[, 2] <- tmp1[, 2] + 1
    tmp1[, 1] <- 1
    tmp1[, 7] <- tmp1[, 7] * alpha
    tmp1 <- tmp1[tmp1[, 7] >= prune_threshold, , drop = FALSE]
    
    S <- combine_states(rbind(tmp0, tmp1))
  }
  
  S <- S[S[, 7] >= prune_threshold, , drop = FALSE]
  if (!nrow(S))
    return(list(LR = numeric(0), prob = numeric(0)))
  
  c1  <- S[, 2]
  T00 <- S[, 3]; T10 <- S[, 4]; T01 <- S[, 5]; T11 <- S[, 6]
  
  LRuc <- vapply(c1, LRuc_count, numeric(1), n = n, p = alpha)
  
  T0   <- T00 + T10
  T1   <- T01 + T11
  pHat <- if (n > 1) T1 / (n - 1) else 0
  num  <- T0 * safe_log(1 - pHat) + T1 * safe_log(pHat)
  
  pi01 <- ifelse((T00 + T01) > 0,
                 T01 / (T00 + T01), 1)
  pi11 <- ifelse((T10 + T11) > 0,
                 T11 / (T10 + T11), 1)
  den  <- T00 * safe_log(1 - pi01) + T01 * safe_log(pi01) +
    T10 * safe_log(1 - pi11) + T11 * safe_log(pi11)
  
  LR <- LRuc - 2 * (num - den)
  keep <- is.finite(LR)
  LR   <- LR[keep]
  S    <- S[keep, , drop = FALSE]
  
  out <- cbind(LR, S[, 7])
  out <- out[order(out[, 1]), , drop = FALSE]
  
  idx <- c(which(diff(out[, 1]) != 0), nrow(out))
  res <- out[idx, , drop = FALSE]
  start <- 1L
  for (i in seq_along(idx)) {
    s <- start; e <- idx[i]
    if (s < e) res[i, 2] <- sum(out[s:e, 2])
    start <- e + 1L
  }
  
  s <- sum(res[, 2])
  if (s == 0)
    return(list(LR = numeric(0), prob = numeric(0)))
  
  res[, 2] <- res[, 2] / s
  list(LR = res[, 1], prob = res[, 2])
}
