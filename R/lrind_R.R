## ------------------------------------------------------------------
##  Pure-R exact distribution of Christoffersen's LR_ind
## ------------------------------------------------------------------
fb_lrind_R <- function(n,
                       alpha           = 0.05,
                       prune_threshold = 1e-15) {
  
  if (n < 2)
    return(list(LR = 0, prob = 1))
  
  combine_states <- function(mat) {
    if (!nrow(mat)) return(mat)
    o   <- do.call(order, lapply(seq_len(5), function(i) mat[, i]))
    mat <- mat[o, , drop = FALSE]
    idx <- c(which(rowSums(abs(diff(mat[, 1:5]))) != 0), nrow(mat))
    out <- mat[idx, , drop = FALSE]
    start <- 1L
    for (i in seq_along(idx)) {
      s <- start; e <- idx[i]
      if (s < e) out[i, 6] <- sum(mat[s:e, 6])
      start <- e + 1L
    }
    out
  }
  
  S <- matrix(c(0, 0, 0, 0, 0, 1 - alpha,
                1, 0, 0, 0, 0, alpha),
              nrow = 2, byrow = TRUE)
  
  for (k in seq_len(n - 1)) {
    S <- S[S[, 6] >= prune_threshold, , drop = FALSE]
    if (!nrow(S)) break
    
    tmp0 <- S
    tmp1 <- S
    
    tmp0[, 2] <- tmp0[, 2] + (tmp0[, 1] == 0)
    tmp0[, 3] <- tmp0[, 3] + (tmp0[, 1] == 1)
    tmp0[, 1] <- 0
    tmp0[, 6] <- tmp0[, 6] * (1 - alpha)
    tmp0 <- tmp0[tmp0[, 6] >= prune_threshold, , drop = FALSE]
    
    tmp1[, 4] <- tmp1[, 4] + (tmp1[, 1] == 0)
    tmp1[, 5] <- tmp1[, 5] + (tmp1[, 1] == 1)
    tmp1[, 1] <- 1
    tmp1[, 6] <- tmp1[, 6] * alpha
    tmp1 <- tmp1[tmp1[, 6] >= prune_threshold, , drop = FALSE]
    
    S <- combine_states(rbind(tmp0, tmp1))
  }
  
  S <- S[S[, 6] >= prune_threshold, , drop = FALSE]
  if (!nrow(S))
    return(list(LR = numeric(0), prob = numeric(0)))
  
  T0  <- S[, 2] + S[, 3]
  T1  <- S[, 4] + S[, 5]
  pH  <- if (n > 1) T1 / (n - 1) else 0
  num <- T0 * safe_log(1 - pH) + T1 * safe_log(pH)
  
  pi01 <- ifelse((S[, 2] + S[, 4]) > 0,
                 S[, 4] / (S[, 2] + S[, 4]), 1)
  pi11 <- ifelse((S[, 3] + S[, 5]) > 0,
                 S[, 5] / (S[, 3] + S[, 5]), 1)
  den  <- S[, 2] * safe_log(1 - pi01) + S[, 4] * safe_log(pi01) +
    S[, 3] * safe_log(1 - pi11) + S[, 5] * safe_log(pi11)
  
  LR <- -2 * (num - den)
  
  keep <- is.finite(LR)
  LR   <- LR[keep]
  S    <- S[keep, , drop = FALSE]
  if (!length(LR))
    return(list(LR = numeric(0), prob = numeric(0)))
  
  out <- cbind(LR, S[, 6])
  out <- out[order(out[, 1]), , drop = FALSE]
  
  idx   <- c(which(diff(out[, 1]) != 0), nrow(out))
  final <- out[idx, , drop = FALSE]
  start <- 1L
  for (i in seq_along(idx)) {
    s <- start; e <- idx[i]
    if (s < e) final[i, 2] <- sum(out[s:e, 2])
    start <- e + 1L
  }
  
  s <- sum(final[, 2])
  if (s == 0)
    return(list(LR = numeric(0), prob = numeric(0)))
  
  final[, 2] <- final[, 2] / s
  list(LR = final[, 1], prob = final[, 2])
}
