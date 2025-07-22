## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ExactVaRTest)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(bench)
library(dplyr)
library(tidyr)
library(purrr)
library(knitr)

n_vec     <- c(50, 100, 250, 500, 750, 1000)
alpha_vec <- c(0.01, 0.025, 0.05)

grid <- expand.grid(n = n_vec, alpha = alpha_vec, KEEP.OUT.ATTRS = FALSE)

bench_one <- function(n, alpha, fun) {
  bm <- bench::mark(
    fun(n = n, alpha = alpha),
    iterations = 3,
    check = FALSE
  )
  tibble(
    n        = n,
    alpha    = alpha,
    median_s = as.numeric(bm$median)
  )
}

timings_ind <- pmap_dfr(grid, bench_one, fun = lr_ind_dist)

timings_cc  <- pmap_dfr(grid, bench_one, fun = lr_cc_dist)

fmt_wide <- function(df) {
  df %>%
    mutate(alpha = sprintf("alpha = %.3f", alpha)) %>%
    pivot_wider(names_from = alpha, values_from = median_s) %>%
    arrange(n)
}

table_ind <- fmt_wide(timings_ind)
table_cc  <- fmt_wide(timings_cc)

kable(
  table_ind,
  digits = 3,
  col.names = c("n", "α = 0.010", "α = 0.025", "α = 0.050"),
  caption   = "Median run-time (seconds) for lr_ind_dist()"
)

kable(
  table_cc,
  digits = 3,
  col.names = c("n", "α = 0.010", "α = 0.025", "α = 0.050"),
  caption   = "Median run-time (seconds) for lr_cc_dist()"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(xts)

alpha  <- 0.01
window <- 250

local_file <- "inst/extdata/GSPC_close.rds"
file_path  <- if (file.exists(local_file)) local_file else
              system.file("extdata", "GSPC_close.rds", package = "ExactVaRTest")

px  <- readRDS(file_path)
ret <- diff(log(px), lag = 1)

VaR <- rollapply(
  ret, window,
  function(r) quantile(r, alpha, na.rm = TRUE),
  by.column = FALSE, align = "right"
)

keep <- complete.cases(ret, VaR)
ret  <- ret[keep]
VaR  <- coredata(VaR[keep])

x <- ifelse(coredata(ret) < VaR, 1L, 0L)

cat("Series length :", length(x), "\n")
cat("Exception rate:", mean(x), "\n\n")

bench_res <- bench::mark(
  LR_ind = backtest_lr(x, alpha, "ind"),
  LR_cc  = backtest_lr(x, alpha, "cc"),
  iterations = 10,
  check      = FALSE
)

suppressWarnings(
  print(bench_res[, c("expression", "median")])
)

res_ind <- backtest_lr(x, alpha, "ind")
res_cc  <- backtest_lr(x, alpha, "cc")

cat("\n--- Independence test ---\n")
print(res_ind)

cat("\n--- Conditional-coverage test ---\n")
print(res_cc)

## ----message=FALSE, warning=FALSE---------------------------------------------
n      <- 250
alpha  <- 0.01

d_ind <- lr_ind_dist(n, alpha)
d_cc <- lr_cc_dist(n, alpha)
d_cc$LR   <- d_cc$LR_cc    
d_cc$prob <- d_cc$prob_cc

oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1)) 

# LR_ind vs χ²(1)
curve(pchisq(x, df = 1), 0, 20, lty = 2, ylab = "CDF",
      xlab = "LR_ind statistic", main = "LR_ind")
lines(stepfun(d_ind$LR, c(0, cumsum(d_ind$prob))), pch = 19)
legend("bottomright", c("Chi-sq(1)", "Exact"), lty = c(2, 1), bty = "n")

# LR_cc vs χ²(2)
curve(pchisq(x, df = 2), 0, 30, lty = 2, ylab = "CDF",
      xlab = "LR_cc statistic", main = "LR_cc")
lines(stepfun(d_cc$LR, c(0, cumsum(d_cc$prob))), pch = 19)
legend("bottomright", c("Chi-sq(2)", "Exact"), lty = c(2, 1), bty = "n")

par(oldpar) 

## ----message=FALSE, warning=FALSE---------------------------------------------
n     <- 250
alpha <- 0.01
set.seed(1)

# LR_cc
p.chi.cc <- replicate(
  1000,
  ExactVaRTest::lr_cc_stat(rbinom(n, 1, alpha), alpha) > qchisq(.95, 2)
)
p.exact.cc <- replicate(
  1000,
  {
    x <- rbinom(n, 1, alpha)
    ExactVaRTest::backtest_lr(x, alpha, "cc")$pval < 0.05
  }
)

# LR_ind
p.chi.ind <- replicate(
  1000,
  ExactVaRTest::lr_ind_stat(rbinom(n, 1, alpha), alpha) > qchisq(.95, 1)
)
p.exact.ind <- replicate(
  1000,
  {
    x <- rbinom(n, 1, alpha)
    ExactVaRTest::backtest_lr(x, alpha, "ind")$pval < 0.05
  }
)

data.frame(
  Test = c("LR_cc Chi^2", "LR_cc Exact", "LR_ind Chi^2", "LR_ind Exact"),
  Size = c(mean(p.chi.cc), mean(p.exact.cc), mean(p.chi.ind), mean(p.exact.ind))
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

alpha <- 0.01
win   <- 250

local_file <- "inst/extdata/GSPC_close.rds"
file_path  <- if (file.exists(local_file)) local_file else
              system.file("extdata", "GSPC_close.rds", package = "ExactVaRTest")

px  <- readRDS(file_path)
px  <- px[index(px) >= "2012-01-01"]

ret <- diff(log(px))

VaR <- rollapply(
  ret, win,
  function(r) quantile(r, alpha, na.rm = TRUE),
  by.column = FALSE, align = "right"
)

keep <- complete.cases(ret, VaR)
r  <- coredata(ret[keep])
v  <- coredata(VaR[keep])
exc <- ifelse(r < v, 1L, 0L)

n_seg <- length(exc) - win + 1
ind_exact <- ind_chi <- cc_exact <- cc_chi <- numeric(n_seg)

for (i in seq_len(n_seg)) {
  seg <- exc[i:(i + win - 1)]
  ind_stat <- ExactVaRTest::lr_ind_stat(seg, alpha)
  cc_stat  <- ExactVaRTest::lr_cc_stat(seg, alpha)
  ind_exact[i] <- ExactVaRTest::pval_lr_ind(ind_stat, win, alpha)
  cc_exact[i]  <- ExactVaRTest::pval_lr_cc(cc_stat,  win, alpha)
  ind_chi[i]   <- pchisq(ind_stat, df = 1, lower.tail = FALSE)
  cc_chi[i]    <- pchisq(cc_stat,  df = 2, lower.tail = FALSE)
}

df <- tibble(
  idx = seq_len(n_seg),
  ind_exact, ind_chi, cc_exact, cc_chi
) %>%
  pivot_longer(cols = -idx,
               names_to  = c("test", "method"),
               names_pattern = "(ind|cc)_(.*)",
               values_to = "pval") %>%
  mutate(test = ifelse(test == "ind", "LRind", "LRcc"),
         method = ifelse(method == "exact", "exact", "chi-sq"))

ggplot(df, aes(idx, pval, colour = method)) +
  geom_step() +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "red") +
  facet_wrap(~ test, ncol = 1, scales = "free_x") +
  scale_colour_manual(values = c("chi-sq" = "red", "exact" = "cyan4")) +
  labs(x = "window start index", y = "p-value",
       title = "Rolling 250-day p-values (α = 1%)") +
  theme_bw()

## ----example------------------------------------------------------------------
library(ExactVaRTest)

n_set      <- c(250, 500, 750, 1000)
alpha_set  <- c(0.005, 0.01, 0.025, 0.05)
gamma_set  <- c(0.90, 0.95, 0.99)

q_lr <- function(d, g) d$LR[which(cumsum(d$prob) >= g)[1L]]

tbl <- expand.grid(n = n_set,
                   alpha = alpha_set,
                   gamma = gamma_set,
                   KEEP.OUT.ATTRS = FALSE,
                   stringsAsFactors = FALSE)

tbl$crit_ind <- mapply(function(n, a, g)
  q_lr(lr_ind_dist(n, a, prune_threshold = 1e-15), g),
  tbl$n, tbl$alpha, tbl$gamma, SIMPLIFY = TRUE)

tbl$crit_cc <- mapply(function(n, a, g) {
  d <- lr_cc_dist(n, a, prune_threshold = 1e-15)
  d$LR   <- d$LR_cc      
  d$prob <- d$prob_cc   
  q_lr(d, g)
}, tbl$n, tbl$alpha, tbl$gamma, SIMPLIFY = TRUE)

print(tbl, digits = 6)

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

