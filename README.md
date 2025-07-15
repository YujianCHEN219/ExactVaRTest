
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExactVaRTest

Fast exact finite-sample back-testing for Value-at-Risk (VaR) models in
R.

<!-- badges: start -->
<!-- badges: end -->

ExactVaRTest implements dynamic programming algorithms (C++ backend with
a pure-R fallback) that give the exact finite-sample distributions and
p-values of Christoffersen’s (1998) independence (IND) and
conditional-coverage (CC) tests for Value-at-Risk (VaR) exception
series.

A one-shot helper `backtest_lr()` returns the LR statistic, its exact
p-value, and a reject / fail-to-reject decision for a chosen
significance level.

## Installation

You can install the development version of ExactVaRTest from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("YujianCHEN219/ExactVaRTest")
```

## Example

``` r
library(ExactVaRTest)

set.seed(42)
x <- rbinom(300, 1, 0.05)          # synthetic exception series (0 = OK, 1 = VaR breach)

bt <- backtest_lr(x, alpha = 0.05, type = "cc")  # LR_cc back-test
print(bt)
#> Exact finite-sample back-test
#> --------------------------------
#> Test           : Conditional-coverage (LR_cc)
#> Sample size    : 300
#> Model alpha    : 0.0500
#> Signif. level  : 0.0500
#> LR statistic   : 1.6798
#> Exact p-value  : 0.5094
#> Decision       : fail to reject at 5.00% level
```

## Main features

Exact LR distributions and p-values for any sample size n (≈ 2 000 in
seconds).

C++ implementation (`Rcpp`) with automatic fallback to a pure-R
reference version.

High-level helper `backtest_lr()` with printing: statistic, p-value,
decision, and all parameters in one object.

MIT-licensed, minimal dependencies (`Rcpp`, `stats`).

## Acknowledgements

I greatly appreciate Christian Francq, Christophe Hurlin, and
Jean-Michel Zakoian’s guidance and support.

In particular, Christian Francq generously shared the initial idea;
without his help, this package would not exist.

## License

MIT © 2025 Yujian Chen
