
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExactVaRTest

Fast exact finite-sample back-testing for Value-at-Risk (VaR) models in
R.

<!-- badges: start -->
<!-- badges: end -->

**ExactVaRTest** implements a forward dynamic programming algorithm (C++
backend with a pure-R fallback) that gives the exact finite-sample
distributions and p-values of Christoffersen’s (1998) independence (IND)
and conditional-coverage (CC) tests for Value-at-Risk (VaR) exception
series. In particular, it corrects the severe size distortions from
which the usual asymptotic $\chi^2$ approximation suffers in small
samples and under extreme coverage rates.

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
x <- rbinom(300, 1, 0.03)          # synthetic 0/1 exception series

bt <- backtest_lr(x, alpha = 0.05, type = "cc")  # exact LR_cc back-test
print(bt)
#> Exact finite-sample back-test
#> --------------------------------
#> Test           : Conditional-coverage (LR_cc)
#> Sample size    : 300
#> Model alpha    : 0.0500
#> Signif. level  : 0.0500
#> LR statistic   : 5.8882
#> Exact p-value  : 0.0442
#> Decision       : REJECT null at 5.00% level
```

## Main features

Exact LR distributions and p-values for any sample size $n$; for
$n \le 2\,000$ the computation finishes in milliseconds to a few
seconds.

C++ implementation (`Rcpp`) with automatic fallback to a pure-R
reference version.

`backtest_lr()` returns LR statistic, exact p-value, and reject /
fail-to-reject decision in one single call.

Minimal dependencies (`Rcpp`, `stats`), works on macOS, Linux, and
Windows.

## References

1.  Christoffersen, P. F. (1998). *Evaluating interval forecasts.*
    International economic review, 841-862.
2.  Mehta, C. R., Patel, N. R., & Gray, R. (1985). *Computing an exact
    confidence interval for the common odds ratio in several 2× 2
    contingency tables.* Journal of the American Statistical
    Association, 80(392), 969-973.

## Acknowledgements

I greatly appreciate Christian Francq, Christophe Hurlin, and
Jean-Michel Zakoian’s guidance and support.

In particular, Christian Francq generously shared the initial idea;
without his help, this package would not exist.

## License

This package is free and open source software, licensed under GPL.
