% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{backtest_lr}
\alias{backtest_lr}
\title{Exact finite‑sample back‑test for a VaR exception series}
\usage{
backtest_lr(
  x,
  alpha = 0.05,
  type = c("uc", "ind", "cc"),
  sig = 0.05,
  prune_threshold = 1e-15
)
}
\arguments{
\item{x}{0/1 exception series.}

\item{alpha}{Exception probability.}

\item{type}{\code{"uc"}, \code{"ind"} or \code{"cc"}.}

\item{sig}{Significance level (default \code{0.05}).}

\item{prune_threshold}{Passed to the dynamic‑programming engine.}
}
\value{
An object of class \code{"ExactVaRBacktest"}.
}
\description{
Exact finite‑sample back‑test for a VaR exception series
}
\examples{
set.seed(123)
x <- rbinom(250, 1, 0.01)
backtest_lr(x, alpha = 0.01, type = "uc")
}
