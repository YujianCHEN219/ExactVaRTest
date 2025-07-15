# ExactVaRTest 0.1.0 (2025-07-15)

## New features
* `backtest_lr()` — one-shot exact back-test that returns statistic,  
  exact finite-sample p-value and reject/fail decision in a single call.
* Custom S3 printer for `"ExactVaRBacktest"` objects with tidy output.

## Improvements
* Added pure-R implementations `lr_ind_stat()` and `lr_cc_stat()`;  
  wrappers automatically fall back to them when the C++ backend is absent.
* Examples now run in milliseconds for *n* ≤ 300; C++ engine remains the default.

## Internal
* First feature-complete release; bumps version to **0.1.0**.
