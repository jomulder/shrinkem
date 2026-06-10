# shrinkem 0.3.0
* Date: 2026-06-10

## Submission comment
This is a minor update. The internal Gibbs samplers have been
reimplemented in C++ (Rcpp/RcppArmadillo) for performance; the
user-facing interface and output are unchanged.

## Test environments
* Local OS X Tahoe 26.2, R 4.5.3, clang-1400.0.29.202, GNU Fortran 12.2.0
* rhub check: Windows Server 2022, R-devel, 64 bit
* rhub check: Windows Server 2022, R version 4.4.0
* rhub check: ubuntu-latest, R-devel

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs

