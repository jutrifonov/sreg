## Release summary

This is a substantial update from CRAN version 2.0.2. Version 2.1.0:

* corrects multi-arm large-strata variance estimation under individual- and
  cluster-level assignment;
* corrects point and variance estimation for small-strata cluster designs;
* corrects the weighting and variance formula for mixed cluster designs;
* applies supplied covariates to both components of individual- and
  cluster-level mixed designs, with an informative identification check;
* extends mixed-design estimation to general k-tuples through the optional
  `k` argument;
* extends `sreg.rgen()` to mixed designs and to customizable large-strata
  individual-level designs; and
* substantially expands the documentation, examples, vignette, and automated
  tests.

The corrected estimators can change numerical results for the affected
multi-arm and cluster-randomized designs. Existing function arguments retain
their previous meanings; the new arguments are optional.

## Test environments

* Local: macOS 26.3 (arm64), R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

The note is local and concerns optional tools used to validate the HTML
manual: the installed HTML Tidy is not recent enough, and the optional `V8`
package is unavailable. Consequently, HTML validation and math-rendering
checks were skipped. The PDF manual, examples, tests, and vignettes all pass.

## Reverse dependencies

There are currently no reverse dependencies on CRAN.
