# sreg 2.0.2

# sreg 2.0.1
* CRAN release of the first stable version of sreg 2.0
# sreg 2.0.0
* Major redesign of the package to support **small strata designs** (e.g., matched pairs and n-tuples), including correct estimators under both **individual-level** and **cluster-level** treatment assignment.
* Added full support for **mixed designs** combining small and big strata, with appropriate estimators implemented.
* Introduced a new **S3 plot method** (`plot.sreg`) for visualizing estimated treatment effects and confidence intervals for objects of class `sreg`.
* Multiple bug fixes and internal improvements for stability and consistency.

# sreg 1.0.1.9000 (development version)
* Ongoing development version.

# sreg 1.0.1
* Fixed a bug in the `sreg` function that caused it to return output for the unadjusted estimator instead of the adjusted estimator when `X` contained a single covariate.  
* Minor improvements and bug fixes.

# sreg 1.0.0
* Initial CRAN release.
