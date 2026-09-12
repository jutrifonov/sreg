# Development version

* Fixed cluster-ID alignment in adjusted and unadjusted large-strata
  estimation, including the large component of mixed designs. Reordering
  observations now preserves estimates, standard errors, tests and intervals,
  with either supplied or inferred cluster sizes.

# sreg 2.1.0

## Estimation and inference

* Corrected the large-strata variance estimator for experiments with multiple
  active treatment arms under both individual- and cluster-level assignment.
  The estimator now includes the contribution of clusters or individuals
  assigned to active arms other than the arm being compared with control.
* Corrected the small-strata cluster estimator. Point estimation now uses
  expanded cluster outcomes, cluster-level covariate means, and normalization
  by the mean represented cluster size. The corresponding variance estimator
  now incorporates the common random denominator and the contributions of all
  treatment arms.
* Corrected mixed-design inference under cluster-level assignment. Component
  estimates are now weighted by their shares of the represented individual
  population, computed from `Ng`, and the variance estimator includes the
  variability of these random population shares.
* Covariates supplied for a mixed design are now used in both its small- and
  large-strata components for individual- and cluster-level assignment. A
  targeted error explains when the large-strata component cannot identify the
  requested treatment-by-stratum adjustment and recommends reducing the
  covariate set or using `X = NULL`.
* Improved HC1 handling in degenerate multi-treatment settings by reverting to
  the unadjusted variance estimate when the finite-sample correction is
  undefined.

## Design support and interface

* Added the optional `k` argument to `sreg()`. It validates the common stratum
  size in uniform small-strata designs and identifies the small-stratum size in
  general mixed designs, extending mixed-design support beyond matched pairs
  and triplets to general k-tuples.
* Improved automatic design classification, validation messages, and warnings
  for small- and mixed-strata designs under both individual- and cluster-level
  assignment.
* Standardized user-facing output to use the term "large strata" rather than
  "big strata."

## Data generation

* Extended `sreg.rgen()` to generate mixed designs through the new
  `mixed.strata` and `n.small` arguments for both individual- and cluster-level
  assignment.
* Added optional stratum-specific allocation probabilities, stratum effects,
  and treatment effects through `allocation.probs`, `stratum.effects`, and
  `treatment.effects.by.stratum` for large-strata individual-level designs.
* Clarified that `n` counts clusters when `cluster = TRUE`, strengthened input
  validation, and corrected large-strata cluster generation so that
  `is.cov = FALSE` no longer returns covariate columns.

## Documentation and maintenance

* Substantially expanded the function documentation, examples, README, and
  introductory vignette to cover large-, small-, mixed-, and cluster-randomized
  designs and the S3 print and plot methods.
* Corrected references and documented the structure of returned objects,
  cluster-size handling, cluster-level covariate aggregation, and mixed-design
  adjustment behavior.
* Expanded the automated test suite for multi-arm variance estimation,
  general k-tuple and mixed designs, cluster estimators, data generation,
  design classification, and adjustment diagnostics.

# sreg 2.0.2

# sreg 2.0.1
* CRAN release of the first stable version of sreg 2.0
# sreg 2.0.0
* Major redesign of the package to support **small strata designs** (e.g., matched pairs and n-tuples), including correct estimators under both **individual-level** and **cluster-level** treatment assignment.
* Added full support for **mixed designs** combining small and large strata, with appropriate estimators implemented.
* Introduced a new **S3 plot method** (`plot.sreg`) for visualizing estimated treatment effects and confidence intervals for objects of class `sreg`.
* Multiple bug fixes and internal improvements for stability and consistency.

# sreg 1.0.1.9000 (development version)
* Ongoing development version.

# sreg 1.0.1
* Fixed a bug in the `sreg` function that caused it to return output for the unadjusted estimator instead of the adjusted estimator when `X` contained a single covariate.  
* Minor improvements and bug fixes.

# sreg 1.0.0
* Initial CRAN release.
