![Static Badge](https://img.shields.io/badge/sreg-0.5.8(dev)-green?logo=GitHub) ![Static Badge](https://img.shields.io/badge/CRAN-Coming%20Soon!-orange?logo=R)

# Stratified Randomized Experiments
This repository houses the 'sreg' package for R, offering a toolkit for estimating average treatment effects (ATEs) in stratified randomized experiments. The package is designed to accommodate scenarios with multiple treatments and cluster-level treatment assignments, and accomodates optimal linear covariate adjustment based on baseline observable characteristics. The package computes estimators and standard errors based on Bugni, Canay, Shaikh (2018), Bugni, Canay, Shaikh, Tabord-Meehan (2023), and Jiang, Linton, Tang, Zhang (2023).

## Installation
The latest version can be installed using `devtools`. The official `CRAN` release will be available soon.
````{r}
library(devtools)
install_github("yurytrifonov/sreg")
````
```{r}
Downloading GitHub repo yurytrifonov/sreg@HEAD
── R CMD build ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/private/var/folders/mp/06gjwr8j56zdp5j2vgdkd4z40000gq/T/RtmpZh7j1Y/remotesfbf765906644/yurytrifonov-sreg-91d11dc/DESCRIPTION’ ...
─  preparing ‘sreg’:
✔  checking DESCRIPTION meta-information
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘sreg_0.5.8.tar.gz’
   
* installing *source* package ‘sreg’ ...
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (sreg)
```
````{r}
library(sreg)
````
```{r}
  ____  ____  _____ ____      Stratified Randomized
 / ___||  _ \| ____/ ___|     Experiments
 \___ \| |_) |  _|| |  _  
  ___) |  _ <| |__| |_| |  
 |____/|_| \_\_____\____| version 0.5.8
                           
```

## The function `sreg()`
Estimates the ATE(s) and the corresponding standard error(s) for a (collection of) treatment(s) relative to a control.
### Syntax
````{r}
sreg(Y, S = NULL, D, G.id = NULL, Ng = NULL, X = NULL, HC1 = TRUE)
````
### Arguments
- **`Y` -** a numeric vector/matrix/data frame of the observed outcomes;
- **`S` -** a numeric vector/matrix/data frame of strata indicators; if `NULL` then the estimation is performed assuming no stratification;
- **`D` -** a numeric  vector/matrix/data frame of treatments indexed by $\\{0, 1, 2, \ldots\\}$, where `D = 0` denotes the control;
- **`G.id` -** a numeric vector/matrix/data frame of cluster indicators; if `NULL` then estimation is performed assuming treatment is assigned at the individual level;
- **`Ng` -** a numeric vector/matrix/data frame of cluster sizes; if `NULL` then `Ng` is assumed to be equal to the number of available observations in every cluster;
- **`X` -** a data frame with columns representing the covariate values for every observation; if `NULL` then the estimator without linear adjustments is applied [^*];
- **`HC1` -** a `TRUE/FALSE` logical argument indicating whether the small sample correction should be applied to the variance estimator.
[^*]: *Note: sreg cannot use individual-level covariates for covariate adjustment in cluster-randomized experiments. Any individual-level covariates will be aggregated to their cluster-level averages.*

### Data Structure
Here we provide an example of a data frame that can be used with `sreg`.
````{r}
|       Y      | S | D | G.id | Ng |     x_1    |      x_2      |
|--------------|---|---|------|----|------------|---------------|
| -0.57773576  | 2 | 0 |  1   | 10 |  1.5597899 |  0.03023334   |
|  1.69495638  | 2 | 0 |  1   | 10 |  1.5597899 |  0.03023334   |
|  2.02033740  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
|  1.22020493  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
|  1.64466086  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
| -0.32365109  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
|  2.21008191  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
| -2.25064316  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
|  0.37962312  | 4 | 2 |  2   | 30 |  0.8747419 | -0.77090031   |
````
### Value
#### Summary
`sreg` prints a *"Stata-style"* table containing the ATE estimates, corresponding standard errors, $t$-statistics, $p$-values, $95\\%$ asymptotic confidence intervals, and significance indicators for different levels $\alpha$. The example of the printed output is provided below.
```{r}
Saturated Model Estimation Results under CAR with clusters and linear adjustments
Observations: 30000 
Clusters: 1000 
Number of treatments: 2 
Number of strata: 4 
Covariates used in linear adjustments: x_1, x_2
---
Coefficients:
      Tau   As.se   T-stat P-value CI.left(95%) CI.right(95%) Significance
1 0.01614 0.04513  0.35753 0.72069     -0.07232        0.1046             
2 0.78642 0.04642 16.94263 0.00000      0.69545        0.8774          ***
---
Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1 ` ` 1
```
#### Return Value
The function returns an object of class `sreg` that is a list containing the following elements:
- **`tau.hat` -**  a $1 \times |\mathcal A|$ vector of ATE estimates, where $|\mathcal A|$ represents the number of treatments;
- **`se.rob` -** a $1 \times |\mathcal A|$ vector of standard errors estimates, where $|\mathcal A|$ represents the number of treatments;
- **`t.stat` -** a $1 \times |\mathcal A|$ vector of $t$-statistics, where $|\mathcal A|$ represents the number of treatments;
- **`p.value` -** a $1 \times |\mathcal A|$ vector of corresponding $p$-values, where $|\mathcal A|$ represents the number of treatments;
- **`CI.left` -** a $1 \times |\mathcal A|$ vector of the left bounds of the $95\\%$ as. confidence interval;
- **`CI.right` -** a $1 \times |\mathcal A|$ vector of the right bounds of the $95\\%$ as. confidence interval;
- **`data` -** an original data of the form `data.frame(Y, S, D, G.id, Ng, X)`;
- **`lin.adj` -** a data frame representing the covariates that were used in implementing linear adjustments.




