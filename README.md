![Static Badge](https://img.shields.io/badge/sreg%20-%200.5.8%20(dev)%20-%20orrange?style=plastic)
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
### Data Example
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



