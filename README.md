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


