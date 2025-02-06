## R CMD check results

0 errors | 0 warnings | 0 note

## Submission comments
- after the initial submission, the invalid file URL in README.md was corrected by using a fully specified URL instead.
- after the second submission, the description has been corrected
- after the second submission, all the print() were eliminated. Instead the S3 method has been added to print the results.
- after the second submission, print() in R/dgp_obs.r and R/dgp_po.r has been substituted with stop().
- after the third submission, the print.sreg() method has been exported with the description and examples provided. The problem with examples using unexported functions has been resolved.
- after the third submission, all the additional linebreaks in the description have been eliminated due to the comments received.
- after the fourth submission, the \value{No return value, called for side effects} tag has been added to the print.sreg.Rd file. 

## sreg 1.0.1
 -  Fixed a bug in the `sreg` function that caused it to return output for the unadjusted estimator instead of the adjusted estimator when `X` contained a single covariate.  
 -  Minor improvements and bug fixes.
