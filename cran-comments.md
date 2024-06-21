## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Submission comments
- after the initial submission, the invalid file URL in README.md was corrected by using a fully specified URL instead.
- after the second submission, the description has corrected
- after the second submission, all the print() were eliminated. Instead the S3 method has been added to print the results.
- after the second submission, print() in R/dgp_obs.r and R/dgp_po.r has been substituted with stop().
- after the third submission, the print.sreg() method has been exported with the description and examples provided. The problem with examples using unexported functions has been resolved.
- after the third submission, all the additional linebreaks in the description have been eliminated due to the comments received.