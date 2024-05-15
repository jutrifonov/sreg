#-------------------------------------------------------------------------
# %#     Function to check if covariates are the same within each cluster
#-------------------------------------------------------------------------
check.cluster <- function(df)
#-------------------------------------------------------------------------
{
  cov.names <- names(df)[-1]
  group.by.cov <- rlang::syms(c("G.id"))

  cov.same <- df %>%
    group_by(!!!group.by.cov) %>%
    summarize(across(all_of(cov.names), ~ n_distinct(.) == 1)) %>%
    ungroup() %>%
    summarize(all.same = all(c_across(all_of(cov.names))))

  return(cov.same$all.same)
}
check.cluster.lvl <- function(G.id, S, D, Ng) {
  dta.check.lvl <- tibble(G.id, S, D, Ng)
  # Check each cluster
  unique.clusters <- unique(dta.check.lvl$G.id)
  for (cluster in unique.clusters) {
    subset.dta <- dta.check.lvl[dta.check.lvl$G.id == cluster, ]

    suppressWarnings({
      # Check if all entries for S, D, and Ng are the same within the cluster
      if (length(unique(subset.dta$S)) > 1 ||
        length(unique(subset.dta$D)) > 1 ||
        length(unique(subset.dta$Ng)) > 1) {
        stop("Error: The values for S, D, and Ng must be consistent within each cluster (i.e., S, D, and Ng are cluster-level variables). Please verify that there are no discrepancies at the individual level within any cluster.")
      }
    })
  }
}
