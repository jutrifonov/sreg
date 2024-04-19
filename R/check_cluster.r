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