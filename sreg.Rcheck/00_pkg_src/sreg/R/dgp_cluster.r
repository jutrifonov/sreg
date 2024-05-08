#-------------------------------------------------------------------
# %#   Cluster sizes generation
#-------------------------------------------------------------------
gen.cluster.sizes <- function(G, max.support)
#-------------------------------------------------------------------
{
  sample <- 10 * (rbbinom(G, max.support, alpha = 1, beta = 1) + 1)
  return(sample)
}
