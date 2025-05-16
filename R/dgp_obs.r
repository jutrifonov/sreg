#-------------------------------------------------------------------
# %#     Generates observed outcomes,
# %#     by taking as input the potential outcomes,
# %#     matrix of strata assignments, pi.vec, and
# %#     number of treatments
#----------------------------------------------------------------------
dgp.obs.sreg <- function(baseline, I.S, pi.vec, n.treat, is.cov = TRUE)
#----------------------------------------------------------------------
{
  if (n.treat != length(pi.vec)) {
    stop("The number of treatments doesn't match the length of vector pi.vec.")
  }
  num.strata <- ncol(I.S)
  n <- length(baseline$Y.0)
  A <- cbind(rep(0, n))
  l.seq <- num.strata / 2

  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  pi.matr.w <- pi.matr * pi.vec

  for (k in 1:num.strata)
  {
    index <- which(I.S[, k] == 1)
    ns <- length(index)

    A[index] <- gen.treat.sreg(pi.matr.w, ns, k)
  }

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Y.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.sreg(n.treat)
  Y.obs <- eval(parse(text = formula))

  if (is.cov == TRUE) {
    ret.list <- list(
      "Y" = Y.obs,
      "D" = A,
      "X" = baseline$X
    )
  } else {
    ret.list <- list(
      "Y" = Y.obs,
      "D" = A
    )
  }
  return(ret.list)
}
#-------------------------------------------------------------------
dgp.obs.creg <- function(baseline, I.S, pi.vec, n.treat)
#-------------------------------------------------------------------
{
  if (n.treat != length(pi.vec)) {
    stop("The number of treatments doesn't match the length of vector pi.vec.")
  }
  num.strata <- ncol(I.S)
  n <- baseline$G
  A <- cbind(rep(0, n))
  l.seq <- num.strata / 2

  pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  pi.matr.w <- pi.matr * pi.vec

  for (k in 1:num.strata)
  {
    index <- which(I.S[, k] == 1)
    ns <- length(index)

    A[index] <- gen.treat.creg(pi.matr.w, ns, k)
  }
  strata.set <- data.frame(I.S)
  strata.set$S <- max.col(strata.set)
  cluster.indicator <- baseline$cl.id
  G.seq <- seq(c(1:baseline$G))
  data.short <- data.frame(
    "cl.id" = G.seq, A, S = strata.set$S, Ng = baseline$Ng,
    baseline$X
  )
  data.long <- data.frame("cl.id" = cluster.indicator)
  merged.data <- merge(data.long, data.short, by = "cl.id")
  length(merged.data$A)
  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[5:ncol(merged.data)]
  Ng <- merged.data$Ng

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula))

  ret.list <- list(
    "Y"           = Y.obs,
    "D"           = A,
    "S"           = S,
    "Z.2"         = baseline$Z.g.2,
    "X"           = X,
    "Ng"          = Ng,
    "G.id"        = cluster.indicator,
    "cl.lvl.data" = data.short
  )
  return(ret.list)
}


dgp.obs.sreg.ss <- function(dgp_list,
                                 n.treat,
                                 k,
                                 treat_sizes,
                                 badmatch = FALSE,
                                 poutcome = FALSE) {
  # Basic checks
  stopifnot(length(treat_sizes) == (n.treat + 1))
  if (sum(treat_sizes) != k) {
    stop("Sum of treat_sizes must equal k.")
  }

  # Number of total arms is n.treat + 1, labeled 0,1,...,n.treat
  all_treat_levels <- 0:n.treat

  # -- Convert the list from dgp_po into a workable data.frame --
  #    We assume that your dgp_list has elements named "Y.0","Y.1",...,"Y.n.treat", plus W, plus X (if is.cov=TRUE).
  #    Make a data.frame that has columns for each potential outcome and your stratification variable(s).

  # First, figure out how many observations:
  n <- length(dgp_list[[paste0("Y.", 0)]]) # length of Y.0
  # Build a DF row-by-row:
  #  Potential Outcomes:
  Y_mat <- sapply(all_treat_levels, function(a) dgp_list[[paste0("Y.", a)]])
  colnames(Y_mat) <- paste0("Y.", all_treat_levels)

  # For convenience, also extract W and (optionally) X
  W <- dgp_list$W
  # If you have a data.frame X, you can splice it in; or if it's just one vector x_1, do the same.
  # Here, assume 'X' is either a single numeric or a data.frame.  We’ll just keep it as-is.
  X <- dgp_list$X

  # Combine into one data.frame:
  # We'll store W in a single column called "W".  If X is a data.frame, you may want to cbind.
  df <- data.frame(Y_mat, W = W)
  # If X is a data.frame with multiple columns, do: df <- cbind(df, X).
  # If it’s just a single vector, do:
  if (is.data.frame(X)) {
    df <- cbind(df, X)
  } else {
    # treat X as a single column
    df[["X"]] <- X
  }

  # For matching, we need a single covariate (or a score) to order or "bad‐match" by.
  # If you want to match by W, we’ll define `match_var = df[["W"]]`.
  # If you want to match by X, define `match_var = df[["X"]]`.
  # Adjust as you like.
  match_var <- df[["W"]] # or df[["X"]], etc.

  # Sort or rearrange for block creation
  if (!badmatch) {
    # Good-match approach: just sort ascending by 'match_var'
    sort_index <- order(match_var)
    df_sorted <- df[sort_index, ]
  } else {
    # "Bad matching": pair smallest with largest, etc. generalized to block size k
    # The simplest approach for pairs is to reorder: first is smallest, second is largest, third is 2nd smallest, fourth is 2nd largest, etc.
    # For block size > 2, you can do an analogous approach (though it gets trickier).
    # Here is a simple generalization that reverses half of the data row‐wise:
    o <- order(match_var) # ascending
    df_asc <- df[o, ]
    # Now we can build "bad" blocks by taking half from the front, half from the back, interleaving, etc.
    # For simplicity, we can do something like: 1st block = 1, N, 2, N-1, ...
    # The code below tries a standard “pair up extremes” approach, repeated for blocks of size k.
    # If k=2, it just does the classical pairs.  If k=3 or 4, you might want to adapt logic further.
    # This skeleton tries to keep it straightforward:
    N <- nrow(df_asc)
    # vector of indices from front/back
    half <- floor(N / 2)
    front_inds <- 1:half
    back_inds <- (N):(half + 1)
    # interleave front and back
    new_order <- as.vector(rbind(front_inds, back_inds))
    # if N is not an even multiple of k, watch out for leftover
    if (length(new_order) < N) {
      leftover <- setdiff(seq_len(N), new_order)
      new_order <- c(new_order, leftover)
    }
    df_sorted <- df_asc[new_order, ]
  }

  # Now we assign blocks of size k.  We'll loop block by block.
  N <- nrow(df_sorted)
  n.blocks <- N %/% k # integer number of full blocks
  if (n.blocks * k != N) {
    stop("The total number of units (N) is not divisible by the block size (k). Please adjust the sample size or choose a different block size for the small strata design.")
  }
  # Prepare container for assigned treatments and block labels
  A <- rep(NA, N) # the assigned treatment
  block_id <- rep(NA, N)

  # For each block, randomly permute the treat_sizes among the k units
  # so that exactly treat_sizes[1] get treat=0, treat_sizes[2] get treat=1, etc.
  # all_treat_levels = c(0,1,...,n.treat)
  for (j in seq_len(n.blocks)) {
    # indices for block j
    these_inds <- ((j - 1) * k + 1):(j * k)

    # We create a vector of length k with exactly treat_sizes[a] copies of 'a', for a in all_treat_levels
    block_treat_vector <- unlist(mapply(function(tt, tsize) rep(tt, tsize),
      tt = all_treat_levels,
      tsize = treat_sizes
    ))
    # Randomly permute that vector so we don't always assign in the same pattern
    block_treat_vector <- sample(block_treat_vector)

    # Assign them
    A[these_inds] <- block_treat_vector
    # Record the block
    block_id[these_inds] <- j
  }

  # Now figure out the observed outcome, y_i, by picking from Y.0..Y.n.treat
  # We'll do it row by row.  For row i, the assigned treatment is A[i].  So the outcome is df_sorted[i, Y.(A[i])]
  # We can do that with indexing or by a small loop:
  y_obs <- numeric(N)
  # column names for the potential outcomes:
  pot_outcome_names <- paste0("Y.", all_treat_levels)
  for (i in seq_len(N)) {
    treat_i <- A[i]
    # which column?  "Y.0" => column 1, "Y.1" => column 2, etc.
    col_idx <- which(all_treat_levels == treat_i)
    y_obs[i] <- df_sorted[i, pot_outcome_names[col_idx]]
  }

  # Final data frame
  # Minimal columns: (y, A, W, X, block)
  out_df <- data.frame(
    y = y_obs,
    A = A,
    W = df_sorted[["W"]],
    block = block_id
  )
  # If your df_sorted includes multiple X columns (like x_1, x_2, etc.), you can add them:
  out_df <- cbind(out_df, df_sorted[, c("x_1", "x_2")])
  # if ("x_1" %in% names(df_sorted)) {
  #  out_df[["X"]] <- df_sorted[["x_1"]] # NB CORRECT FOR MULTIPLE COVARIATES!!! Working on it... No flexibility needed for dgp functions!
  # }
  # Optionally attach all potential outcomes for debugging/evaluation
  if (poutcome) {
    out_df <- cbind(
      out_df,
      df_sorted[, pot_outcome_names, drop = FALSE]
    )
  }

  # Return
  rownames(out_df) <- NULL
  return(out_df)
}

#-------------------------------------------------------------------
dgp.obs.creg.ss <- function(baseline, n.treat, k, treat_sizes)
#-------------------------------------------------------------------
{
  # num.strata <- ncol(I.S)
  n <- baseline$G
  A <- cbind(rep(0, n))
  # Number of total arms is n.treat + 1, labeled 0,1,...,n.treat
  all_treat_levels <- 0:n.treat
  # l.seq <- num.strata / 2

  # pi.matr <- matrix(1, ncol = num.strata, nrow = n.treat)
  # pi.matr.w <- pi.matr * pi.vec

  # for (k in 1:num.strata)
  # {
  #  index <- which(I.S[, k] == 1)
  #  ns <- length(index)

  #   A[index] <- gen.treat.creg(pi.matr.w, ns, k)
  # }
  Y_mat <- sapply(all_treat_levels, function(a) baseline[[paste0("Yig.", a)]])
  colnames(Y_mat) <- paste0("Y.", all_treat_levels)

  cluster.indicator <- baseline$cl.id
  W <- baseline$Z.g.2
  X <- baseline$X

  G.seq <- seq(c(1:baseline$G))

  # Just on a cluster level cl.id, Ng, W, X
  data.short <- data.frame(
    "cl.id" = G.seq, Ng = baseline$Ng,
    W, X
  )

  # print(data.short)

  # Stratification. We use the cluster level data to first assign strata
  match_var <- data.short[["W"]] # or df[["X"]], etc.

  sort_index <- order(match_var)
  data.short_sorted <- data.short[sort_index, ]

  # Now we assign blocks of size k.  We'll loop block by block.
  N <- nrow(data.short_sorted)
  n.blocks <- N %/% k # integer number of full blocks
  if (n.blocks * k != N) {
    stop("The total number of units (N) is not divisible by the block size (k). Please adjust the sample size or choose a different block size for the small strata design.")
  }

  # Prepare container for assigned treatments and block labels
  A <- rep(NA, N) # the assigned treatment
  block_id <- rep(NA, N)

  # For each block, randomly permute the treat_sizes among the k units
  # so that exactly treat_sizes[1] get treat=0, treat_sizes[2] get treat=1, etc.
  # all_treat_levels = c(0,1,...,n.treat)
  for (j in seq_len(n.blocks)) {
    # indices for block j
    these_inds <- ((j - 1) * k + 1):(j * k)

    # We create a vector of length k with exactly treat_sizes[a] copies of 'a', for a in all_treat_levels
    block_treat_vector <- unlist(mapply(function(tt, tsize) rep(tt, tsize),
      tt = all_treat_levels,
      tsize = treat_sizes
    ))
    # Randomly permute that vector so we don't always assign in the same pattern
    block_treat_vector <- sample(block_treat_vector)

    # Assign them
    A[these_inds] <- block_treat_vector
    # Record the block
    block_id[these_inds] <- j
  }

  data.short <- data.frame(
    A,
    S = block_id, data.short_sorted
  )
  # print(data.short)
  # print(block_id)
  # print(A)

  # print(data.short_sorted)
  # short data
  # data.short <- data.frame(
  #  "cl.id" = G.seq, A, S = strata.set$S, Ng = baseline$Ng,
  #  baseline$X
  # )

  data.long <- data.frame("cl.id" = cluster.indicator)
  # print(data.long)
  # individual level data frame without Y! Just cl.id, Ng, W, X
  merged.data <- merge(data.long, data.short, by = "cl.id")
  # print(head(merged.data, 200)) # Works! But we need to append potential outcomes!

  A <- merged.data$A
  S <- merged.data$S
  X <- merged.data[6:ncol(merged.data)]
  Ng <- merged.data$Ng

  for (a in 0:n.treat)
  {
    assign(paste("Y.", a, sep = ""), baseline[[paste("Yig.", a, sep = "")]])
  }
  # print(Y.1)
  formula <- gen.rubin.formula.creg(n.treat)
  Y.obs <- eval(parse(text = formula)) # observed outcomes for every individual
  # print(Y.obs)

  ret.list <- list(
    "Y"           = Y.obs,
    "D"           = A,
    "S"           = S,
    "W"           = baseline$Z.g.2,
    "X"           = X,
    "Ng"          = Ng,
    "G.id"        = cluster.indicator,
    "cl.lvl.data" = data.short
  )
  return(ret.list)
}
