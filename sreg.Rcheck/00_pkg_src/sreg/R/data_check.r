check.data.types <- function(Y, S, D, G.id, Ng, X) {
  # Create a list of non-NULL variables
  non.null.vars <- list(Y, S, D, G.id, Ng, X)[sapply(list(Y, S, D, G.id, Ng, X), Negate(is.null))]

  # Check types only for non-NULL variables
  all.correct.types <- all(
    sapply(non.null.vars, function(var) {
      is.matrix(var) || is.numeric(var) || is.data.frame(var)
    })
  )

  if (!all.correct.types) {
    stop("At least one non-NULL input variable has a different type than matrix, numeric vector, or data frame.")
  }
}
check.integers <- function(S, D, G.id, Ng) {
  non.null.vars <- list(S = S, D = D, G.id = G.id, Ng = Ng)

  for (var.name in names(non.null.vars)) {
    if (!is.null(non.null.vars[[var.name]]) &&
      !all(is.na(non.null.vars[[var.name]]) |
        as.integer(non.null.vars[[var.name]]) == non.null.vars[[var.name]])) {
      stop(paste("Variable", var.name, "must contain only integer values."))
    }
  }
}
