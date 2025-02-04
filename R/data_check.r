check.data.types <- function(Y, S, D, G.id, Ng, X) {
  non.null.vars <- list(Y, S, D, G.id, Ng, X)[sapply(list(Y, S, D, G.id, Ng, X), Negate(is.null))]
  all.correct.types <- all(
    sapply(non.null.vars, function(var) {
      is.matrix(var) || is.numeric(var) || is.data.frame(var)
    })
  )

  if (!all.correct.types) {
    stop("Error: At least one non-NULL input variable has a different type than matrix, numeric vector, or data frame.")
  }
}
check.integers <- function(S, D, G.id, Ng) {
  non.null.vars <- list(S = S, D = D, G.id = G.id, Ng = Ng)

  for (var.name in names(non.null.vars)) {
    var <- non.null.vars[[var.name]]

    if (!is.null(var)) {
      if (is.data.frame(var)) {
        if (!all(sapply(var, function(col) all(is.na(col) | as.integer(col) == col)))) {
          stop(paste("Error: Variable", var.name, "must contain only integer values."))
        }
      } else {
        if (!all(is.na(var) | as.integer(var) == var)) {
          stop(paste("Error: Variable", var.name, "must contain only integer values."))
        }
      }
    }
  }
}
check.range <- function(var, range.min = NULL, range.max = NULL) {
  if (is.data.frame(var)) {
    for (col in names(var)) {
      data <- var[[col]]
      current_range_min <- if (is.null(range.min)) min(data, na.rm = TRUE) else range.min
      current_range_max <- if (is.null(range.max)) max(data, na.rm = TRUE) else range.max

      missing_values <- setdiff(seq(current_range_min, current_range_max), data)
      if (length(missing_values) > 0) {
        stop(paste0(
          "Error: There are skipped values in the range of ", col, ": ",
          toString(missing_values), ". ",
          "Variables S and D must not contain any skipped values within the range. For example, if min(S) = 1 and max(S) = 3, then S should encompass values 1, 2, and 3."
        ))
      }
    }
  } else {
    current_range_min <- if (is.null(range.min)) min(var, na.rm = TRUE) else range.min
    current_range_max <- if (is.null(range.max)) max(var, na.rm = TRUE) else range.max

    missing_values <- setdiff(seq(current_range_min, current_range_max), var)
    if (length(missing_values) > 0) {
      stop(paste0(
        "Error: There are skipped values in the range of ", deparse(substitute(var)), ": ",
        toString(missing_values), ". ",
        "Variables S and D must not contain any skipped values within the range. For example, if min(S) = 1 and max(S) = 3, then S should encompass values 1, 2, and 3."
      ))
    }
  }
}
boolean.check <- function(var) {
  is.boolean <- function(x) {
    is.logical(x) && length(x) == 1 && !is.na(x)
  }
  if (!is.boolean(var)) {
    stop("Error: The value of HC must be either TRUE or FALSE. A non-boolean value was provided.")
  }
}

check.within.stratatreatment.variation <- function(data) {
  covariate_columns <- names(data)[-(1:2)]

  variation_check <- data %>%
    group_by(.data$S, .data$D) %>%
    summarise(across(all_of(covariate_columns), ~ n_distinct(.) > 1, .names = "check_{.col}"),
              .groups = "drop")
  
  all_variation <- variation_check %>%
    summarise(across(starts_with("check_"), all)) %>%
    unlist()
  
  all(all_variation)
}

check.within.strata.variation <- function(data) {

  covariate_columns <- names(data)[-(1:2)]

  variation_check <- data %>%
    group_by(.data$S) %>%
    summarise(across(all_of(covariate_columns), ~ n_distinct(.) > 1, .names = "check_{.col}"),
              .groups = "drop") 
  
  all_variation <- variation_check %>%
    summarise(across(starts_with("check_"), all)) %>%
    unlist()
  
  all(all_variation)
}
