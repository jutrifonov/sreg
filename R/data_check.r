# Suppress NSE notes during R CMD check
utils::globalVariables(c("size", "freq"))

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
boolean.check.ss <- function(var) {
  is.boolean <- function(x) {
    is.logical(x) && length(x) == 1 && !is.na(x)
  }
  if (!is.boolean(var)) {
    stop("Error: The value of small.strata must be either TRUE or FALSE. A non-boolean value was provided.")
  }
}

check.within.stratatreatment.variation <- function(data) {
  # covariate_columns <- names(data)[-(1:2)]

  # Ensure the columns exist
  stopifnot(all(c("S", "D") %in% names(data)))

  # Exclude known non-covariate columns
  covariate_columns <- setdiff(names(data), c("S", "D", "G.id"))


  variation_check <- data %>%
    group_by(.data$S, .data$D) %>%
    summarise(across(all_of(covariate_columns), ~ n_distinct(.) > 1, .names = "check_{.col}"),
      .groups = "drop"
    )

  all_variation <- variation_check %>%
    summarise(across(starts_with("check_"), all)) %>%
    unlist()

  all(all_variation)
}

check.within.strata.variation <- function(data) {
  # Ensure the columns exist
  stopifnot(all(c("S", "D") %in% names(data)))

  # Exclude known non-covariate columns
  covariate_columns <- setdiff(names(data), c("S", "D", "G.id"))

  variation_check <- data %>%
    group_by(.data$S) %>%
    summarise(across(all_of(covariate_columns), ~ n_distinct(.) > 1, .names = "check_{.col}"),
      .groups = "drop"
    )

  all_variation <- variation_check %>%
    summarise(across(starts_with("check_"), all)) %>%
    unlist()

  all(all_variation)
}

# ------------------------------------------------------------------
#  Classify strata as "small" or "big" in a mixed design
# ------------------------------------------------------------------
# design.classifier <- function(data, S, G.id = NULL, keep.size = FALSE, warn = TRUE) {
#   ## 1. Capture column names (works with bare names or strings)
#   S_name <- if (is.character(substitute(S))) substitute(S) else deparse(substitute(S))
#   G_name <- if (!missing(G.id)) {
#     if (is.character(substitute(G.id))) substitute(G.id) else deparse(substitute(G.id))
#   } else {
#     NULL
#   }

#   ## 2. Compute stratum sizes
#   if (!is.null(G_name)) {
#     # -- clustered: count *clusters* per stratum
#     cluster_strata <- dplyr::distinct(data, .data[[S_name]], .data[[G_name]])
#     strata_sizes   <- dplyr::count(cluster_strata, .data[[S_name]], name = "size")
#   } else {
#     # -- individual: count *units* per stratum
#     strata_sizes   <- dplyr::count(data, .data[[S_name]], name = "size")
#   }

#   ## 3. If every stratum has the same size, nothing to classify
#   if (length(unique(strata_sizes$size)) == 1) {
#     if (!keep.size) strata_sizes$size <- NULL  # cosmetic
#     return(dplyr::left_join(data, strata_sizes, by = S_name))
#   }

#   ## 4. Find the modal size and classify
#   modal_size <- strata_sizes %>%
#     dplyr::count(size) %>%
#     dplyr::arrange(dplyr::desc(n)) %>%
#     dplyr::pull(size) %>%
#     .[[1]]

#   strata_sizes <- strata_sizes %>%
#     dplyr::mutate(stratum_type = ifelse(size == modal_size, "small", "big"))

#   if (!keep.size) strata_sizes$size <- NULL

#   ## 5. Merge back and optionally warn
#   out <- dplyr::left_join(data, strata_sizes, by = S_name)

#   if (warn)
#     warning("Mixed design detected: strata of varying cluster/unit counts. ",
#             "Weighted estimators will be used.", call. = FALSE)

#   return(out)
# }

design.classifier <- function(data, S, G.id = NULL, keep.size = FALSE, warn = TRUE, small.strata = TRUE) {
  S_name <- if (is.character(substitute(S))) substitute(S) else deparse(substitute(S))
  G_name <- if (!missing(G.id)) {
    if (is.character(substitute(G.id))) substitute(G.id) else deparse(substitute(G.id))
  } else {
    NULL
  }

  # ------------------------------------------------------------------
  # If the user passed duplicate column names (e.g., two "Ng" columns),
  # keep the last occurrence and drop the rest.  This preserves the
  # X columns (last ones) and ensures that downstream code referencing
  # X_names and "Ng" will work, while preventing dplyr joins
  # from failing due to non‑unique names.
  # ------------------------------------------------------------------
  # Keep the *last* occurrence of any duplicated column name (so X columns survive)
  data <- data[, !duplicated(names(data), fromLast = TRUE), drop = FALSE]

  if (!small.strata) {
    # Always compute strata sizes for warning purposes
    if (!is.null(G_name)) {
      cluster_strata <- dplyr::distinct(data, .data[[S_name]], .data[[G_name]])
      strata_sizes <- dplyr::count(cluster_strata, .data[[S_name]], name = "size")
    } else {
      strata_sizes <- dplyr::count(data, .data[[S_name]], name = "size")
    }

    unique_sizes <- unique(strata_sizes$size)
    n_strata <- nrow(strata_sizes)

    # Count frequencies of each size
    size_counts <- strata_sizes %>%
      count(size, name = "count") %>%
      mutate(freq = count / n_strata)

    # Check for small strata that meet the 25% rule
    small_modal_sizes <- size_counts %>%
      filter(size <= 3, freq >= 0.25) %>%
      arrange(desc(count))

    if (warn && nrow(small_modal_sizes) > 0) {
      warning("At least 25% of strata are small, but small.strata = FALSE. If the experimental design includes small strata (e.g., matched pairs or triplets), then setting small.strata = FALSE may lead to invalid standard errors. If strata sizes vary, the design may be mixed. In that case, setting small.strata = TRUE will apply estimators suitable for such mixed designs.", call. = FALSE)
    }

    return(data)
  }

  if (!is.null(G_name)) {
    cluster_strata <- dplyr::distinct(data, .data[[S_name]], .data[[G_name]])
    strata_sizes <- dplyr::count(cluster_strata, .data[[S_name]], name = "size")
  } else {
    strata_sizes <- dplyr::count(data, .data[[S_name]], name = "size")
  }

  unique_sizes <- unique(strata_sizes$size)
  n_strata <- nrow(strata_sizes)

  if (length(unique_sizes) == 1) {
    strata_sizes$stratum_type <- "small"
    if (!keep.size) strata_sizes$size <- NULL
    # Remove duplicate names in the join data (for safety)
    data_renamed <- data
    data_renamed <- data_renamed[, !duplicated(names(data_renamed), fromLast = TRUE), drop = FALSE]
    return(dplyr::left_join(data_renamed, strata_sizes, by = S_name))
  }

  # Count frequencies of each size
  size_counts <- strata_sizes %>%
    count(size, name = "count") %>%
    mutate(freq = count / n_strata)

  # Filter for small strata that meet the 25% rule
  small_modal_sizes <- size_counts %>%
    filter(size <= 3, freq >= 0.25) %>%
    arrange(desc(count))

  if (nrow(small_modal_sizes) == 0) {
    stop("Invalid input: Either all strata are large or too few strata qualify as 'small' to proceed with small.strata = TRUE. Please set small.strata = FALSE.

Note: This error may also occur if G.id is NULL despite the data being cluster-randomized. If your design involves clusters, make sure to specify G.id correctly.")
  }

  modal_size <- small_modal_sizes$size[1]

  # Classify
  strata_sizes <- strata_sizes %>%
    mutate(stratum_type = ifelse(size == modal_size, "small", "big"))

  if (!keep.size) strata_sizes$size <- NULL
  # Remove duplicate names in the join data (for safety)
  data_renamed <- data
  data_renamed <- data_renamed[, !duplicated(names(data_renamed), fromLast = TRUE), drop = FALSE]
  out <- dplyr::left_join(data_renamed, strata_sizes, by = S_name)

  if (warn && any(strata_sizes$stratum_type == "big")) {
    warning("Mixed design detected: at least 25% of strata are small. Weighted estimators will be used.", call. = FALSE)
  }

  return(out)
}
