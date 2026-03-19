#' Check Register Data Consistency
#'
#' Apply a set of consistency checks to person-year register data and return
#' the cleaned data alongside the IDs removed by each check.
#'
#' @param data A data.frame in long (person-year) format.
#' @param id_col Column name for individual identifiers.
#' @param year_col Column name for calendar year.
#' @param firstimmig_col Column name for first immigration year.
#' @param death_col Column name for death indicator (0/1).
#' @param emig_col Column name for emigration indicator (0/1).
#' @param immig_col Column name for immigration indicator (0/1).
#' @param reimmig_col Column name for re-immigration indicator (0/1).
#' @param year_beginning First year in the study window.
#' @param final_year Last year in the study window.
#'
#' @return A list with components:
#' \describe{
#'   \item{data}{Filtered data.frame after removing inconsistent IDs.}
#'   \item{removed}{Named list of ID vectors removed at each check.}
#' }
#' @export
oc2_check_register_data <- function(
  data,
  id_col = "id",
  year_col = "year",
  firstimmig_col = "firstimmig",
  death_col = "death",
  emig_col = "emig",
  immig_col = "immig",
  reimmig_col = "reimmig",
  year_beginning,
  final_year
) {
  stopifnot(is.data.frame(data))
  cols <- c(id_col, year_col, firstimmig_col, death_col, emig_col, immig_col, reimmig_col)
  missing_cols <- cols[!cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (missing(year_beginning) || missing(final_year)) {
    stop("year_beginning and final_year must be provided.")
  }

  id <- data[[id_col]]
  year <- data[[year_col]]
  firstimmig <- data[[firstimmig_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  data <- data[order(id, year), , drop = FALSE]

  removed <- list()

  # Check 1: first immigration after first observed year
  min_year <- tapply(year, id, min, na.rm = TRUE)
  first_by_id <- tapply(firstimmig, id, function(x) x[1])
  check1_ids <- names(min_year)[(min_year - first_by_id) > 0]
  if (length(check1_ids) > 0) {
    removed$check1 <- check1_ids
    data <- data[!(id %in% check1_ids), , drop = FALSE]
  }

  # Check 2: single record with death or emigration
  counts <- tapply(id, id, length)
  death_sum <- tapply(death, id, sum, na.rm = TRUE)
  emig_sum <- tapply(emig, id, sum, na.rm = TRUE)
  check2_ids <- names(counts)[counts == 1 & (death_sum == 1 | emig_sum == 1)]
  if (length(check2_ids) > 0) {
    removed$check2 <- check2_ids
    data <- data[!(id %in% check2_ids), , drop = FALSE]
  }

  id <- data[[id_col]]
  year <- data[[year_col]]
  firstimmig <- data[[firstimmig_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  # Check 3: disappear before final_year without death/emig
  max_year <- tapply(year, id, max, na.rm = TRUE)
  death_sum <- tapply(death, id, sum, na.rm = TRUE)
  emig_sum <- tapply(emig, id, sum, na.rm = TRUE)
  check3_ids <- names(max_year)[death_sum == 0 & emig_sum == 0 & max_year < final_year]
  if (length(check3_ids) > 0) {
    removed$check3 <- check3_ids
    data <- data[!(id %in% check3_ids), , drop = FALSE]
  }

  # Check 4: records after death year
  death_year <- tapply(ifelse(death == 1, year, NA), id, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else min(x)
  })
  data$oc2_data_end <- ifelse(!is.na(death_year[as.character(id)]),
                              death_year[as.character(id)],
                              final_year)
  data <- data[year <= data$oc2_data_end, , drop = FALSE]
  data$oc2_data_end <- NULL

  id <- data[[id_col]]
  year <- data[[year_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  # Check 5: gap years for those without reimmigration
  reimmig_sum <- tapply(reimmig, id, sum, na.rm = TRUE)
  no_reimmig_ids <- names(reimmig_sum)[reimmig_sum == 0]
  check5_ids <- character(0)
  for (cid in no_reimmig_ids) {
    years <- year[id == cid]
    if (length(years) > 1) {
      if (mean(diff(sort(years)), na.rm = TRUE) != 1) {
        check5_ids <- c(check5_ids, cid)
      }
    }
  }
  if (length(check5_ids) > 0) {
    removed$check5 <- unique(check5_ids)
    data <- data[!(id %in% check5_ids), , drop = FALSE]
  }

  id <- data[[id_col]]
  year <- data[[year_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  # Check 6: leave more than once without return (immig - emig not 0 or 1)
  immig_sum <- tapply(immig, id, sum, na.rm = TRUE)
  emig_sum <- tapply(emig, id, sum, na.rm = TRUE)
  dif_mig <- immig_sum - emig_sum
  check6_ids <- names(dif_mig)[!(dif_mig %in% c(0, 1))]
  if (length(check6_ids) > 0) {
    removed$check6 <- check6_ids
    data <- data[!(id %in% check6_ids), , drop = FALSE]
  }

  id <- data[[id_col]]
  year <- data[[year_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  # Check 7: leave and return in same year
  check7_ids <- unique(id[reimmig == 1 & emig == 1])
  if (length(check7_ids) > 0) {
    removed$check7 <- check7_ids
    data <- data[!(id %in% check7_ids), , drop = FALSE]
  }

  # Check 8: immig and emig same year
  check8_ids <- unique(id[immig == 1 & emig == 1])
  if (length(check8_ids) > 0) {
    removed$check8 <- check8_ids
    data <- data[!(id %in% check8_ids), , drop = FALSE]
  }

  id <- data[[id_col]]
  year <- data[[year_col]]
  death <- data[[death_col]]
  emig <- data[[emig_col]]
  immig <- data[[immig_col]]
  reimmig <- data[[reimmig_col]]

  # Check 9: activity between emig and reimmig
  leave_year <- tapply(ifelse(emig == 1, year, NA), id, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else min(x)
  })
  return_year <- tapply(ifelse(reimmig == 1, year, NA), id, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else min(x)
  })
  check9_ids <- character(0)
  for (cid in unique(id)) {
    ly <- leave_year[as.character(cid)]
    ry <- return_year[as.character(cid)]
    if (!is.na(ly) && !is.na(ry)) {
      if (any(year[id == cid] > ly & year[id == cid] < ry)) {
        check9_ids <- c(check9_ids, cid)
      }
    }
  }
  if (length(check9_ids) > 0) {
    removed$check9 <- unique(check9_ids)
    data <- data[!(id %in% check9_ids), , drop = FALSE]
  }

  # Check 10: immig and death same year
  check10_ids <- unique(id[immig == 1 & death == 1])
  if (length(check10_ids) > 0) {
    removed$check10 <- check10_ids
    data <- data[!(id %in% check10_ids), , drop = FALSE]
  }

  # Check 11: re-register after last emig without immig
  emig_year <- ifelse(emig == 1, year, 0)
  immig_year <- ifelse(immig == 1, year, 0)
  max_year <- tapply(year, id, max, na.rm = TRUE)
  max_emig <- tapply(emig_year, id, max, na.rm = TRUE)
  max_immig <- tapply(immig_year, id, max, na.rm = TRUE)
  check11_ids <- names(max_year)[max_emig > max_immig & max_emig < max_year]
  if (length(check11_ids) > 0) {
    removed$check11 <- check11_ids
    data <- data[!(id %in% check11_ids), , drop = FALSE]
  }

  list(data = data, removed = removed)
}

#' Prepare Register Variables for Modeling
#'
#' Standardize and clean register indicators (e.g., convert to binary,
#' replace missing values, and zero-out registers on emigration/death).
#'
#' @param data A data.frame in long (person-year) format.
#' @param register_cols Character vector of register indicator columns.
#' @param emig_col Column name for emigration indicator (0/1).
#' @param death_col Column name for death indicator (0/1).
#' @param binary_rules Optional named list describing binary conversions.
#'   Each entry should include `source` and optionally `threshold`.
#' @param na_to_zero Logical; if TRUE, replace NA with 0 in register columns.
#' @param zero_on_emig Logical; if TRUE, set register columns to 0 on emigration.
#' @param zero_on_death Logical; if TRUE, set register columns to 0 on death.
#'
#' @return A data.frame with cleaned register columns.
#' @export
oc2_prepare_register_data <- function(
  data,
  register_cols,
  emig_col = "emig",
  death_col = "death",
  binary_rules = NULL,
  na_to_zero = TRUE,
  zero_on_emig = TRUE,
  zero_on_death = TRUE
) {
  stopifnot(is.data.frame(data))
  if (missing(register_cols) || length(register_cols) == 0) {
    stop("register_cols must be provided.")
  }
  missing_cols <- register_cols[!register_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing register columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(binary_rules)) {
    for (new_name in names(binary_rules)) {
      rule <- binary_rules[[new_name]]
      if (is.null(rule$source)) {
        stop("binary_rules entries must include a 'source' field.")
      }
      if (!rule$source %in% names(data)) {
        stop("binary_rules source column not found: ", rule$source)
      }
      threshold <- if (!is.null(rule$threshold)) rule$threshold else 0
      data[[new_name]] <- ifelse(data[[rule$source]] > threshold, 1, 0)
      if (!new_name %in% register_cols) {
        register_cols <- c(register_cols, new_name)
      }
    }
  }

  if (na_to_zero) {
    for (col in register_cols) {
      data[[col]][is.na(data[[col]])] <- 0
    }
    if (emig_col %in% names(data)) {
      data[[emig_col]][is.na(data[[emig_col]])] <- 0
    }
    if (death_col %in% names(data)) {
      data[[death_col]][is.na(data[[death_col]])] <- 0
    }
  }

  if (zero_on_emig && emig_col %in% names(data)) {
    idx <- which(data[[emig_col]] == 1)
    if (length(idx) > 0) {
      data[idx, register_cols] <- 0
    }
  }

  if (zero_on_death && death_col %in% names(data)) {
    idx <- which(data[[death_col]] == 1)
    if (length(idx) > 0) {
      data[idx, register_cols] <- 0
    }
  }

  data
}

#' Build Observation Inputs for model_BLB
#'
#' Create observation matrices and combinations for the BLB model from
#' prepared register data. This helper is intended for internal use.
#'
#' @param data A data.frame in long (person-year) format.
#' @param id_col Column name for individual identifiers.
#' @param year_col Column name for calendar year.
#' @param register_cols Character vector of register indicator columns.
#' @param covariate_cols Character vector of covariate columns to include in
#'   observation combinations.
#' @param reimmig_col Column name for re-immigration indicator (0/1).
#' @param death_col Column name for death indicator (0/1).
#' @param emig_col Column name for emigration indicator (0/1).
#' @param year_beginning First year in the study window.
#' @param final_year Last year in the study window.
#' @param combins Optional matrix of combination indices used for unobserved years.
#'
#' @return A list with elements `y_matrix`, `X`, and `num_combos`.
#' @keywords internal
oc2_prepare_model_BLB_inputs <- function(
  data,
  id_col = "id",
  year_col = "year",
  register_cols,
  covariate_cols,
  reimmig_col = "reimmig",
  death_col = "death",
  emig_col = "emig",
  year_beginning,
  final_year,
  combins = NULL
) {
  stopifnot(is.data.frame(data))
  if (missing(register_cols) || length(register_cols) == 0) {
    stop("register_cols must be provided.")
  }
  if (missing(covariate_cols)) {
    covariate_cols <- character(0)
  }

  all_cols <- c(id_col, year_col, register_cols, covariate_cols, reimmig_col, death_col, emig_col)
  missing_cols <- all_cols[!all_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (missing(year_beginning) || missing(final_year)) {
    stop("year_beginning and final_year must be provided.")
  }

  binary_cols <- c(register_cols, covariate_cols)
  grid_list <- rep(list(c(1, 0)), length(binary_cols))
  X <- expand.grid(grid_list)
  names(X) <- binary_cols
  X$row_id <- seq_len(nrow(X))

  num_combos <- nrow(X)

  X0 <- X
  X0[[reimmig_col]] <- 0
  X0$y <- X0$row_id

  X1 <- X
  X1[[reimmig_col]] <- 1
  X1$y <- X1$row_id + num_combos + 2

  X_all <- rbind(X0, X1)
  X_all$row_id <- NULL

  merge_cols <- c(register_cols, covariate_cols, reimmig_col)
  data <- merge(data, X_all, by = merge_cols, all.x = TRUE)

  if (any(is.na(data$y))) {
    stop("Some rows could not be matched to an observation combination. Check register/covariate coding.")
  }

  data$y[data[[death_col]] == 1] <- num_combos + 1
  data$y[data[[emig_col]] == 1] <- num_combos + 2

  data_wide <- reshape(
    data[, c(id_col, year_col, "y")],
    idvar = id_col,
    timevar = year_col,
    direction = "wide"
  )

  year_cols <- grep("^y\\.", names(data_wide), value = TRUE)
  years <- as.integer(sub("y\\.", "", year_cols))
  year_cols <- year_cols[order(years)]
  data_wide <- data_wide[, c(id_col, year_cols), drop = FALSE]

  y_matrix <- as.matrix(data_wide[, -1, drop = FALSE])

  if (!is.null(combins)) {
    if (!is.matrix(combins)) {
      stop("combins must be a matrix if provided.")
    }
    if (nrow(combins) != nrow(y_matrix)) {
      stop("combins must have the same number of rows as y_matrix.")
    }
    unobs_rows <- apply(X[, register_cols, drop = FALSE], 1, function(r) all(r == 0))
    cov_comb <- X[unobs_rows, covariate_cols, drop = FALSE]

    for (i in seq_len(nrow(y_matrix))) {
      first_year <- min(data[[year_col]][data[[id_col]] == data_wide[[id_col]][i]])
      for (t in (first_year - year_beginning + 1):ncol(y_matrix)) {
        if (is.na(y_matrix[i, t])) {
          cov <- cov_comb[1, , drop = FALSE]
          idx <- which(apply(cov_comb, 1, function(row) all(row == cov)))
          y_matrix[i, t] <- combins[i, t - 1] * 2^length(register_cols)
        }
      }
    }
  }

  list(
    y_matrix = y_matrix,
    X = X[, binary_cols, drop = FALSE],
    num_combos = num_combos
  )
}
