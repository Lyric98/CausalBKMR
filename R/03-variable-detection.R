#' @file 03-variable-detection.R
#' @title Variable pattern detection for g-BKMR
#' @description Automatically detects variable structures in prepared g-BKMR data.
#' Identifies exposure variables, time-dependent covariates, and their naming patterns.

#' Detect variable patterns in g-BKMR data
#'
#' @description Automatically detects the structure of exposure and
#' time-dependent covariate variables in prepared g-BKMR data. This function
#' is essential for the automated analysis pipeline.
#'
#' @param data Data frame containing the variables to analyze. Should be in
#'   g-BKMR format (wide format with proper variable naming).
#' @param T Integer. Number of time points in the study.
#'
#' @return A list containing the detected variable structure:
#' \describe{
#'   \item{p}{Integer. Number of exposure variables per time point}
#'   \item{Ldim}{Integer. Number of time-dependent covariates per time point}
#'   \item{td_covariate_names}{Character vector. Names of time-dependent covariates}
#'   \item{detected_pattern}{Character. Pattern used for detection}
#'   \item{baseline_vars}{Character vector. Baseline covariates}
#'   \item{baseline_td_vars}{Character vector. Deprecated; always empty}
#'   \item{td_vars_by_time}{List. Time-dependent variables for each time point}
#' }
#'
#' @details
#' The function recognizes post-baseline time-dependent covariates named
#' \code{name_1}, \code{name_2}, ..., \code{name_(T-1)}.
#'
#' For exposure variables, it looks for the pattern: logM1_0, logM2_0, logM1_1, logM2_1, etc.
#'
#' @examples
#' \dontrun{
#' # Create test data in g-BKMR format
#' test_data <- data.frame(
#'   id = 1:100,
#'   baseline_1 = rbinom(100, 1, 0.5),
#'   logM1_0 = rnorm(100, 0, 1),
#'   logM2_0 = rnorm(100, 0, 1),
#'   logM1_1 = rnorm(100, 0, 1),
#'   logM2_1 = rnorm(100, 0, 1),
#'   bmi_1 = rnorm(100, 25, 3),
#'   bp_1 = rnorm(100, 120, 15),
#'   Y = rnorm(100, 0, 1)
#' )
#'
#' # Detect variable patterns
#' detection_result <- detect_variable_patterns(test_data, T = 2)
#'
#' # View results
#' print(detection_result)
#' }
#'
#' @export
detect_variable_patterns <- function(data, T) {

  p <- length(grep("^logM\\d+_0$", names(data)))
  if (p == 0L) stop("No exposure columns matching logM*_0 were detected.")

  for (t in seq.int(0L, T - 1L)) {
    expected_exposures <- paste0("logM", seq_len(p), "_", t)
    missing_exposures <- setdiff(expected_exposures, names(data))
    if (length(missing_exposures) > 0L) {
      stop("Missing exposure columns: ", paste(missing_exposures, collapse = ", "))
    }
  }

  baseline_vars <- grep("^baseline_[0-9]+$", names(data), value = TRUE)
  if (length(baseline_vars) == 0L) {
    stop("No baseline covariates detected. Expected columns named baseline_1, baseline_2, ...")
  }

  td_covariate_names <- character(0)
  detected_pattern <- "none"

  if (T > 1) {
    td_vars_t1 <- grep("_1$", names(data), value = TRUE)
    td_vars_t1 <- td_vars_t1[!grepl("^logM\\d+_1$", td_vars_t1)]
    td_vars_t1 <- td_vars_t1[!grepl("^baseline_[0-9]+$", td_vars_t1)]
    td_vars_t1 <- setdiff(td_vars_t1, c("Y", "id"))

    td_covariate_names <- sub("_1$", "", td_vars_t1)
    detected_pattern <- if (length(td_covariate_names) > 0L) "postbaseline_underscore" else "none"
  }

  Ldim <- length(td_covariate_names)
  td_vars_by_time <- list()

  if (Ldim > 0 && T > 1) {
    for (t in seq_len(T - 1L)) {
      td_vars_t <- paste0(td_covariate_names, "_", t)
      missing_td <- setdiff(td_vars_t, names(data))
      if (length(missing_td) > 0L) {
        stop("Cannot detect time-dependent covariates for time ", t,
             ". Missing: ", paste(missing_td, collapse = ", "),
             ". Available variables: ", paste(names(data), collapse = ", "))
      }

      td_vars_by_time[[paste0("t", t)]] <- td_vars_t
    }
  }

  result <- list(
    p = p,
    Ldim = Ldim,
    td_covariate_names = td_covariate_names,
    detected_pattern = detected_pattern,
    baseline_vars = baseline_vars,
    baseline_td_vars = character(0),
    td_vars_by_time = td_vars_by_time
  )

  # Add informative class for printing
  class(result) <- c("gbkmr_detection", "list")

  return(result)
}

#' Print method for gbkmr_detection objects
#'
#' @param x A gbkmr_detection object
#' @param ... Additional arguments (not used)
#'
#' @return Invisible x
#' @export
print.gbkmr_detection <- function(x, ...) {
  cat("g-BKMR Variable Detection Results\n")
  cat("================================\n")
  cat("Number of exposures per time point (p):", x$p, "\n")
  cat("Number of time-dependent covariates per time point (Ldim):", x$Ldim, "\n")

  if (length(x$td_covariate_names) > 0) {
    cat("Time-dependent covariate names:", paste(x$td_covariate_names, collapse = ", "), "\n")
  }

  cat("Detected naming pattern:", x$detected_pattern, "\n")

  if (length(x$baseline_vars) > 0) {
    cat("Baseline covariates:", paste(x$baseline_vars, collapse = ", "), "\n")
  }

  if (length(x$td_vars_by_time) > 0) {
    cat("\nTime-dependent variables by time point:\n")
    for (t in names(x$td_vars_by_time)) {
      cat("  Time", t, ":", paste(x$td_vars_by_time[[t]], collapse = ", "), "\n")
    }
  }

  invisible(x)
}
