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
#'   \item{baseline_td_vars}{Character vector. Baseline time-dependent variables}
#'   \item{td_vars_by_time}{List. Time-dependent variables for each time point}
#' }
#'
#' @details
#' The function recognizes several naming patterns for time-dependent covariates:
#' \itemize{
#'   \item **known_with_underscore**: bmi_0, bp_0, bmi_1, bp_1, etc.
#'   \item **known_ending_zero**: bmi0, bp0, bmi1, bp1, etc.
#'   \item **generated_format**: waist0_1, waist0_2, waist1_1, waist1_2, etc.
#'   \item **generic_format**: td_covariate1_0, td_covariate2_0, etc.
#' }
#'
#' For exposure variables, it looks for the pattern: logM1_0, logM2_0, logM1_1, logM2_1, etc.
#'
#' @examples
#' \dontrun{
#' # Create test data in g-BKMR format
#' test_data <- data.frame(
#'   id = 1:100,
#'   sex = rbinom(100, 1, 0.5),
#'   bmi_0 = rnorm(100, 25, 3),
#'   bp_0 = rnorm(100, 120, 15),
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

  # Detect number of exposures per time point (p).
  p <- length(grep("^logM[0-9]+_0$", names(data)))

  # Current prepared data uses underscore-suffixed TD covariates:
  #   name_0, name_1, ..., name_(T-1)
  # Permit digits in the base name (for example L1_0 or bmi2_0),
  # but exclude mixture, outcome, id, sex, and generated baseline columns.
  baseline_td_vars <- grep("_0$", names(data), value = TRUE)
  baseline_td_vars <- baseline_td_vars[
    !grepl("^logM[0-9]+_0$", baseline_td_vars) &
      !baseline_td_vars %in% c("Y", "id", "sex") &
      !grepl("^baseline_[0-9]+$", baseline_td_vars)
  ]

  if (length(baseline_td_vars) == 0) {
    Ldim <- 0L
    td_covariate_names <- character(0)
    detected_pattern <- "none"
    td_vars_by_time <- list()
    warning("No baseline time-dependent covariates detected. Setting Ldim = 0.")
  } else {
    Ldim <- length(baseline_td_vars)
    td_covariate_names <- sub("_0$", "", baseline_td_vars)
    detected_pattern <- "underscore_suffix"
    td_vars_by_time <- list()

    if (T > 1) {
      for (t in seq_len(T - 1)) {
        td_vars_t <- paste0(td_covariate_names, "_", t)
        missing_t <- setdiff(td_vars_t, names(data))

        if (length(missing_t) > 0) {
          stop("Cannot detect time-dependent covariates for time ", t,
               ". Missing expected variable(s): ", paste(missing_t, collapse = ", "),
               ". Available variables: ", paste(names(data), collapse = ", "))
        }

        td_vars_by_time[[paste0("t", t)]] <- td_vars_t
      }
    }
  }

  result <- list(
    p = p,
    Ldim = Ldim,
    td_covariate_names = td_covariate_names,
    detected_pattern = detected_pattern,
    baseline_td_vars = baseline_td_vars,
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

  if (length(x$baseline_td_vars) > 0) {
    cat("Baseline time-dependent variables:", paste(x$baseline_td_vars, collapse = ", "), "\n")
  }

  if (length(x$td_vars_by_time) > 0) {
    cat("\nTime-dependent variables by time point:\n")
    for (i in seq_along(x$td_vars_by_time)) {
      time_label <- names(x$td_vars_by_time)[i]
      if (is.null(time_label) || is.na(time_label) || time_label == "") {
        time_label <- as.character(i)
      } else {
        time_label <- sub("^t", "", time_label)
      }
      cat("  Time", time_label, ":",
          paste(x$td_vars_by_time[[i]], collapse = ", "), "\n")
    }
  }

  invisible(x)
}
