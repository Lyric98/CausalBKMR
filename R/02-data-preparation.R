# Data preparation and formatting for g-BKMR

#' Prepare user matrices for g-BKMR analysis
#'
#' @description Converts user-provided matrices (Y, Z, X) into the wide-format
#' data structure required for g-BKMR analysis. Supports both continuous and
#' binary time-dependent covariates with enhanced input validation.
#'
#' @param Y Numeric vector. Outcome variable (length n).
#' @param Z Numeric matrix. Mixture exposure matrix (n x (Adim x T)).
#' @param X Numeric matrix. Covariate matrix (n x (Ldim x (T - 1) + baseline_covs)).
#' @param time_points Integer. Number of time points (T).
#' @param mixture_components Integer. Number of mixture components per time point (Adim).
#' @param td_covariates Integer. Number of time-dependent covariates per time point (Ldim).
#' @param baseline_covariates Integer. Number of baseline covariates (default: 1).
#' @param td_covariate_names Character vector. Names for time-dependent covariates (optional).
#' @param log_transform_mixtures Logical. Whether to log-transform mixture exposures (default: TRUE).
#' @param validate_input Logical. Whether to validate input dimensions (default: TRUE).
#'
#' @return A data frame in g-BKMR format with proper variable naming and metadata.
#'
#' @details
#' The function expects matrices organized as follows:
#' \itemize{
#'   \item Z matrix: Mixtures in chronological order (Mix1_T0, Mix2_T0, ..., MixAdim_T0, Mix1_T1, ...)
#'   \item X matrix: (TD_Cov1_T1, TD_Cov2_T1, ..., TD_CovLdim_T1, ..., TD_CovLdim_T(T-1), Baseline1, Baseline2, ...)
#' }
#'
#' The output data frame has the following structure:
#' \itemize{
#'   \item baseline_1, baseline_2, ...: Baseline covariates
#'   \item logM1_0, logM2_0, ...: Mixture exposures at time 0
#'   \item logM1_1, logM2_1, ...: Mixture exposures at time 1
#'   \item td_covariate1_1, td_covariate2_1, ...: Time-dependent covariates at time 1
#'   \item Y: Outcome variable
#'   \item id: Subject identifier
#' }
#'
#' @examples
#' \dontrun{
#' # Generate test data
#' n <- 200
#' Y <- rnorm(n)
#' Z <- matrix(rlnorm(n * 6), nrow = n, ncol = 6)  # 2 metals x 3 time points
#' X <- matrix(rnorm(n * 6), nrow = n, ncol = 6)   # 2 TD covs x 2 post-baseline times + 2 baseline
#'
#' # Prepare data
#' prepared_data <- prepare_gbkmr_data(
#'   Y = Y, Z = Z, X = X,
#'   time_points = 3,
#'   mixture_components = 2,
#'   td_covariates = 2,
#'   baseline_covariates = 2,
#'   td_covariate_names = c("bmi", "bp")
#' )
#'
#' # Check the structure
#' str(prepared_data)
#' head(prepared_data)
#' }
#'
#' @importFrom stats sd
#' @export
prepare_gbkmr_data <- function(
    Y,                     # Outcome vector (length n)
    Z,                     # Mixture exposure matrix (n x (Adim * T))
    X,                     # Covariate matrix (n x (Ldim * (T - 1) + baseline_covs))
    time_points,           # Number of time points (T)
    mixture_components,    # Number of mixture components per time point (Adim)
    td_covariates,         # Number of time-dependent covariates per time point (Ldim)
    baseline_covariates = 1,  # Number of baseline covariates
    td_covariate_names = NULL,  # User-provided names for time-dependent covariates
    log_transform_mixtures = TRUE,   # Whether to log-transform mixtures
    validate_input = TRUE            # Whether to validate input dimensions
) {

  # Get dimensions
  n <- length(Y)
  T <- time_points
  Adim <- mixture_components
  Ldim <- td_covariates
  n_baseline <- baseline_covariates

  # Convert X matrix to numeric if it contains non-numeric columns
  if (!is.numeric(X)) {
    X_numeric <- matrix(nrow = nrow(X), ncol = ncol(X))
    for (j in 1:ncol(X)) {
      col_data <- X[, j]
      col_numeric <- suppressWarnings(as.numeric(col_data))
      if (sum(is.na(col_numeric)) > sum(is.na(col_data))) {
        col_numeric <- as.numeric(factor(col_data))
      }
      X_numeric[, j] <- col_numeric
    }
    X <- X_numeric
  }

  if (!is.numeric(Z)) {
    Z <- apply(Z, 2, as.numeric)
  }

  # Input validation
  if (validate_input) {
    validate_user_matrices(Y, Z, X, T, Adim, Ldim, n_baseline)
  }

  # Determine naming strategy for time-dependent covariates
  use_user_names <- !is.null(td_covariate_names) && length(td_covariate_names) == Ldim

  if (use_user_names) {
    cat("Using user-provided time-dependent covariate names: ", paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    td_covariate_names <- if (Ldim > 0) paste0("td_covariate", seq_len(Ldim)) else character(0)
    cat("Using generic time-dependent covariate names: ", paste(td_covariate_names, collapse = ", "), "\n")
  }

  cat("Converting user matrices to g-BKMR format...\n")
  cat("Data dimensions: n =", n, ", T =", T, ", Adim =", Adim, ", Ldim =", Ldim, "\n")

  # Initialize the data frame
  df <- data.frame(id = 1:n)

  if (n_baseline < 1) {
    stop("At least one baseline covariate is required; no mock baseline covariate is created.")
  }

  baseline_start_col <- ncol(X) - n_baseline + 1
  baseline_data <- X[, baseline_start_col:ncol(X), drop = FALSE]
  for (i in seq_len(n_baseline)) {
    df[[paste0("baseline_", i)]] <- baseline_data[, i]
  }

  # Add mixture exposures in chronological order (logM1_0, logM2_0, ..., logMAdim_T-1)
  for (t in seq.int(0L, T - 1L)) {
    for (a in seq_len(Adim)) {
      col_idx <- t * Adim + a
      var_name <- paste0("logM", a, "_", t)

      mixture_data <- Z[, col_idx]

      # Log-transform if requested
      if (log_transform_mixtures) {
        if (any(mixture_data <= 0, na.rm = TRUE)) {
          positive_values <- mixture_data[is.finite(mixture_data) & mixture_data > 0]
          if (length(positive_values) == 0L) {
            stop("Cannot log-transform mixture column ", var_name,
                 ": it contains no positive finite values. ",
                 "Set log_transform_mixtures = FALSE or fix the exposure scale.")
          }
          min_positive <- min(positive_values)
          shift_value <- abs(min(mixture_data, na.rm = TRUE)) + min_positive * 0.001
          mixture_data <- log(mixture_data + shift_value)
          if (any(!is.finite(mixture_data))) {
            stop("Log transformation produced non-finite values for mixture column ",
                 var_name, ". Check for missing or infinite exposure values.")
          }
          warning("Some mixture values were <= 0. Added constant ", round(shift_value, 6), " before log transformation.")
        } else {
          mixture_data <- log(mixture_data)
          if (any(!is.finite(mixture_data))) {
            stop("Log transformation produced non-finite values for mixture column ",
                 var_name, ". Check for missing or infinite exposure values.")
          }
        }
      }

      df[[var_name]] <- mixture_data
    }
  }

  if (T > 1 && Ldim > 0) {
    for (t in seq_len(T - 1L)) {
      for (l in seq_len(Ldim)) {
        if (use_user_names) {
          var_name <- paste0(td_covariate_names[l], "_", t)
        } else {
          var_name <- paste0("td_covariate", l, "_", t)
        }

        col_idx <- (t - 1L) * Ldim + l

        if (col_idx <= (ncol(X) - n_baseline)) {
          df[[var_name]] <- X[, col_idx]
        } else {
          stop("Error in indexing time-dependent covariates. Check matrix dimensions.")
        }
      }
    }
  }

  # Add outcome
  df$Y <- Y

  ordered_cols <- paste0("baseline_", seq_len(n_baseline))

  # Add mixture exposures in chronological order
  for (t in seq.int(0L, T - 1L)) {
    for (a in seq_len(Adim)) {
      ordered_cols <- c(ordered_cols, paste0("logM", a, "_", t))
    }
  }

  # Add time-dependent covariates in chronological order
  if (T > 1 && Ldim > 0) {
    for (t in seq_len(T - 1L)) {
      if (use_user_names) {
        ordered_cols <- c(ordered_cols, paste0(td_covariate_names, "_", t))
      } else {
        ordered_cols <- c(ordered_cols, paste0("td_covariate", seq_len(Ldim), "_", t))
      }
    }
  }

  # Add outcome and id
  ordered_cols <- c(ordered_cols, "Y", "id")

  # Reorder the dataframe
  df <- df[, ordered_cols]

  # Add metadata
  attr(df, "data_info") <- list(
    n = n,
    time_points = T,
    mixture_components = Adim,
    td_covariates = Ldim,
    baseline_covariates = n_baseline,
    td_covariate_names = td_covariate_names,
    use_user_names = use_user_names,
    log_transformed = log_transform_mixtures,
    source = "user_matrices"
  )

  # Print summary
  cat("[OK] Data conversion successful!\n")
  cat("[OK] Generated", nrow(df), "x", ncol(df), "data frame\n")
  if (use_user_names) {
    cat("[OK] Using user covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  } else {
    cat("[OK] Using generic covariate names:", paste(td_covariate_names, collapse = ", "), "\n")
  }
  cat("[OK] Ready for g-BKMR analysis\n\n")

  return(df)
}
