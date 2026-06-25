test_that("prepare_gbkmr_data returns correct dimensions", {
  set.seed(1)
  n <- 50; Tp <- 3; Adim <- 2; Ldim <- 1; nb <- 1
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * Adim * Tp)) + 0.01, n, Adim * Tp)
  X <- matrix(rnorm(n * (Ldim * (Tp - 1) + nb)), n, Ldim * (Tp - 1) + nb)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = Tp, mixture_components = Adim,
    td_covariates = Ldim, baseline_covariates = nb,
    td_covariate_names = "waist", log_transform_mixtures = TRUE)

  expect_equal(nrow(dat), n)
  # cols = baseline_1 + logM1_0..logM2_2 (6) + waist_1 + waist_2 + Y + id = 11
  expect_equal(ncol(dat), 11)
  expect_true("Y" %in% names(dat))
  expect_true("id" %in% names(dat))
  expect_true("baseline_1" %in% names(dat))
  expect_false("waist_0" %in% names(dat))
})

test_that("prepare_gbkmr_data uses generic names when none provided", {
  set.seed(1)
  n <- 30
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 4)) + 0.01, n, 4)
  X <- matrix(rnorm(n * 2), n, 2)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 2, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    log_transform_mixtures = TRUE)

  expect_true("td_covariate1_1" %in% names(dat))
  expect_false("td_covariate1_0" %in% names(dat))
})

test_that("prepare_gbkmr_data handles log_transform_mixtures = FALSE", {
  set.seed(1)
  n <- 30
  Y <- rnorm(n)
  Z <- matrix(rnorm(n * 4), n, 4)  # can have negatives
  X <- matrix(rnorm(n * 2), n, 2)

  dat <- prepare_gbkmr_data(Y, Z, X,
    time_points = 2, mixture_components = 2,
    td_covariates = 1, baseline_covariates = 1,
    log_transform_mixtures = FALSE)

  # raw values should be stored directly
  expect_equal(nrow(dat), n)
})

test_that("prepare_gbkmr_data rejects mixture columns with no positive values before log transform", {
  set.seed(1)
  n <- 10
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 4)) + 0.01, n, 4)
  Z[, 2] <- 0
  X <- matrix(rnorm(n * 2), n, 2)

  expect_error(
    prepare_gbkmr_data(Y, Z, X,
      time_points = 2, mixture_components = 2,
      td_covariates = 1, baseline_covariates = 1,
      log_transform_mixtures = TRUE),
    "contains no positive finite values"
  )
})

test_that("prepare_gbkmr_data rejects missing baseline covariates", {
  set.seed(1)
  n <- 20
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 4)) + 0.01, n, 4)
  X <- matrix(rnorm(n), n, 1)

  expect_error(
    prepare_gbkmr_data(Y, Z, X,
      time_points = 2, mixture_components = 2,
      td_covariates = 1, baseline_covariates = 0),
    "baseline_covariates"
  )
})
