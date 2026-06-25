test_that("validate_user_matrices accepts valid inputs", {
  set.seed(1)
  n <- 50
  Y <- rnorm(n)
  Z <- matrix(rnorm(n * 6), n, 6)
  X <- matrix(rnorm(n * 3), n, 3)

  expect_invisible(validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 1))
})

test_that("validate_user_matrices rejects non-numeric Y", {
  Y <- letters[1:10]
  Z <- matrix(1, 10, 6)
  X <- matrix(1, 10, 3)

  expect_error(validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 1),
               "numeric")
})

test_that("validate_user_matrices rejects wrong Z dimensions", {
  set.seed(1)
  n <- 30
  Y <- rnorm(n)
  Z <- matrix(rnorm(n * 4), n, 4)  # should be 6
  X <- matrix(rnorm(n * 3), n, 3)

  expect_error(validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 1),
               "dimension")
})

test_that("validate_user_matrices rejects wrong X dimensions", {
  set.seed(1)
  n <- 30
  Y <- rnorm(n)
  Z <- matrix(rnorm(n * 6), n, 6)
  X <- matrix(rnorm(n * 2), n, 2)  # should be 3

  expect_error(validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 1),
               "dimension")
})

test_that("validate_user_matrices warns on NA values", {
  set.seed(1)
  n <- 30
  Y <- rnorm(n); Y[1] <- NA
  Z <- matrix(rnorm(n * 6), n, 6)
  X <- matrix(rnorm(n * 3), n, 3)

  expect_warning(validate_user_matrices(Y, Z, X, T = 3, Adim = 2, Ldim = 1, n_baseline = 1),
                 "missing")
})

test_that("validate_user_matrices rejects zero baseline covariates", {
  n <- 10
  Y <- rnorm(n)
  Z <- matrix(rnorm(n * 4), n, 4)
  X <- matrix(rnorm(n), n, 1)

  expect_error(
    validate_user_matrices(Y, Z, X, T = 2, Adim = 2, Ldim = 1, n_baseline = 0),
    "baseline_covariates"
  )
})
