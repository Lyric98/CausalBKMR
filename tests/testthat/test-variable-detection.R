test_that("detect_variable_patterns finds user-named covariates", {
  set.seed(1)
  n <- 50
  dat <- data.frame(
    id = 1:n, baseline_1 = rbinom(n, 1, 0.5),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n), waist_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 2)
  expect_equal(det$Ldim, 1)
  expect_equal(det$td_covariate_names, "waist")
  expect_equal(det$baseline_vars, "baseline_1")
  expect_s3_class(det, "gbkmr_detection")
})

test_that("detect_variable_patterns finds generic covariates", {
  set.seed(1)
  n <- 30
  dat <- data.frame(
    id = 1:n, baseline_1 = rbinom(n, 1, 0.5),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n),
    td_covariate1_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 2)
  expect_equal(det$Ldim, 1)
  expect_equal(det$detected_pattern, "postbaseline_underscore")
})

test_that("detect_variable_patterns counts exposures correctly", {
  set.seed(1)
  n <- 20
  dat <- data.frame(
    id = 1:n, baseline_1 = rbinom(n, 1, 0.5),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n), logM3_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n), logM3_1 = rnorm(n),
    bmi_1 = rnorm(n), bp_1 = rnorm(n),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 3)
  expect_equal(det$Ldim, 2)
})

test_that("detect_variable_patterns handles digits in TD covariate names", {
  set.seed(1)
  n <- 20
  dat <- data.frame(
    id = 1:n, baseline_1 = rbinom(n, 1, 0.5),
    logM1_0 = rnorm(n), logM2_0 = rnorm(n),
    logM1_1 = rnorm(n), logM2_1 = rnorm(n), L1_1 = rbinom(n, 1, 0.5),
    Y = rnorm(n)
  )

  det <- detect_variable_patterns(dat, T = 2)
  expect_equal(det$p, 2)
  expect_equal(det$td_covariate_names, "L1")
  expect_equal(det$td_vars_by_time$t1, "L1_1")
})

test_that("detect_variable_patterns rejects missing baseline covariates", {
  n <- 20
  dat <- data.frame(
    id = 1:n,
    logM1_0 = rnorm(n), logM1_1 = rnorm(n),
    waist_1 = rnorm(n), Y = rnorm(n)
  )

  expect_error(detect_variable_patterns(dat, T = 2), "baseline")
})
