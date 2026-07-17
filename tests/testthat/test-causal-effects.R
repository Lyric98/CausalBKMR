test_that("causal effect helpers preserve exposure-specific quantiles", {
  exposure_data <- cbind(
    logM1_0 = 1:8,
    logM2_0 = seq(10, 80, by = 10)
  )
  raw <- list(
    meta = list(exposure_names = colnames(exposure_data)),
    gcomp_state = list(exposure_data = exposure_data)
  )

  observed <- causalBKMR:::.gbkmr_exposure_quantiles(raw, c(0.25, 0.75))
  expect_equal(unname(observed[, "logM1_0"]),
               as.numeric(quantile(exposure_data[, "logM1_0"], c(0.25, 0.75))))
  expect_equal(unname(observed[, "logM2_0"]),
               as.numeric(quantile(exposure_data[, "logM2_0"], c(0.25, 0.75))))
})

test_that("posterior summaries use draw-wise credible intervals", {
  draws <- c(-2, -1, 0, 1, 4)
  observed <- causalBKMR:::.gbkmr_draw_interval(draws, 0.80)

  expect_equal(unname(observed["estimate"]), mean(draws))
  expect_equal(unname(observed["lower"]), unname(quantile(draws, 0.10)))
  expect_equal(unname(observed["upper"]), unname(quantile(draws, 0.90)))
})

test_that("bivariate pair helper supports directed matrix layout", {
  selected <- c("M1", "M2", "M3")
  directed <- causalBKMR:::.gbkmr_pairs(selected, NULL, ordered = TRUE)
  compact <- causalBKMR:::.gbkmr_pairs(selected, NULL, ordered = FALSE)

  expect_equal(nrow(directed), 6L)
  expect_false(any(directed$focal == directed$conditional))
  expect_true(all(c("M1 | M2", "M2 | M1") %in%
                    paste(directed$focal, directed$conditional, sep = " | ")))
  expect_equal(nrow(compact), 3L)
})

test_that("causal plotting commands use full longitudinal interventions", {
  skip_if_not_installed("bkmr")
  skip_if_not_installed("ggplot2")
  set.seed(91)
  n <- 40
  Y <- rnorm(n)
  Z <- matrix(abs(rnorm(n * 4)) + 0.01, n, 4)
  X <- cbind(rnorm(n), rbinom(n, 1, 0.5))

  dat <- prepare_gbkmr_data(
    Y, Z, X,
    time_points = 2,
    mixture_components = 2,
    td_covariates = 1,
    baseline_covariates = 1,
    td_covariate_names = "waist",
    log_transform_mixtures = TRUE
  )
  fit <- suppressWarnings(gbkmr_run(
    data = dat,
    time_points = 2,
    iter = 60,
    sel = c(40, 50, 60),
    n = 30,
    K = 1,
    n_knots = 8,
    engine = "bkmr",
    verbose = FALSE
  ))
  exposure_names <- fit$raw_results$meta$exposure_names
  focal <- exposure_names[1]
  conditional <- exposure_names[2]
  exposure_data <- fit$raw_results$gcomp_state$exposure_data
  q_value <- function(name, probability) {
    unname(quantile(exposure_data[, name], probability))
  }

  overall <- gbkmr_causal_overall(
    fit, quantiles = c(0.25, 0.75), K = 1, seed = 12
  )
  univariate <- gbkmr_causal_univariate(
    fit, exposures = focal, quantiles = c(0.25, 0.75), K = 1, seed = 12
  )
  bivariate <- gbkmr_causal_bivariate(
    fit,
    exposures = c(focal, conditional),
    layout = "matrix",
    quantiles = c(0.25, 0.75),
    K = 1,
    seed = 12
  )
  iqr <- gbkmr_causal_iqr(
    fit, exposures = focal, K = 1, seed = 12
  )
  interaction <- gbkmr_causal_interaction(
    fit, exposures = focal, K = 1, seed = 12
  )

  results <- list(overall, univariate, bivariate, iqr, interaction)
  expect_true(all(vapply(results, inherits, logical(1), "gbkmr_causal_effects")))
  expect_true(all(vapply(results, function(x) inherits(x$plot, "ggplot"), logical(1))))
  expect_true(all(vapply(results, function(x) nrow(x$draws) > 0, logical(1))))

  expect_equal(overall$summary$estimate[overall$summary$quantile == 0.25], 0)
  nonfocal <- setdiff(exposure_names, focal)
  expect_true(all(vapply(nonfocal, function(name) {
    all(univariate$interventions[[name]] == q_value(name, 0.50))
  }, logical(1))))

  remaining <- setdiff(exposure_names, c(focal, conditional))
  expect_true(all(vapply(remaining, function(name) {
    all(bivariate$interventions[[name]] == q_value(name, 0.50))
  }, logical(1))))
  expect_equal(sort(unique(bivariate$interventions$conditional_quantile)),
               c(0.25, 0.50, 0.75))
  expect_equal(nrow(unique(bivariate$summary[c("focal", "conditional")])), 2L)
  expect_s3_class(bivariate$plot$facet, "FacetGrid")
  expect_identical(bivariate$settings$layout, "matrix")

  iqr_low_background <- iqr$interventions$background_quantile == 0.25
  expect_true(all(vapply(nonfocal, function(name) {
    all(iqr$interventions[[name]][iqr_low_background] == q_value(name, 0.25))
  }, logical(1))))

  low <- iqr$draws[iqr$draws$background_quantile == 0.25, ]
  high <- iqr$draws[iqr$draws$background_quantile == 0.75, ]
  low <- low[order(low$.draw), ]
  high <- high[order(high$.draw), ]
  interaction_draws <- interaction$draws[order(interaction$draws$.draw), ]
  expect_equal(interaction_draws$interaction_contrast, high$effect - low$effect)
  expect_equal(
    interaction$summary$lower,
    unname(quantile(interaction_draws$interaction_contrast, 0.025))
  )
  expect_equal(
    interaction$summary$upper,
    unname(quantile(interaction_draws$interaction_contrast, 0.975))
  )
})
