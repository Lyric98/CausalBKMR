#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath("tools/generate-causal-plot-examples.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
local_library <- file.path(repo_root, ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

for (file in sort(list.files(file.path(repo_root, "R"),
                             pattern = "[.]R$", full.names = TRUE))) {
  source(file)
}

if (!requireNamespace("bkmr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Packages 'bkmr' and 'ggplot2' are required.")
}

set.seed(20260716)
n <- 180L
T_points <- 2L
p <- 3L

sex <- stats::rbinom(n, 1, 0.5)
age_z <- as.numeric(scale(stats::rnorm(n, 55, 9)))

correlation <- matrix(0.45, p, p)
diag(correlation) <- 1
log_a0 <- 1 + 0.55 * (matrix(stats::rnorm(n * p), n, p) %*%
                        chol(correlation))
x10 <- log_a0[, 1] - 1
x20 <- log_a0[, 2] - 1

waist_1 <- 0.60 * x10 - 0.35 * x20 + 0.45 * x10 * x20 +
  0.25 * sex + 0.15 * age_z + stats::rnorm(n, 0, 0.55)

visit_two_noise <- matrix(stats::rnorm(n * p), n, p) %*% chol(correlation)
log_a1 <- 1 + 0.50 * sweep(log_a0, 2, 1) +
  0.18 * waist_1 + 0.38 * visit_two_noise
x11 <- log_a1[, 1] - 1
x21 <- log_a1[, 2] - 1
x31 <- log_a1[, 3] - 1

Y <- 1.10 * x10 - 0.55 * x20 + 0.70 * x11 + 0.35 * x31^2 +
  1.15 * x10 * x20 + 0.55 * x11 * x21 + 0.40 * waist_1 +
  0.30 * sex + 0.20 * age_z + stats::rnorm(n, 0, 0.85)

Z <- cbind(exp(log_a0), exp(log_a1))
X <- cbind(waist_1, sex, age_z)
prepared <- prepare_gbkmr_data(
  Y = Y,
  Z = Z,
  X = X,
  time_points = T_points,
  mixture_components = p,
  td_covariates = 1,
  baseline_covariates = 2,
  td_covariate_names = "waist",
  log_transform_mixtures = TRUE
)

plot_sel <- seq(720, 1200, by = 16)
fit <- suppressWarnings(gbkmr_run(
  data = prepared,
  outcome = "Y",
  outcome_type = "continuous",
  time_points = T_points,
  currind = 2026,
  sel = plot_sel,
  K = 25,
  iter = 1200,
  n_knots = 25,
  engine = "bkmr",
  verbose = FALSE
))

plot_K <- 10L
overall <- gbkmr_causal_overall(
  fit,
  quantiles = seq(0.30, 0.90, by = 0.05),
  reference = 0.25,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

display_exposures <- c("logM1_0", "logM2_0", "logM1_1", "logM2_1")
univariate <- gbkmr_causal_univariate(
  fit,
  exposures = display_exposures,
  quantiles = seq(0.10, 0.90, by = 0.10),
  reference = 0.25,
  background = 0.50,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

pairs <- data.frame(
  focal = c("logM1_0", "logM1_1", "logM2_0"),
  conditional = c("logM2_0", "logM2_1", "logM3_0")
)
bivariate <- gbkmr_causal_bivariate(
  fit,
  exposures = unique(unlist(pairs)),
  pairs = pairs,
  quantiles = seq(0.10, 0.90, by = 0.10),
  conditional_quantiles = c(0.25, 0.50, 0.75),
  reference = 0.25,
  background = 0.50,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

iqr_effects <- gbkmr_causal_iqr(
  fit,
  contrast_quantiles = c(0.25, 0.75),
  background_quantiles = c(0.25, 0.50, 0.75),
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

interactions <- gbkmr_causal_interaction(
  fit,
  contrast_quantiles = c(0.25, 0.75),
  low_background = 0.25,
  high_background = 0.75,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

output_dir <- file.path(repo_root, "vignettes", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_example <- function(result, name, width, height) {
  ggplot2::ggsave(
    file.path(output_dir, paste0(name, ".png")),
    result$plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  paper_plot <- result$plot +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(plot.title = ggplot2::element_blank())
  ggplot2::ggsave(
    file.path(output_dir, paste0(name, "-paper.pdf")),
    paper_plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

save_example(overall, "causal-overall", 7.0, 4.4)
save_example(univariate, "causal-univariate", 8.0, 5.4)
save_example(bivariate, "causal-bivariate", 9.0, 4.8)
save_example(iqr_effects, "causal-iqr", 7.5, 5.3)
save_example(interactions, "causal-interaction", 7.5, 5.3)

message("Wrote causal plot examples to ", output_dir)
