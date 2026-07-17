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
n <- 500L
T_points <- 3L
p <- 5L

sex <- stats::rbinom(n, 1, 0.5)
age_z <- as.numeric(scale(stats::rnorm(n, 55, 9)))

correlation <- matrix(0.20, p, p)
diag(correlation) <- 1
log_a0 <- 1 + 0.50 * (matrix(stats::rnorm(n * p), n, p) %*%
                        chol(correlation))
centered_a0 <- sweep(log_a0, 2, 1)

waist_1 <- 0.55 * centered_a0[, 1] - 0.30 * centered_a0[, 2] +
  0.35 * centered_a0[, 1] * centered_a0[, 2] +
  0.25 * sex + 0.15 * age_z + stats::rnorm(n, 0, 0.55)

visit_two_noise <- matrix(stats::rnorm(n * p), n, p) %*% chol(correlation)
log_a1 <- 1 + 0.20 * centered_a0 + 0.16 * waist_1 +
  0.46 * visit_two_noise
centered_a1 <- sweep(log_a1, 2, 1)

waist_2 <- 0.55 * waist_1 + 0.45 * centered_a1[, 1] -
  0.30 * centered_a1[, 3] +
  0.35 * centered_a1[, 1] * centered_a1[, 2] +
  0.20 * sex + 0.10 * age_z + stats::rnorm(n, 0, 0.55)

visit_three_noise <- matrix(stats::rnorm(n * p), n, p) %*% chol(correlation)
log_a2 <- 1 + 0.20 * centered_a1 + 0.15 * waist_2 +
  0.46 * visit_three_noise
centered_a2 <- sweep(log_a2, 2, 1)

all_exposures <- cbind(centered_a0, centered_a1, centered_a2)
interaction_weights <- c(
  1, 1, 1, -1, -1,
  1, -1, 1, -1, 2,
  1, -1, -1, -1, 1
)
interaction_score <- drop(all_exposures %*% interaction_weights)

Y <- 20.00 * interaction_score +
  8.00 * (interaction_score^2 - mean(interaction_score^2)) +
  0.35 * waist_2 + 0.30 * sex + 0.20 * age_z +
  stats::rnorm(n, 0, 0.90)

Z <- cbind(exp(log_a0), exp(log_a1), exp(log_a2))
X <- cbind(waist_1, waist_2, sex, age_z)
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

plot_sel <- seq(1200, 2000, by = 40)
fit <- suppressWarnings(gbkmr_run(
  data = prepared,
  outcome = "Y",
  outcome_type = "continuous",
  time_points = T_points,
  currind = 2026,
  sel = plot_sel,
  K = 15,
  iter = 2000,
  n_knots = 40,
  engine = "bkmr",
  verbose = FALSE
))

plot_K <- 5L
overall <- gbkmr_causal_overall(
  fit,
  quantiles = seq(0.30, 0.90, by = 0.05),
  reference = 0.25,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

univariate <- gbkmr_causal_univariate(
  fit,
  quantiles = seq(0.10, 0.90, by = 0.10),
  reference = 0.25,
  background = 0.50,
  K = plot_K,
  sel = plot_sel,
  seed = 2026
)

bivariate <- lapply(seq_len(T_points) - 1L, function(time_index) {
  gbkmr_causal_bivariate(
    fit,
    time_points = time_index,
    layout = "matrix",
    quantiles = seq(0.10, 0.90, by = 0.10),
    conditional_quantiles = c(0.25, 0.50, 0.75),
    reference = 0.25,
    background = 0.50,
    K = plot_K,
    sel = plot_sel,
    seed = 2026
  )
})

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
save_example(univariate, "causal-univariate", 9.0, 8.0)
for (visit_index in seq_along(bivariate)) {
  save_example(
    bivariate[[visit_index]],
    paste0("causal-bivariate-visit", visit_index),
    10.0,
    8.2
  )
}
save_example(iqr_effects, "causal-iqr", 8.0, 8.0)
save_example(interactions, "causal-interaction", 8.0, 8.0)

message("Wrote causal plot examples to ", output_dir)
