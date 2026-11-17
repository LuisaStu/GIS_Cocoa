# Compare Ghana-wide major-road and urban-centre distances under four
# fixed-count nulls on the same national 1 km candidate grid.

source("R/01_setup.R")

candidate_path <- file.path(processed_dir, "ghana_1km_candidate_points.gpkg")
if (!file.exists(candidate_path)) {
  stop("Missing Ghana candidate grid; run R/12_ghana_wide_null_model.R first.")
}
candidate_points <- st_read(candidate_path, quiet = TRUE)
if (!isTRUE(st_crs(candidate_points) == st_crs(2136))) {
  stop("Ghana candidate grid must use EPSG:2136.")
}
required_fields <- c(
  "x", "y",
  "cocoa_positive", "environmental_weight", "environmental_weight_40km",
  "mean_temp_c", "log_precip", "elevation_m",
  "nearest_major_road_distance_m", "nearest_town_distance_m"
)
if (!all(required_fields %in% names(candidate_points))) {
  stop("Ghana candidate grid is missing predictors, scores, or distances.")
}

candidate_data <- st_drop_geometry(candidate_points)
candidate_n <- nrow(candidate_data)
observed_n <- sum(candidate_data$cocoa_positive)
if (observed_n == 0 || observed_n >= candidate_n) {
  stop("Ghana mapped-positive count is not valid for fixed-count sampling.")
}
if (any(!is.finite(candidate_data$nearest_major_road_distance_m)) ||
    any(!is.finite(candidate_data$nearest_town_distance_m))) {
  stop("Ghana road and town distances must be finite.")
}

fit_in_sample <- mgcv::bam(
  presence_response ~
    s(mean_temp_c, k = 5) +
    s(log_precip, k = 5) +
    s(elevation_m, k = 5),
  data = candidate_data,
  family = stats::binomial(),
  method = "fREML",
  select = TRUE,
  discrete = TRUE,
  nthreads = 2
)
in_sample_score <- as.numeric(stats::predict(
  fit_in_sample,
  newdata = candidate_data,
  type = "response"
))
if (any(!is.finite(in_sample_score))) {
  stop("The in-sample comparison model produced non-finite scores.")
}
in_sample_score <- pmax(in_sample_score, .Machine$double.eps)

auc_score <- function(observed, predicted) {
  positive_count <- as.double(sum(observed == 1L))
  negative_count <- as.double(sum(observed == 0L))
  if (positive_count == 0 || negative_count == 0) {
    return(NA_real_)
  }
  ranks <- rank(predicted, ties.method = "average")
  (sum(ranks[observed == 1L]) -
     positive_count * (positive_count + 1) / 2) /
    (positive_count * negative_count)
}

fold_diagnostics <- readr::read_csv(
  file.path(tables_dir, "ghana_environmental_cv_diagnostics.csv"),
  show_col_types = FALSE
)
crossfit_performance <- data.frame(
  block_size_m = c(20000, 40000),
  pooled_oof_auc = c(
    auc_score(candidate_data$cocoa_positive, candidate_data$environmental_score),
    auc_score(candidate_data$cocoa_positive, candidate_data$environmental_score_40km)
  ),
  pooled_oof_brier_score = c(
    mean((candidate_data$environmental_score - candidate_data$cocoa_positive)^2),
    mean((candidate_data$environmental_score_40km - candidate_data$cocoa_positive)^2)
  ),
  pooled_oof_log_loss = c(
    -mean(
      candidate_data$cocoa_positive *
        log(pmax(candidate_data$environmental_score, 1e-12)) +
        (1 - candidate_data$cocoa_positive) *
        log1p(-pmin(candidate_data$environmental_score, 1 - 1e-12))
    ),
    -mean(
      candidate_data$cocoa_positive *
        log(pmax(candidate_data$environmental_score_40km, 1e-12)) +
        (1 - candidate_data$cocoa_positive) *
        log1p(-pmin(candidate_data$environmental_score_40km, 1 - 1e-12))
    )
  ),
  mean_valid_fold_auc = vapply(c(20000, 40000), function(block_size) {
    mean(
      fold_diagnostics$heldout_auc[
        fold_diagnostics$block_size_m == block_size
      ],
      na.rm = TRUE
    )
  }, numeric(1)),
  folds_with_valid_auc = vapply(c(20000, 40000), function(block_size) {
    sum(
      is.finite(fold_diagnostics$heldout_auc) &
        fold_diagnostics$block_size_m == block_size
    )
  }, integer(1)),
  folds_without_positives = vapply(c(20000, 40000), function(block_size) {
    sum(fold_diagnostics$block_size_m == block_size &
          fold_diagnostics$heldout_positive_cells == 0L)
  }, integer(1)),
  stringsAsFactors = FALSE
)
write_csv(
  crossfit_performance,
  file.path(tables_dir, "ghana_environmental_cv_summary.csv")
)

weight_sets <- list(
  spatial_block_20km_crossfit = candidate_data$environmental_weight,
  spatial_block_40km_crossfit = candidate_data$environmental_weight_40km,
  in_sample_environmental = in_sample_score / sum(in_sample_score)
)
for (model_name in names(weight_sets)) {
  weights <- weight_sets[[model_name]]
  if (length(weights) != candidate_n ||
      any(!is.finite(weights)) ||
      any(weights <= 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("Invalid fixed-count sampling weights for ", model_name, ".")
  }
}

weighted_sample_without_replacement <- function(probabilities, size) {
  if (length(probabilities) != candidate_n ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0) ||
      abs(sum(probabilities) - 1) > 1e-8 ||
      size < 1 ||
      size >= length(probabilities)) {
    stop("Invalid input to weighted national sampling.")
  }
  uniforms <- pmax(stats::runif(candidate_n), .Machine$double.xmin)
  log_keys <- log(probabilities) - log(-log(uniforms))
  selected <- order(log_keys, decreasing = TRUE)[seq_len(size)]
  if (length(selected) != size || anyDuplicated(selected) > 0) {
    stop("Weighted national draw failed fixed-count/unique-cell validation.")
  }
  selected
}

observed_road_median <- median(
  candidate_data$nearest_major_road_distance_m[candidate_data$cocoa_positive]
)
observed_town_median <- median(
  candidate_data$nearest_town_distance_m[candidate_data$cocoa_positive]
)
n_simulations <- 999L
model_names <- c(names(weight_sets), "uniform_region")
primary_model <- "spatial_block_20km_crossfit"
# Each eligible 1 km cell contributes 1 km² of grid support to its quadrat.
quadrat_width_m <- 25000
grid_to_m <- 0.304799710181509
quadrat_x <- floor(candidate_data$x * grid_to_m / quadrat_width_m)
quadrat_y <- floor(candidate_data$y * grid_to_m / quadrat_width_m)
quadrat_key <- paste(quadrat_x, quadrat_y)
representatives <- which(!duplicated(quadrat_key))
quadrat_id <- match(quadrat_key, quadrat_key[representatives])
quadrat_n <- length(representatives)
quadrat_area <- tabulate(quadrat_id, nbins = quadrat_n)
quadrat_observed <- tabulate(
  quadrat_id[candidate_data$cocoa_positive], nbins = quadrat_n
)
quadrat_draws <- matrix(0L, nrow = n_simulations, ncol = quadrat_n)
distance_breaks <- c(-Inf, 5000, 10000, 20000, 50000, Inf)
distance_labels <- c("0-5 km", "5-10 km", "10-20 km", "20-50 km", "50+ km")
profile_data <- list(
  nearest_major_road = candidate_data$nearest_major_road_distance_m,
  nearest_town_city = candidate_data$nearest_town_distance_m
)
profile_bins <- lapply(profile_data, function(distance) {
  cut(
    distance,
    breaks = distance_breaks,
    labels = distance_labels,
    right = TRUE,
    include.lowest = TRUE
  )
})
profile_candidate_counts <- lapply(profile_bins, function(bin) {
  tabulate(as.integer(bin), nbins = length(distance_labels))
})
profile_observed_counts <- lapply(profile_bins, function(bin) {
  tabulate(
    as.integer(bin[candidate_data$cocoa_positive]),
    nbins = length(distance_labels)
  )
})
profile_draw_counts <- array(
  NA_integer_,
  dim = c(
    n_simulations,
    length(model_names),
    length(distance_labels),
    length(profile_data)
  ),
  dimnames = list(
    NULL,
    model_names,
    distance_labels,
    names(profile_data)
  )
)
road_draws <- matrix(
  NA_real_,
  nrow = n_simulations,
  ncol = length(model_names),
  dimnames = list(NULL, model_names)
)
town_draws <- road_draws
inclusion_counts <- lapply(model_names, function(model) integer(candidate_n))
names(inclusion_counts) <- model_names

set.seed(20261118)
for (simulation in seq_len(n_simulations)) {
  selected_by_model <- vector("list", length(model_names))
  names(selected_by_model) <- model_names
  for (model_name in names(weight_sets)) {
    selected_by_model[[model_name]] <- weighted_sample_without_replacement(
      weight_sets[[model_name]],
      observed_n
    )
  }
  selected_by_model$uniform_region <- sample.int(
    candidate_n,
    size = observed_n,
    replace = FALSE
  )
  quadrat_draws[simulation, ] <- tabulate(
    quadrat_id[selected_by_model[[primary_model]]], nbins = quadrat_n
  )

  for (model_name in model_names) {
    selected <- selected_by_model[[model_name]]
    road_draws[simulation, model_name] <- median(
      candidate_data$nearest_major_road_distance_m[selected]
    )
    town_draws[simulation, model_name] <- median(
      candidate_data$nearest_town_distance_m[selected]
    )
    for (statistic_name in names(profile_data)) {
      profile_draw_counts[
        simulation,
        model_name,
        ,
        statistic_name
      ] <- tabulate(
        as.integer(profile_bins[[statistic_name]][selected]),
        nbins = length(distance_labels)
      )
    }
    inclusion_counts[[model_name]][selected] <-
      inclusion_counts[[model_name]][selected] + 1L
  }
}
stopifnot(all(rowSums(quadrat_draws) == observed_n))
quadrat_quantiles <- apply(
  sweep(quadrat_draws, 2, quadrat_area, "/"),
  2, stats::quantile, probs = c(0.025, 0.5, 0.975)
)
write_csv(
  data.frame(
    x_m = (quadrat_x[representatives] + 0.5) * quadrat_width_m,
    y_m = (quadrat_y[representatives] + 0.5) * quadrat_width_m,
    eligible_grid_area_km2 = quadrat_area,
    observed_intensity = quadrat_observed / quadrat_area,
    simulated_lower = quadrat_quantiles[1, ],
    simulated_median = quadrat_quantiles[2, ],
    simulated_upper = quadrat_quantiles[3, ]
  ),
  file.path(tables_dir, "ghana_quadrat_intensity.csv")
)
if (any(!is.finite(road_draws)) ||
    any(!is.finite(town_draws)) ||
    any(vapply(inclusion_counts, sum, numeric(1)) !=
          observed_n * n_simulations)) {
  stop("National simulations failed finite-distance or fixed-count validation.")
}

summarize_statistic <- function(model_name, statistic_name, observed, draws) {
  quantiles <- unname(stats::quantile(draws, c(0.025, 0.5, 0.975)))
  data.frame(
    null_model = model_name,
    statistic = statistic_name,
    observed_median_m = observed,
    simulated_median_2_5pct = quantiles[[1]],
    simulated_median_50pct = quantiles[[2]],
    simulated_median_97_5pct = quantiles[[3]],
    observed_relative_to_null = observed / quantiles[[2]],
    lower_tail_monte_carlo_p = (1 + sum(draws <= observed)) / (n_simulations + 1),
    simulation_count = n_simulations,
    fixed_observation_count = observed_n,
    candidate_cell_count = candidate_n,
    stringsAsFactors = FALSE
  )
}

summary_rows <- lapply(model_names, function(model_name) {
  rbind(
    summarize_statistic(
      model_name,
      "nearest_major_road",
      observed_road_median,
      road_draws[, model_name]
    ),
    summarize_statistic(
      model_name,
      "nearest_town_city",
      observed_town_median,
      town_draws[, model_name]
    )
  )
})
simulation_summary <- do.call(rbind, summary_rows)
write_csv(
  simulation_summary,
  file.path(tables_dir, "ghana_simulation_summary.csv")
)

simulation_draws <- do.call(rbind, lapply(model_names, function(model_name) {
  data.frame(
    simulation_id = seq_len(n_simulations),
    null_model = model_name,
    median_major_road_distance_m = road_draws[, model_name],
    median_town_city_distance_m = town_draws[, model_name],
    stringsAsFactors = FALSE
  )
}))
write_csv(
  simulation_draws,
  file.path(tables_dir, "ghana_simulation_draws.csv")
)

profile_draws <- do.call(rbind, lapply(names(profile_data), function(statistic_name) {
  do.call(rbind, lapply(model_names, function(model_name) {
    simulated_counts <- profile_draw_counts[
      ,
      model_name,
      ,
      statistic_name
    ]
    data.frame(
      null_model = model_name,
      statistic = statistic_name,
      simulation_id = rep(seq_len(n_simulations), each = length(distance_labels)),
      distance_bin = rep(distance_labels, times = n_simulations),
      candidate_cell_count = rep(
        profile_candidate_counts[[statistic_name]],
        times = n_simulations
      ),
      observed_positive_cells = rep(
        profile_observed_counts[[statistic_name]],
        times = n_simulations
      ),
      observed_intensity = rep(
        profile_observed_counts[[statistic_name]] /
          profile_candidate_counts[[statistic_name]],
        times = n_simulations
      ),
      simulated_positive_cells = as.vector(t(simulated_counts)),
      simulated_intensity = as.vector(t(simulated_counts)) /
        rep(profile_candidate_counts[[statistic_name]], times = n_simulations),
      stringsAsFactors = FALSE
    )
  }))
}))
if (nrow(profile_draws) !=
    n_simulations * length(model_names) *
      length(distance_labels) * length(profile_data) ||
    any(!is.finite(profile_draws$simulated_intensity))) {
  stop("Distance-binned intensity profile validation failed.")
}
write_csv(
  profile_draws,
  file.path(tables_dir, "ghana_distance_binned_intensity_draws.csv")
)

candidate_data$observed_presence <- as.integer(candidate_data$cocoa_positive)
for (model_name in model_names) {
  candidate_data[[paste0(model_name, "_simulated_inclusion")]] <-
    inclusion_counts[[model_name]] / n_simulations
}
intensity_columns <- c(
  "cell", "x", "y", "cocoa_fraction", "observed_presence",
  paste0(model_names, "_simulated_inclusion")
)
intensity_data <- candidate_data[intensity_columns]
intensity_data$x_m <- intensity_data$x * 0.304799710181509
intensity_data$y_m <- intensity_data$y * 0.304799710181509
write_csv(
  intensity_data,
  file.path(tables_dir, "ghana_intensity_comparison.csv")
)

message("Ghana-wide 999-draw simulations complete.")
print(simulation_summary)
