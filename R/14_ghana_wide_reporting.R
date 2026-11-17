# Validate national simulation outputs and build Ghana-wide result figures.

source("R/01_setup.R")

candidate_path <- file.path(processed_dir, "ghana_1km_candidate_points.gpkg")
if (!file.exists(candidate_path)) {
  stop("Missing Ghana candidate grid; run R/12_ghana_wide_null_model.R first.")
}
candidate_points <- st_read(candidate_path, quiet = TRUE)
if (!isTRUE(st_crs(candidate_points) == st_crs(2136))) {
  stop("Ghana candidate grid must use EPSG:2136.")
}
candidate_data <- st_drop_geometry(candidate_points)
required_candidate_fields <- c(
  "cocoa_positive", "environmental_score", "environmental_score_40km"
)
if (!all(required_candidate_fields %in% names(candidate_data))) {
  stop("Ghana candidate grid is missing mapped outcomes or cross-fitted scores.")
}

candidate_n <- nrow(candidate_data)
observed_n <- sum(candidate_data$cocoa_positive)
if (candidate_n != 240294L || observed_n != 57131L) {
  stop(
    "Unexpected Ghana candidate/positive counts: ",
    candidate_n, " candidates and ", observed_n, " positives."
  )
}

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

fold_diagnostics <- read_csv(
  file.path(tables_dir, "ghana_environmental_cv_diagnostics.csv"),
  show_col_types = FALSE
)
performance_rows <- lapply(c(20000, 40000), function(block_size) {
  score <- if (block_size == 20000) {
    candidate_data$environmental_score
  } else {
    candidate_data$environmental_score_40km
  }
  observed <- as.integer(candidate_data$cocoa_positive)
  selected_folds <- fold_diagnostics[
    fold_diagnostics$block_size_m == block_size,
    ,
    drop = FALSE
  ]
  if (nrow(selected_folds) != 5L) {
    stop("Expected five spatial-fold diagnostics at ", block_size, " m.")
  }
  clipped <- pmin(pmax(score, 1e-12), 1 - 1e-12)
  data.frame(
    block_size_m = block_size,
    pooled_oof_auc = auc_score(observed, score),
    pooled_oof_brier_score = mean((score - observed)^2),
    pooled_oof_log_loss = -mean(
      observed * log(clipped) + (1 - observed) * log1p(-clipped)
    ),
    mean_valid_fold_auc = mean(selected_folds$heldout_auc, na.rm = TRUE),
    valid_auc_folds = sum(is.finite(selected_folds$heldout_auc)),
    folds_without_positives = sum(selected_folds$heldout_positive_cells == 0L),
    fold_auc_range = paste(
      round(range(selected_folds$heldout_auc, na.rm = TRUE), 3),
      collapse = " to "
    ),
    stringsAsFactors = FALSE
  )
})
crossfit_performance <- do.call(rbind, performance_rows)
write_csv(
  crossfit_performance,
  file.path(tables_dir, "ghana_environmental_cv_summary.csv")
)

admin <- st_read(
  file.path(raw_dir, "ghana.gpkg"),
  layer = "gis_osm_adminareas_a_free",
  quiet = TRUE
)
ghana <- admin[
  admin$fclass == "national" & admin$name == "Ghana",
  ,
  drop = FALSE
]
if (nrow(ghana) != 1L) {
  stop("Expected one Ghana national boundary for the coverage audit.")
}
ghana_union <- st_union(st_transform(st_make_valid(ghana), 2136))

osm_path <- file.path(raw_dir, "ghana.gpkg")
roads <- st_read(
  osm_path,
  query = paste(
    "SELECT * FROM gis_osm_roads_free",
    "WHERE fclass IN ('motorway', 'trunk', 'primary', 'secondary')"
  ),
  quiet = TRUE
)
roads <- st_collection_extract(
  st_intersection(st_transform(roads, 2136), ghana_union),
  "LINESTRING", warn = FALSE
)
places <- st_read(osm_path, layer = "gis_osm_places_free", quiet = TRUE)
places <- places[places$fclass %in% c("city", "town", "national_capital"), ]
places <- st_filter(st_transform(places, 2136), ghana_union, .predicate = st_within)
cocoa_cells <- candidate_data[candidate_data$cocoa_positive, c("x", "y")]
# Simplify lines for display at national scale. The distance analysis uses full geometry.
display_boundary <- st_simplify(ghana_union, dTolerance = units::set_units(150, "m"))
display_roads <- st_simplify(roads, dTolerance = units::set_units(100, "m"))

save_opening_map <- function() {
  grDevices::png(file.path(figures_dir, "ghana_cocoa_access_map.png"),
                 width = 6.5, height = 7.5, units = "in", res = 250,
                 type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(2.7, 0.3, 0.3, 0.3), xpd = NA)
  plot(display_boundary, col = "#F5F3EC", border = "#555555", lwd = 0.6)
  graphics::points(cocoa_cells$x, cocoa_cells$y, col = "#408653", pch = 15, cex = 0.16)
  plot(st_geometry(display_roads), add = TRUE, col = "#A16D4D", lwd = 0.55)
  plot(st_geometry(places), add = TRUE, col = "#252249", pch = 16, cex = 0.35)
  bounds <- graphics::par("usr")
  graphics::legend(
    mean(bounds[1:2]), bounds[3], xjust = 0.5, yjust = 1,
    legend = c("Mapped cocoa", "Major road", "Urban centre"),
    col = c("#408653", "#A16D4D", "#252249"), pch = c(15, NA, 16),
    lty = c(NA, 1, NA), horiz = TRUE, bty = "n"
  )
}
save_opening_map()

cocoa_raster <- terra::rast(file.path(raw_dir, "cocoa_map_th_065.tif"))
cocoa_extent <- as.vector(terra::ext(cocoa_raster))
cocoa_bbox <- st_as_sfc(st_bbox(
  c(
    xmin = cocoa_extent[["xmin"]],
    ymin = cocoa_extent[["ymin"]],
    xmax = cocoa_extent[["xmax"]],
    ymax = cocoa_extent[["ymax"]]
  ),
  crs = st_crs(4326)
))
coverage_geometry <- st_transform(cocoa_bbox, 2136)
boundary_area_km2 <- as.numeric(units::set_units(st_area(ghana_union), "km^2"))
coverage_area_km2 <- as.numeric(units::set_units(
  st_area(st_intersection(ghana_union, coverage_geometry)),
  "km^2"
))
outside_area_km2 <- as.numeric(units::set_units(
  st_area(st_difference(ghana_union, coverage_geometry)),
  "km^2"
))
coverage_summary <- data.frame(
  boundary_area_km2 = boundary_area_km2,
  raster_bbox_intersection_km2 = coverage_area_km2,
  outside_raster_bbox_km2 = outside_area_km2,
  raster_bbox_intersection_pct = 100 * coverage_area_km2 / boundary_area_km2,
  outside_raster_bbox_pct = 100 * outside_area_km2 / boundary_area_km2,
  cocoa_raster_west_lon = cocoa_extent[["xmin"]],
  cocoa_raster_east_lon = cocoa_extent[["xmax"]],
  boundary_source = "Geofabrik Ghana OSM national administrative polygon",
  stringsAsFactors = FALSE
)

quadrats <- read_csv(
  file.path(tables_dir, "ghana_quadrat_intensity.csv"), show_col_types = FALSE
)
quadrat_surfaces <- rbind(
  transform(quadrats, panel = "Observed intensity", intensity = observed_intensity),
  transform(quadrats, panel = "Simulated median", intensity = simulated_median)
)
quadrat_surfaces$panel <- factor(
  quadrat_surfaces$panel, levels = c("Observed intensity", "Simulated median")
)
quadrats$departure <- factor(
  case_when(
    quadrats$observed_intensity > quadrats$simulated_upper ~ "Above interval",
    quadrats$observed_intensity < quadrats$simulated_lower ~ "Below interval",
    TRUE ~ "Within interval"
  ),
  levels = c("Below interval", "Within interval", "Above interval")
)
ring_coordinates <- st_coordinates(st_cast(display_boundary, "MULTIPOLYGON"))
boundary_m <- data.frame(
  x_m = ring_coordinates[, "X"] * 0.304799710181509,
  y_m = ring_coordinates[, "Y"] * 0.304799710181509,
  ring = interaction(ring_coordinates[, "L1"], ring_coordinates[, "L2"])
)
quadrat_theme <- theme_void() + theme(
  plot.title = element_text(face = "bold"), legend.position = "bottom"
)
quadrat_map <- ggplot(quadrat_surfaces, aes(x_m, y_m, fill = intensity)) +
  geom_tile(width = 25000, height = 25000) +
  geom_path(data = boundary_m, aes(x_m, y_m, group = ring),
            inherit.aes = FALSE, linewidth = 0.25) +
  facet_wrap(~panel) +
  coord_equal() +
  scale_fill_viridis_c(limits = c(0, 1), name = "Positive cells per eligible km²") +
  quadrat_theme
ggsave(file.path(figures_dir, "ghana_quadrat_intensity.png"), quadrat_map,
       width = 8, height = 5, dpi = 250, bg = "white")
quadrat_departure <- ggplot(quadrats, aes(x_m, y_m, fill = departure)) +
  geom_tile(width = 25000, height = 25000) +
  geom_path(data = boundary_m, aes(x_m, y_m, group = ring),
            inherit.aes = FALSE, linewidth = 0.25) +
  coord_equal() +
  scale_fill_manual(values = c("#2166AC", "#EDEDED", "#B2182B"), drop = FALSE) +
  labs(fill = NULL) +
  quadrat_theme
ggsave(file.path(figures_dir, "ghana_quadrat_intervals.png"), quadrat_departure,
       width = 6, height = 5, dpi = 250, bg = "white")
write_csv(
  coverage_summary,
  file.path(tables_dir, "ghana_source_coverage_summary.csv")
)

simulation_summary <- read_csv(
  file.path(tables_dir, "ghana_simulation_summary.csv"),
  show_col_types = FALSE
)
simulation_draws <- read_csv(
  file.path(tables_dir, "ghana_simulation_draws.csv"),
  show_col_types = FALSE
)
profile_draws <- read_csv(
  file.path(tables_dir, "ghana_distance_binned_intensity_draws.csv"),
  show_col_types = FALSE
)
intensity <- read_csv(
  file.path(tables_dir, "ghana_intensity_comparison.csv"),
  show_col_types = FALSE
)
model_names <- c(
  "spatial_block_20km_crossfit",
  "spatial_block_40km_crossfit",
  "in_sample_environmental",
  "uniform_region"
)
if (nrow(simulation_summary) != 8L ||
    nrow(simulation_draws) != 999L * length(model_names) ||
    nrow(profile_draws) != 999L * length(model_names) * 5L * 2L ||
    nrow(intensity) != candidate_n ||
    any(simulation_summary$simulation_count != 999L) ||
    any(simulation_summary$candidate_cell_count != candidate_n) ||
    any(simulation_summary$fixed_observation_count != observed_n) ||
    !all(model_names %in% unique(simulation_draws$null_model))) {
  stop("National simulation outputs have unexpected dimensions or models.")
}
profile_keys <- c("null_model", "statistic", "simulation_id", "distance_bin")
if (anyDuplicated(profile_draws[profile_keys]) > 0) {
  stop("Distance-binned intensity output contains duplicate draw/bin rows.")
}
profile_totals <- profile_draws |>
  dplyr::group_by(null_model, statistic, simulation_id) |>
  dplyr::summarise(
    simulated_positive_cells = sum(simulated_positive_cells),
    observed_positive_cells = sum(observed_positive_cells),
    candidate_cell_count = sum(candidate_cell_count),
    .groups = "drop"
  )
if (nrow(profile_totals) != 999L * length(model_names) * 2L ||
    any(profile_totals$simulated_positive_cells != observed_n) ||
    any(profile_totals$observed_positive_cells != observed_n) ||
    any(profile_totals$candidate_cell_count != candidate_n)) {
  stop("Distance-binned intensity draws failed fixed-count or bin-total validation.")
}

qa_rows <- lapply(model_names, function(model_name) {
  model_draws <- simulation_draws[
    simulation_draws$null_model == model_name,
    ,
    drop = FALSE
  ]
  inclusion_name <- paste0(model_name, "_simulated_inclusion")
  if (!inclusion_name %in% names(intensity)) {
    stop("Missing simulated inclusion surface for ", model_name, ".")
  }
  inclusion_sum <- sum(intensity[[inclusion_name]])
  finite_road <- all(is.finite(model_draws$median_major_road_distance_m))
  finite_town <- all(is.finite(model_draws$median_town_city_distance_m))
  unique_draws <- length(unique(model_draws$simulation_id))
  if (nrow(model_draws) != 999L ||
      unique_draws != 999L ||
      !finite_road ||
      !finite_town ||
      abs(inclusion_sum - observed_n) > 1e-5 ||
      any(!is.finite(intensity[[inclusion_name]])) ||
      any(intensity[[inclusion_name]] < 0 | intensity[[inclusion_name]] > 1)) {
    stop("Simulation QA failed for ", model_name, ".")
  }
  data.frame(
    null_model = model_name,
    simulation_count = nrow(model_draws),
    unique_simulation_ids = unique_draws,
    candidate_cell_count = candidate_n,
    fixed_observation_count = observed_n,
    inclusion_probability_sum = inclusion_sum,
    finite_road_and_town_statistics = finite_road && finite_town,
    stringsAsFactors = FALSE
  )
})
simulation_qa <- do.call(rbind, qa_rows)
write_csv(simulation_qa, file.path(tables_dir, "ghana_simulation_qa.csv"))

profile_summary <- profile_draws |>
  dplyr::group_by(null_model, statistic, distance_bin) |>
  dplyr::summarise(
    lower_2_5pct = stats::quantile(simulated_intensity, 0.025),
    median_intensity = stats::median(simulated_intensity),
    upper_97_5pct = stats::quantile(simulated_intensity, 0.975),
    observed_intensity = dplyr::first(observed_intensity),
    candidate_cell_count = dplyr::first(candidate_cell_count),
    observed_outside_pointwise_interval =
      dplyr::first(observed_intensity) < lower_2_5pct |
      dplyr::first(observed_intensity) > upper_97_5pct,
    .groups = "drop"
  )
write_csv(
  profile_summary,
  file.path(tables_dir, "ghana_distance_binned_intensity_summary.csv")
)
profile_summary$statistic <- factor(
  profile_summary$statistic,
  levels = c("nearest_major_road", "nearest_town_city"),
  labels = c("Major-road distance", "Town/city distance")
)
profile_summary$distance_bin <- factor(
  profile_summary$distance_bin,
  levels = c("0-5 km", "5-10 km", "10-20 km", "20-50 km", "50+ km")
)
profile_summary$null_model <- factor(
  profile_summary$null_model,
  levels = model_names,
  labels = c(
    "20 km cross-fit",
    "40 km cross-fit",
    "In-sample environmental",
    "Uniform region"
  )
)
profile_plot <- ggplot(
  profile_summary,
  aes(
    x = distance_bin,
    y = median_intensity,
    group = null_model,
    color = null_model,
    fill = null_model
  )
) +
  geom_ribbon(
    aes(ymin = lower_2_5pct, ymax = upper_97_5pct),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  geom_line(
    aes(y = observed_intensity, group = 1),
    color = "black",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(aes(y = observed_intensity), color = "black", size = 1.2) +
  facet_grid(statistic ~ null_model) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Observed and simulated distance-binned cocoa intensity",
    subtitle = paste(
      "Cell-positive share by distance band; shaded bands are pointwise 95% intervals",
      "from 999 fixed-count simulations. Dashed black line is observed."
    ),
    x = "Distance from access feature",
    y = "Positive-cell share"
  )
ggsave(
  file.path(figures_dir, "ghana_distance_binned_intensity.png"),
  plot = profile_plot,
  dpi = 300,
  width = 14,
  height = 8,
  units = "in"
)

draw_long <- rbind(
  data.frame(
    statistic = "Nearest major road",
    null_model = simulation_draws$null_model,
    distance_km = simulation_draws$median_major_road_distance_m / 1000
  ),
  data.frame(
    statistic = "Nearest town/city proxy",
    null_model = simulation_draws$null_model,
    distance_km = simulation_draws$median_town_city_distance_m / 1000
  )
)
draw_long$null_model <- factor(
  draw_long$null_model,
  levels = model_names,
  labels = c(
    "20 km cross-fit",
    "40 km cross-fit",
    "In-sample environmental",
    "Uniform region"
  )
)
observed_lines <- data.frame(
  statistic = c("Nearest major road", "Nearest town/city proxy"),
  distance_km = c(
    simulation_summary$observed_median_m[
      simulation_summary$statistic == "nearest_major_road"
    ][[1]] / 1000,
    simulation_summary$observed_median_m[
      simulation_summary$statistic == "nearest_town_city"
    ][[1]] / 1000
  )
)
distance_plot <- ggplot(draw_long, aes(x = distance_km, color = null_model)) +
  geom_density(linewidth = 0.8) +
  geom_vline(
    data = observed_lines,
    aes(xintercept = distance_km),
    color = "black",
    linetype = "dashed",
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~statistic, scales = "free_x", ncol = 1) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = "Observed and simulated median proximity",
    subtitle = "Dashed line is observed; each null uses 999 fixed-count draws.",
    x = "Median distance (km)",
    y = "Density",
    color = "Null model"
  )
ggsave(
  file.path(figures_dir, "ghana_simulation_distance_distributions.png"),
  plot = distance_plot,
  dpi = 300,
  width = 9,
  height = 8,
  units = "in"
)

observed_surface <- data.frame(
  x = intensity$x_m,
  y = intensity$y_m,
  surface = "Observed mapped-positive cells",
  value = intensity$observed_presence
)
primary_surface <- data.frame(
  x = intensity$x_m,
  y = intensity$y_m,
  surface = "20 km cross-fit mean simulated inclusion",
  value = intensity$spatial_block_20km_crossfit_simulated_inclusion
)
sensitivity_surface <- data.frame(
  x = intensity$x_m,
  y = intensity$y_m,
  surface = "40 km cross-fit mean simulated inclusion",
  value = intensity$spatial_block_40km_crossfit_simulated_inclusion
)
surface_plot_data <- rbind(
  observed_surface,
  primary_surface,
  sensitivity_surface
)
surface_plot_data$surface <- factor(
  surface_plot_data$surface,
  levels = c(
    "Observed mapped-positive cells",
    "20 km cross-fit mean simulated inclusion",
    "40 km cross-fit mean simulated inclusion"
  )
)
intensity_plot <- ggplot(
  surface_plot_data,
  aes(x = x, y = y, fill = value)
) +
  geom_raster() +
  facet_wrap(~surface, ncol = 1) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08519C",
    limits = c(0, 1),
    name = "Cell value"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Observed cocoa occupancy and simulated first-order intensity",
    subtitle = "Mean cell inclusion across 999 fixed-count simulations; 1 km grid.",
    x = "Easting (m)",
    y = "Northing (m)"
  )
ggsave(
  file.path(figures_dir, "ghana_environmental_intensity_comparison.png"),
  plot = intensity_plot,
  dpi = 300,
  width = 9,
  height = 14,
  units = "in"
)

intensity$residual_20km <- intensity$observed_presence -
  intensity$spatial_block_20km_crossfit_simulated_inclusion
residual_plot <- ggplot(
  intensity,
  aes(x = x_m, y = y_m, fill = residual_20km)
) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Observed -\nmean inclusion"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Observed occupancy minus 20 km environmental-null intensity",
    subtitle = "Positive values mark cells mapped positive but rarely selected by the null.",
    x = "Easting (m)",
    y = "Northing (m)"
  )
ggsave(
  file.path(figures_dir, "ghana_environmental_intensity_residual.png"),
  plot = residual_plot,
  dpi = 300,
  width = 9,
  height = 8,
  units = "in"
)

message("National simulation QA and figures complete.")
print(crossfit_performance)
print(simulation_qa)
