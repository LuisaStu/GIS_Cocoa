# Build a national Ghana 1 km occupancy grid and spatially cross-fitted
# environmental models. Road and town access are excluded from predictors.

source("R/01_setup.R")

national_crs <- 2136
gold_coast_foot_m <- 0.304799710181509
grid_resolution_m <- 1000
block_sizes_m <- c(20000, 40000)
number_of_folds <- 5L
fold_seed <- 20261117

osm_path <- file.path(raw_dir, "ghana.gpkg")
cocoa_path <- file.path(raw_dir, "cocoa_map_th_065.tif")
climate_paths <- c(
  bio1 = file.path(raw_dir, "wc2.1_2.5m_bio_1.tif"),
  bio12 = file.path(raw_dir, "wc2.1_2.5m_bio_12.tif"),
  elevation = file.path(raw_dir, "wc2.1_2.5m_elev.tif")
)
required_paths <- c(osm_path, cocoa_path, climate_paths)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing Ghana-wide inputs: ", paste(missing_paths, collapse = ", "))
}
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("The `mgcv` package is required for environmental cross-fitting.")
}

admin <- st_read(
  osm_path,
  layer = "gis_osm_adminareas_a_free",
  quiet = TRUE
)
ghana <- admin[admin$fclass == "national" & admin$name == "Ghana", , drop = FALSE]
if (nrow(ghana) != 1L) {
  stop("Expected one Ghana national boundary in the Geofabrik administrative layer.")
}
ghana <- st_make_valid(st_transform(ghana, national_crs))
ghana_union <- st_union(ghana)
ghana_area_km2 <- as.numeric(units::set_units(st_area(ghana_union), "km^2"))
if (!is.finite(ghana_area_km2) || ghana_area_km2 < 200000 || ghana_area_km2 > 280000) {
  stop("The Ghana boundary has an implausible area: ", ghana_area_km2, " km^2.")
}

candidate_path <- file.path(processed_dir, "ghana_1km_candidate_points.gpkg")
cocoa_occupancy_path <- file.path(processed_dir, "ghana_1km_cocoa_occupancy.tif")
grid_cache_path <- file.path(processed_dir, "ghana_1km_model_grid.rds")
if (file.exists(candidate_path)) {
  candidate_points <- st_read(candidate_path, quiet = TRUE)
  if (!isTRUE(st_crs(candidate_points) == st_crs(national_crs))) {
    stop("Cached Ghana candidate grid has the wrong CRS.")
  }
  candidate_data <- st_drop_geometry(candidate_points)
  required_cached <- c(
    "cell", "x", "y", "cocoa_fraction", "bio1_tenths_c",
    "bio12_mm", "elevation_m", "nearest_major_road_distance_m",
    "nearest_town_distance_m"
  )
  if (!all(required_cached %in% names(candidate_data))) {
    stop("Cached Ghana candidate grid is incomplete; remove it only after inspection.")
  }
  message("Reusing the cached Ghana-wide candidate grid.")
} else if (file.exists(grid_cache_path)) {
  candidate_data <- readRDS(grid_cache_path)
  message("Reusing the rasterized Ghana model grid.")
} else {
  region_vect <- terra::vect(ghana_union)
  region_wgs84 <- st_transform(ghana_union, 4326)
  region_wgs84_vect <- terra::vect(region_wgs84)
  resolution_units <- grid_resolution_m / gold_coast_foot_m
  template <- terra::rast(
    terra::ext(region_vect),
    resolution = resolution_units,
    crs = paste0("EPSG:", national_crs)
  )
  region_mask <- terra::rasterize(region_vect, template, field = 1)
  names(region_mask) <- "in_ghana"

  cocoa_source <- terra::rast(cocoa_path)
  if (file.exists(cocoa_occupancy_path)) {
    cocoa_occupancy <- terra::rast(cocoa_occupancy_path)
  } else {
    cocoa_occupancy <- terra::project(
      cocoa_source,
      template,
      method = "mean",
      filename = cocoa_occupancy_path,
      overwrite = TRUE,
      wopt = list(gdal = "COMPRESS=LZW")
    )
  }
  if (!terra::compareGeom(cocoa_occupancy, template, stopOnError = FALSE)) {
    stop("National cocoa raster did not align to the Ghana analysis grid.")
  }
  source_extent <- terra::ext(cocoa_source)
  source_extent_values <- as.vector(source_extent)
  source_bbox <- st_bbox(c(
    xmin = source_extent_values[["xmin"]],
    ymin = source_extent_values[["ymin"]],
    xmax = source_extent_values[["xmax"]],
    ymax = source_extent_values[["ymax"]]
  ), crs = st_crs(4326))
  source_coverage <- st_transform(st_as_sfc(source_bbox), national_crs)
  rm(cocoa_source)
  invisible(gc())
  coverage_mask <- terra::rasterize(
    terra::vect(source_coverage),
    template,
    field = 1,
    touches = TRUE
  )
  names(coverage_mask) <- "cocoa_coverage"
  cocoa_present <- terra::ifel(is.na(cocoa_occupancy), 0, cocoa_occupancy > 0)
  names(cocoa_present) <- "cocoa_positive"
  rm(cocoa_occupancy)
  invisible(gc())

  project_climate <- function(path, layer_name, method) {
    cache_path <- file.path(
      processed_dir,
      paste0("ghana_1km_", layer_name, ".tif")
    )
    if (file.exists(cache_path)) {
      projected <- terra::rast(cache_path)
      if (!terra::compareGeom(projected, template, stopOnError = FALSE)) {
        stop("Cached Ghana predictor has the wrong grid: ", cache_path)
      }
      names(projected) <- layer_name
      return(projected)
    }
    source_raster <- terra::rast(path)
    cropped <- terra::crop(
      source_raster,
      terra::ext(region_wgs84_vect),
      snap = "out"
    )
    if (terra::ncell(cropped) == 0) {
      stop("Climate/elevation raster does not overlap Ghana: ", path)
    }
    projected <- terra::project(
      cropped,
      template,
      method = method,
      filename = cache_path,
      overwrite = TRUE,
      wopt = list(gdal = "COMPRESS=LZW")
    )
    if (!terra::compareGeom(projected, template, stopOnError = FALSE)) {
      stop("Projected predictor does not align with the Ghana grid: ", path)
    }
    names(projected) <- layer_name
    rm(source_raster, cropped)
    invisible(gc())
    projected
  }

  bio1 <- project_climate(climate_paths[["bio1"]], "bio1_tenths_c", "bilinear")
  message("Prepared national BIO1 layer.")
  bio12 <- project_climate(climate_paths[["bio12"]], "bio12_mm", "bilinear")
  message("Prepared national BIO12 layer.")
  elevation <- project_climate(
    climate_paths[["elevation"]],
    "elevation_m",
    "bilinear"
  )
  message("Prepared national elevation layer.")
  grid_stack <- c(region_mask, coverage_mask, cocoa_present, bio1, bio12, elevation)
  grid_data <- terra::as.data.frame(
    grid_stack,
    xy = TRUE,
    cells = TRUE,
    na.rm = FALSE
  )
  expected_fields <- c(
    "cell", "x", "y", "in_ghana", "cocoa_coverage", "cocoa_positive",
    "bio1_tenths_c", "bio12_mm", "elevation_m"
  )
  if (!all(expected_fields %in% names(grid_data))) {
    stop(
      "National raster stack is missing expected fields: ",
      paste(setdiff(expected_fields, names(grid_data)), collapse = ", ")
    )
  }
  candidate_data <- grid_data[
    !is.na(grid_data$in_ghana) &
      !is.na(grid_data$cocoa_coverage) &
      stats::complete.cases(grid_data[c(
        "bio1_tenths_c", "bio12_mm", "elevation_m", "cocoa_positive"
      )]),
    ,
    drop = FALSE
  ]
  if (nrow(candidate_data) == 0) {
    stop("No Ghana grid cells have complete cocoa coverage and environmental data.")
  }
  saveRDS(candidate_data, grid_cache_path)
  message("Saved Ghana model-grid checkpoint with ", nrow(candidate_data), " cells.")
}

if (any(!is.finite(candidate_data$bio1_tenths_c)) ||
    any(!is.finite(candidate_data$bio12_mm)) ||
    any(!is.finite(candidate_data$elevation_m))) {
  stop("Ghana candidate environmental predictors contain non-finite values.")
}
candidate_data$mean_temp_c <- candidate_data$bio1_tenths_c / 10
candidate_data$log_precip <- log1p(candidate_data$bio12_mm)
candidate_data$cocoa_positive <- as.logical(candidate_data$cocoa_positive)
candidate_data$cocoa_fraction <- as.numeric(candidate_data$cocoa_positive)
candidate_data$presence_response <- as.integer(candidate_data$cocoa_positive)
if (length(unique(candidate_data$presence_response)) != 2L) {
  stop("The Ghana candidate grid must contain mapped positives and negatives.")
}
observed_n <- sum(candidate_data$cocoa_positive)
if (observed_n < 20 || observed_n >= nrow(candidate_data)) {
  stop("The Ghana-wide mapped-positive count is not suitable for fixed-count sampling.")
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

fit_spatial_crossfit <- function(data, block_size_m, seed) {
  message("Starting ", block_size_m, " m spatial cross-fitting.")
  block_size_units <- block_size_m / gold_coast_foot_m
  data$block_x <- floor(data$x / block_size_units)
  data$block_y <- floor(data$y / block_size_units)
  blocks <- unique(data[c("block_x", "block_y")])
  if (nrow(blocks) < number_of_folds) {
    stop("Too few spatial blocks for the ", block_size_m, " m cross-fit.")
  }
  set.seed(seed)
  clusters <- stats::kmeans(
    as.matrix(blocks[c("block_x", "block_y")]),
    centers = number_of_folds,
    nstart = 50
  )
  blocks$spatial_fold <- clusters$cluster
  data$spatial_fold <- blocks$spatial_fold[
    match(
      paste(data$block_x, data$block_y),
      paste(blocks$block_x, blocks$block_y)
    )
  ]
  predictions <- rep(NA_real_, nrow(data))
  fold_diagnostics <- vector("list", number_of_folds)
  fold_terms <- vector("list", number_of_folds)

  for (fold in seq_len(number_of_folds)) {
    training <- data[data$spatial_fold != fold, , drop = FALSE]
    testing_rows <- data$spatial_fold == fold
    testing <- data[testing_rows, , drop = FALSE]
    if (nrow(testing) == 0 ||
        length(unique(training$presence_response)) != 2L) {
      stop(
        "Spatial fold ", fold, " at ", block_size_m,
        " m has an empty test group or training data without both classes."
      )
    }
    message(
      "Fitting ", block_size_m, " m fold ", fold, " of ",
      number_of_folds, "."
    )
    model <- mgcv::bam(
      presence_response ~
        s(mean_temp_c, k = 5) +
        s(log_precip, k = 5) +
        s(elevation_m, k = 5),
      data = training,
      family = stats::binomial(),
      method = "fREML",
      select = TRUE,
      discrete = TRUE,
      nthreads = 2
    )
    fold_prediction <- as.numeric(stats::predict(
      model,
      newdata = testing,
      type = "response"
    ))
    if (any(!is.finite(fold_prediction))) {
      stop("Non-finite held-out predictions at ", block_size_m, " m, fold ", fold, ".")
    }
    predictions[testing_rows] <- fold_prediction
    clipped <- pmin(pmax(fold_prediction, 1e-12), 1 - 1e-12)
    heldout <- testing$presence_response
    fold_diagnostics[[fold]] <- data.frame(
      block_size_m = block_size_m,
      spatial_fold = fold,
      spatial_blocks = sum(blocks$spatial_fold == fold),
      training_cells = nrow(training),
      heldout_cells = nrow(testing),
      heldout_positive_cells = sum(heldout == 1L),
      heldout_negative_cells = sum(heldout == 0L),
      heldout_auc = auc_score(heldout, fold_prediction),
      heldout_brier_score = mean((fold_prediction - heldout)^2),
      heldout_log_loss = -mean(
        heldout * log(clipped) + (1 - heldout) * log(1 - clipped)
      ),
      stringsAsFactors = FALSE
    )
    terms <- as.data.frame(summary(model)$s.table)
    terms$term <- rownames(terms)
    terms$block_size_m <- block_size_m
    terms$spatial_fold <- fold
    rownames(terms) <- NULL
    fold_terms[[fold]] <- terms
    rm(model, training, testing)
    invisible(gc())
  }
  if (any(!is.finite(predictions))) {
    stop("The ", block_size_m, " m cross-fit did not predict every candidate cell.")
  }
  list(
    predictions = predictions,
    diagnostics = do.call(rbind, fold_diagnostics),
    terms = do.call(rbind, fold_terms)
  )
}

crossfit_20km <- fit_spatial_crossfit(candidate_data, block_sizes_m[[1]], fold_seed)
crossfit_40km <- fit_spatial_crossfit(candidate_data, block_sizes_m[[2]], fold_seed)
candidate_data$environmental_score <- pmax(
  crossfit_20km$predictions,
  .Machine$double.eps
)
candidate_data$environmental_weight <- candidate_data$environmental_score /
  sum(candidate_data$environmental_score)
candidate_data$environmental_score_40km <- pmax(
  crossfit_40km$predictions,
  .Machine$double.eps
)
candidate_data$environmental_weight_40km <-
  candidate_data$environmental_score_40km /
  sum(candidate_data$environmental_score_40km)
if (any(!is.finite(candidate_data$environmental_weight)) ||
    any(candidate_data$environmental_weight <= 0) ||
    abs(sum(candidate_data$environmental_weight) - 1) > 1e-8 ||
    any(!is.finite(candidate_data$environmental_weight_40km)) ||
    abs(sum(candidate_data$environmental_weight_40km) - 1) > 1e-8) {
  stop("Ghana environmental simulation weights are invalid.")
}

roads <- st_read(osm_path, layer = "gis_osm_roads_free", quiet = TRUE)
message("Clipping national major-road network.")
major_classes <- c("motorway", "trunk", "primary", "secondary")
roads <- roads[roads$fclass %in% major_classes, , drop = FALSE]
roads <- st_transform(roads, national_crs)
roads <- st_intersection(roads, ghana_union)
roads <- st_collection_extract(roads, "LINESTRING", warn = FALSE)
if (nrow(roads) == 0) {
  stop("No major roads intersect the Ghana boundary.")
}
places <- st_read(osm_path, layer = "gis_osm_places_free", quiet = TRUE)
urban_classes <- c("city", "town", "national_capital")
places <- places[places$fclass %in% urban_classes, , drop = FALSE]
places <- st_transform(places, national_crs)
places <- st_filter(places, ghana_union, .predicate = st_within)
if (nrow(places) == 0) {
  stop("No OSM towns or cities intersect the Ghana boundary.")
}

candidate_points <- st_as_sf(
  candidate_data,
  coords = c("x", "y"),
  remove = FALSE,
  crs = national_crs
)
message("Computing national nearest-major-road distances.")
road_nearest <- st_nearest_feature(candidate_points, roads)
candidate_points$nearest_major_road_distance_m <- as.numeric(units::set_units(
  st_distance(candidate_points, roads[road_nearest, ], by_element = TRUE),
  "m"
))
message("Computing national nearest-town distances.")
town_nearest <- st_nearest_feature(candidate_points, places)
candidate_points$nearest_town_distance_m <- as.numeric(units::set_units(
  st_distance(candidate_points, places[town_nearest, ], by_element = TRUE),
  "m"
))
if (any(!is.finite(candidate_points$nearest_major_road_distance_m)) ||
    any(!is.finite(candidate_points$nearest_town_distance_m))) {
  stop("Ghana-wide road or town distances are incomplete.")
}

candidate_path <- file.path(processed_dir, "ghana_1km_candidate_points.gpkg")
st_write(candidate_points, candidate_path, delete_dsn = TRUE, quiet = TRUE)
candidate_check <- st_read(candidate_path, quiet = TRUE)
if (nrow(candidate_check) != nrow(candidate_points) ||
    sum(candidate_check$cocoa_positive) != observed_n ||
    any(!is.finite(candidate_check$environmental_weight)) ||
    any(!is.finite(candidate_check$nearest_major_road_distance_m)) ||
    any(!is.finite(candidate_check$nearest_town_distance_m))) {
  stop("Ghana candidate GeoPackage read-back validation failed.")
}

diagnostics <- rbind(crossfit_20km$diagnostics, crossfit_40km$diagnostics)
terms <- rbind(crossfit_20km$terms, crossfit_40km$terms)
readr::write_csv(
  diagnostics,
  file.path(tables_dir, "ghana_environmental_cv_diagnostics.csv")
)
readr::write_csv(
  terms,
  file.path(tables_dir, "ghana_environmental_model_terms.csv")
)
summary <- data.frame(
  geographic_scope = "Ghana",
  boundary_source = "Geofabrik Ghana OSM national administrative polygon",
  analysis_crs = "EPSG:2136 (Accra / Ghana National Grid; Gold Coast feet)",
  grid_resolution_m = grid_resolution_m,
  boundary_area_km2 = ghana_area_km2,
  eligible_candidate_cells = nrow(candidate_data),
  observed_positive_cells = observed_n,
  mapped_positive_share = observed_n / nrow(candidate_data),
  cocoa_source_coverage_extent = paste(
    signif(c(
      terra::ext(terra::rast(cocoa_path))$xmin,
      terra::ext(terra::rast(cocoa_path))$xmax
    ), 6),
    collapse = " to "
  ),
  major_road_features = nrow(roads),
  urban_centre_features = nrow(places),
  crossfit_20km_mean_heldout_auc = mean(
    crossfit_20km$diagnostics$heldout_auc,
    na.rm = TRUE
  ),
  crossfit_20km_mean_heldout_brier = mean(
    crossfit_20km$diagnostics$heldout_brier_score
  ),
  crossfit_20km_mean_heldout_log_loss = mean(
    crossfit_20km$diagnostics$heldout_log_loss
  ),
  crossfit_40km_mean_heldout_auc = mean(
    crossfit_40km$diagnostics$heldout_auc,
    na.rm = TRUE
  ),
  crossfit_40km_mean_heldout_brier = mean(
    crossfit_40km$diagnostics$heldout_brier_score
  ),
  crossfit_40km_mean_heldout_log_loss = mean(
    crossfit_40km$diagnostics$heldout_log_loss
  ),
  observed_median_major_road_distance_m = median(
    candidate_points$nearest_major_road_distance_m[candidate_points$cocoa_positive]
  ),
  observed_median_town_distance_m = median(
    candidate_points$nearest_town_distance_m[candidate_points$cocoa_positive]
  ),
  nodata_treatment = "Source NaN is interpreted as no cocoa per researcher instruction; cells outside the raster extent are excluded.",
  stringsAsFactors = FALSE
)
readr::write_csv(summary, file.path(tables_dir, "ghana_null_model_summary.csv"))
message("Ghana-wide shared grid and spatial cross-fits complete.")
message("Eligible 1 km cells: ", nrow(candidate_data))
message("Mapped-positive cells: ", observed_n)
message("20 km held-out AUC: ", round(summary$crossfit_20km_mean_heldout_auc, 3))
message("40 km held-out AUC: ", round(summary$crossfit_40km_mean_heldout_auc, 3))
message("Median distance to major roads (m): ", round(summary$observed_median_major_road_distance_m, 1))
message("Median distance to towns/cities (m): ", round(summary$observed_median_town_distance_m, 1))
