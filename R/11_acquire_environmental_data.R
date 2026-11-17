# Acquire WorldClim 2.1 climate and elevation predictors for the environmental
# occurrence model. Only BIO1, BIO12, and elevation are extracted.

source("R/01_setup.R")

worldclim <- list(
  bio = list(
    url = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip",
    archive = file.path(raw_dir, "wc2.1_2.5m_bio.zip"),
    expected_size = 658405521,
    members = c("wc2.1_2.5m_bio_1.tif", "wc2.1_2.5m_bio_12.tif")
  ),
  elevation = list(
    url = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_elev.zip",
    archive = file.path(raw_dir, "wc2.1_2.5m_elev.zip"),
    expected_size = 17566799,
    members = "wc2.1_2.5m_elev.tif"
  )
)

for (dataset in worldclim) {
  missing_members <- dataset$members[
    !file.exists(file.path(raw_dir, dataset$members))
  ]
  if (length(missing_members) > 0 && !file.exists(dataset$archive)) {
    status <- utils::download.file(
      dataset$url,
      destfile = dataset$archive,
      mode = "wb",
      quiet = FALSE
    )
    if (status != 0) {
      stop("WorldClim download failed: ", dataset$url)
    }
  }
  if (file.exists(dataset$archive)) {
    observed_size <- unname(file.info(dataset$archive)$size)
    if (observed_size != dataset$expected_size) {
      stop(
        "Unexpected WorldClim archive size for ", basename(dataset$archive),
        ": expected ", dataset$expected_size, " bytes, got ", observed_size, "."
      )
    }
  }
  if (length(missing_members) == 0) {
    next
  }
  if (!file.exists(dataset$archive)) {
    stop("Missing WorldClim archive needed to extract: ", paste(missing_members, collapse = ", "))
  }
  members <- utils::unzip(dataset$archive, list = TRUE)$Name
  for (target in missing_members) {
    matched_member <- members[basename(members) == target]
    if (length(matched_member) != 1L) {
      stop(
        "Expected exactly one ", target, " member in ",
        basename(dataset$archive), "; found ", length(matched_member), "."
      )
    }
    extracted_path <- file.path(raw_dir, target)
    if (!file.exists(extracted_path)) {
      utils::unzip(
        dataset$archive,
        files = matched_member,
        exdir = raw_dir,
        overwrite = FALSE
      )
      extracted_member_path <- file.path(raw_dir, matched_member)
      if (!identical(extracted_member_path, extracted_path) &&
          !file.rename(extracted_member_path, extracted_path)) {
        stop("Could not move extracted WorldClim raster to ", extracted_path)
      }
    }
    if (!file.exists(extracted_path) || file.info(extracted_path)$size <= 0) {
      stop("WorldClim raster extraction failed: ", target)
    }
  }
}

rasters <- list(
  bio1 = list(
    path = file.path(raw_dir, "wc2.1_2.5m_bio_1.tif"),
    variable = "BIO1 annual mean temperature (stored in tenths of degrees C)",
    data_year = "1970-2000"
  ),
  bio12 = list(
    path = file.path(raw_dir, "wc2.1_2.5m_bio_12.tif"),
    variable = "BIO12 annual precipitation (mm)",
    data_year = "1970-2000"
  ),
  elevation = list(
    path = file.path(raw_dir, "wc2.1_2.5m_elev.tif"),
    variable = "Elevation (m; SRTM-derived)",
    data_year = "SRTM source"
  )
)

sha256_file <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (.Platform$OS.type == "windows") {
    output <- system2(
      "certutil.exe",
      args = c(
        "-hashfile",
        shQuote(normalizePath(path, winslash = "\\"), type = "cmd"),
        "SHA256"
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  } else {
    output <- system2(
      "sha256sum",
      args = shQuote(path),
      stdout = TRUE,
      stderr = TRUE
    )
  }
  hashes <- regmatches(
    output,
    gregexpr("[[:xdigit:]]{64}", output, perl = TRUE, useBytes = TRUE)
  )
  hashes <- unlist(hashes, use.names = FALSE)
  if (length(hashes) != 1L) {
    stop("Could not compute a SHA-256 checksum for ", path, ".")
  }
  tolower(hashes[[1]])
}

manifest_rows <- lapply(names(rasters), function(variable_name) {
  item <- rasters[[variable_name]]
  raster <- terra::rast(item$path)
  if (terra::nlyr(raster) != 1L || is.na(terra::crs(raster))) {
    stop("WorldClim raster is missing its expected single band or CRS: ", item$path)
  }
  data.frame(
    source = paste0("WORLDCLIM_21_", toupper(variable_name)),
    provider = "WorldClim / University of California, Davis",
    access_date = as.character(Sys.Date()),
    local_path = file.path("data", "raw", basename(item$path)),
    license = "WorldClim data use terms; cite WorldClim 2.1",
    source_url = if (variable_name == "elevation") worldclim$elevation$url else worldclim$bio$url,
    file_size_bytes = unname(file.info(item$path)$size),
    sha256 = sha256_file(item$path),
    crs = terra::crs(raster),
    bounding_box = paste(
      terra::xmin(raster), terra::ymin(raster),
      terra::xmax(raster), terra::ymax(raster),
      sep = ";"
    ),
    dimensions_or_feature_count = paste0(terra::ncol(raster), "x", terra::nrow(raster)),
    fields_or_band_info = item$variable,
    nodata_or_missing = "See GeoTIFF NoData metadata",
    coverage = "Global land coverage; clipped to the Ghana study area during preparation",
    data_year = item$data_year,
    caveats = "WorldClim 2.1; BIO1 and BIO12 are approximately 4.5 km source resolution and are resampled to the 1 km analysis grid.",
    status = "downloaded_validated",
    stringsAsFactors = FALSE
  )
})
archive_names <- c("bio", "elevation")
archive_names <- archive_names[
  vapply(
    worldclim[archive_names],
    function(dataset) file.exists(dataset$archive),
    logical(1)
  )
]
archive_rows <- lapply(archive_names, function(dataset_name) {
  dataset <- worldclim[[dataset_name]]
  data.frame(
    source = paste0("WORLDCLIM_21_", toupper(dataset_name), "_ARCHIVE"),
    provider = "WorldClim / University of California, Davis",
    access_date = as.character(Sys.Date()),
    local_path = file.path("data", "raw", basename(dataset$archive)),
    license = "WorldClim data use terms; cite WorldClim 2.1",
    source_url = dataset$url,
    file_size_bytes = unname(file.info(dataset$archive)$size),
    sha256 = sha256_file(dataset$archive),
    crs = NA_character_,
    bounding_box = NA_character_,
    dimensions_or_feature_count = "ZIP archive",
    fields_or_band_info = paste("Contains selected rasters:", paste(dataset$members, collapse = ", ")),
    nodata_or_missing = NA_character_,
    coverage = "Global",
    data_year = "1970-2000 climate normals or SRTM-derived elevation",
    caveats = "Official WorldClim 2.1 global archive; only selected members are extracted for analysis.",
    status = "downloaded_validated",
    stringsAsFactors = FALSE
  )
})

manifest_path <- file.path(project_root, "data", "manifest.csv")
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
new_rows <- rbind(do.call(rbind, manifest_rows), do.call(rbind, archive_rows))
manifest <- manifest[!grepl("^WORLDCLIM_21_", manifest$source), , drop = FALSE]
readr::write_csv(rbind(manifest, new_rows), manifest_path)

message("WorldClim BIO1, BIO12, and elevation rasters acquired and validated.")
