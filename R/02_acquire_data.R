# 02_acquire_data.R
# Download the ETH Zurich binary threshold raster in verified HTTP byte ranges,
# then validate the local file and write its observed metadata. Byte ranges are
# used because the provider rate-limited direct full-file requests.

source("R/01_setup.R")

cocoa_path <- file.path(raw_dir, "cocoa_map_th_065.tif")
expected_size_bytes <- 119543698
download_url <- paste0(
  "https://www.research-collection.ethz.ch/server/api/core/bitstreams/",
  "2956d638-81e9-4075-837c-5d2fd4820faf/content"
)

download_eth_raster <- function() {
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Install the `curl` R package to download the ETH raster.")
  }

  partial_path <- paste0(cocoa_path, ".part")
  chunk_size <- 8 * 1024^2
  if (!file.exists(partial_path) && !file.create(partial_path)) {
    stop("Could not create the partial ETH download file.")
  }
  offset <- if (file.exists(partial_path)) {
    unname(file.info(partial_path)$size)
  } else {
    0
  }
  if (offset > expected_size_bytes ||
      (offset < expected_size_bytes && offset %% chunk_size != 0)) {
    stop("Partial download has an unexpected length; inspect it before resuming.")
  }

  output <- file(partial_path, open = if (offset == 0) "wb" else "ab")
  output_open <- TRUE
  on.exit(if (output_open) close(output), add = TRUE)
  if (offset < expected_size_bytes) {
    starts <- seq.int(from = offset, to = expected_size_bytes - 1, by = chunk_size)
    for (start in starts) {
      end <- min(start + chunk_size - 1, expected_size_bytes - 1)
      chunk_path <- tempfile(fileext = ".part")
      response <- curl::curl_fetch_disk(
        download_url,
        path = chunk_path,
        handle = curl::new_handle(
          range = paste0(start, "-", end),
          useragent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
        )
      )
      expected_chunk_size <- end - start + 1
      observed_chunk_size <- unname(file.info(chunk_path)$size)
      if (response$status_code != 206 ||
          observed_chunk_size != expected_chunk_size) {
        unlink(chunk_path)
        stop(
          "ETH byte-range request failed for bytes ", start, "-", end,
          " (HTTP ", response$status_code, "; received ",
          observed_chunk_size, " bytes)."
        )
      }

      input <- file(chunk_path, open = "rb")
      bytes <- readBin(input, what = "raw", n = expected_chunk_size)
      close(input)
      writeBin(bytes, output)
      unlink(chunk_path)
      message("Downloaded ", end + 1, " / ", expected_size_bytes, " bytes.")
      Sys.sleep(0.5)
    }
  }
  close(output)
  output_open <- FALSE
  if (unname(file.info(partial_path)$size) != expected_size_bytes) {
    stop("The completed ETH download has an unexpected size.")
  }
  if (!file.rename(partial_path, cocoa_path)) {
    stop("Could not move the validated partial ETH download into place.")
  }
  invisible(cocoa_path)
}

if (!file.exists(cocoa_path)) {
  download_eth_raster()
}

observed_size_bytes <- unname(file.info(cocoa_path)$size)
if (observed_size_bytes != expected_size_bytes) {
  stop(
    "Unexpected ETH raster size: expected ", expected_size_bytes,
    " bytes, got ", observed_size_bytes, "."
  )
}

cocoa <- terra::rast(cocoa_path)
if (terra::nlyr(cocoa) != 1L) {
  stop("Expected a single-band ETH cocoa raster.")
}
if (!terra::same.crs(cocoa, "EPSG:4326")) {
  stop("Unexpected native CRS for ETH raster: ", terra::crs(cocoa))
}

metadata <- data.frame(
  provider = "ETH Zurich",
  title = "Côte d'Ivoire and Ghana cocoa extent derived from Sentinel-2",
  access_date = as.character(Sys.Date()),
  local_path = "data/raw/cocoa_map_th_065.tif",
  license = "CC BY 4.0",
  url = download_url,
  file_size_bytes = observed_size_bytes,
  raster_crs = "EPSG:4326",
  raster_width = terra::ncol(cocoa),
  raster_height = terra::nrow(cocoa),
  raster_bands = terra::nlyr(cocoa),
  raster_datatype = terra::datatype(cocoa),
  raster_nodata = "NaN",
  xmin = terra::xmin(cocoa),
  ymin = terra::ymin(cocoa),
  xmax = terra::xmax(cocoa),
  ymax = terra::ymax(cocoa),
  stringsAsFactors = FALSE
)

readr::write_csv(metadata, file.path(tables_dir, "cocoa_raster_metadata.csv"))
message("ETH cocoa raster metadata validated and saved.")
message("Native CRS: ", terra::crs(cocoa))
message(
  "Dimensions: ", terra::ncol(cocoa), " x ", terra::nrow(cocoa),
  "; file size: ", observed_size_bytes, " bytes."
)
