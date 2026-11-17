# Cocoa Farming and Access to Markets in Ghana

This Georgetown University spatial data science project asks whether mapped cocoa-growing areas in Ghana are closer to major roads and urban centres than expected given their climate and elevation. Road proximity is the primary measure of accessibility. Urban-centre proximity is a secondary proxy for potential market access. The study examines first-order spatial intensity, meaning how the expected occurrence of mapped cocoa varies across locations.

## Data and sources

| Input | Role in the analysis | Source |
| --- | --- | --- |
| Thresholded cocoa map | Approximately 10 m satellite-derived cocoa extent for Ghana and C&ocirc;te d'Ivoire, based on imagery from around 2019-2021. The workflow uses the 0.65-threshold raster. | [ETH Zurich Research Collection](https://doi.org/10.3929/ethz-b-000654400), [Kalischek et al.](https://doi.org/10.1038/s43016-023-00751-8) |
| Roads, settlements and boundary | The November 9, 2026 OpenStreetMap snapshot provides motorways and trunk, primary and secondary roads, plus mapped cities, towns and national capitals. The national administrative polygon defines the study area. | [Geofabrik Ghana](https://download.geofabrik.de/africa/ghana.html), [OpenStreetMap contributors and ODbL](https://www.openstreetmap.org/copyright) |
| Climate and elevation | Annual mean temperature (BIO1), annual precipitation (BIO12) and SRTM-derived elevation at 2.5 arc-minutes, roughly 4.5 km in Ghana. Climate normals cover 1970-2000. | [WorldClim 2.1](https://www.worldclim.org/data/worldclim21.html) |

The [data manifest](data/manifest.csv) records input filenames, versions, licences and checksums.

## Methodology

1. **Align the layers on a 1 km grid.** The workflow uses the Ghana National Grid (EPSG:2136) and converts projection units to metres for reporting. Eligible cells fall inside the working Ghana boundary and the cocoa raster coverage mask and have complete environmental predictors. Each cell with mapped cocoa contributes one observation at its centre. The workflow treats NaN within the raster coverage mask as no mapped cocoa, an assumption awaiting provider confirmation.
2. **Measure proximity.** Compute the straight-line distance from each eligible cell centre to the nearest major road and the nearest qualifying urban centre. The primary statistic is median road distance across cocoa-positive cells. Median urban-centre distance is secondary.
3. **Estimate an environmental baseline.** Fit a binomial generalized additive model with smooth terms for mean temperature, log-transformed annual precipitation and elevation. The primary model groups 20 km spatial blocks into five folds and predicts each fold using a model fitted on the other four. The resulting scores weight how often cells are selected for simulations. Road and settlement distances are excluded from the model.
4. **Simulate and compare.** Generate 999 patterns per baseline, drawing cells without replacement and holding the observed cocoa-positive count fixed. Compare the primary baseline with a 40 km block design, an environmental model fitted to all eligible cells, and uniform sampling over eligible cells. The corrected lower-tail Monte Carlo p-value is (1 + simulated medians at or below the observed median) / 1000. Maps compare observed occupancy with simulated inclusion frequency. An exploratory distance profile compares observed and simulated shares across distance bands with pointwise 95% simulation intervals.

## Interpretation

The national run includes 240,294 eligible cells, of which 57,131 contain mapped cocoa. Median major-road distance was 11.96 km for observed cocoa-positive cells and 13.25 km across primary simulated patterns. The corresponding urban-centre distances were 26.43 km and 28.58 km. Both observed medians fell below all 999 primary simulated medians, yielding a corrected Monte Carlo p-value of 0.001 for each comparison. See the [Quarto manuscript](manuscript/manuscript.qmd) for the full comparisons and model diagnostics.

The results describe spatial association conditional on the mapped cocoa extent, study window and environmental model. A positive cell represents mapped cocoa presence and can contain multiple farms. Travel times and access to verified cocoa purchasing sites would require dedicated data. The OSM snapshot postdates the cocoa imagery, the working boundary needs independent validation, and the study window lacks an independent land-cover mask. Climate layers are coarser than the analysis grid, and each block design has a held-out fold without positive cells. These limitations and possible cocoa-map detection bias affect the strength of the comparison.

## Reproduce the analysis

Install R and Quarto. From the repository root, restore the R environment and acquire the raster inputs.

```sh
Rscript -e "renv::restore()"
Rscript R/02_acquire_data.R
Rscript R/11_acquire_environmental_data.R
```

Extract the [dated Geofabrik archive](https://download.geofabrik.de/africa/ghana-261109-free.gpkg.zip) into data/raw/ and verify the SHA-256 of ghana.gpkg against the manifest. Then run the national workflow.

```sh
Rscript R/12_ghana_wide_null_model.R
Rscript R/13_ghana_wide_simulations.R
Rscript R/14_ghana_wide_reporting.R
quarto render manuscript --to html
```

Intermediate spatial data go to data/processed/. Figures, tables and the rendered manuscript go to outputs/.
