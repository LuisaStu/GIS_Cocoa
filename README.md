# Cocoa Farming and Access to Markets in Ghana

This Georgetown University spatial data science project asks whether mapped cocoa-growing areas in Ghana are closer to major roads and urban centres than expected given their climate and elevation. Road proximity is the primary measure of accessibility. Urban-centre proximity is a secondary proxy for potential market access. The study examines first-order spatial intensity, meaning how the expected occurrence of mapped cocoa varies across locations.

The directional hypothesis is that the median distance to the nearest major road is shorter for mapped cocoa cells than for simulated cells drawn according to the environmental baseline. The secondary comparison asks the same question for the nearest urban centre. The proximity statistic is computed separately for roads and urban centres, so it measures each relationship on its own.

## Data and sources

| Input | Role in the analysis | Source |
| --- | --- | --- |
| Thresholded cocoa map | Approximately 10 m satellite-derived cocoa extent for Ghana and C&ocirc;te d'Ivoire, based on imagery from around 2019-2021. The workflow uses the 0.65-threshold raster. | [ETH Zurich Research Collection](https://doi.org/10.3929/ethz-b-000654400), [Kalischek et al.](https://doi.org/10.1038/s43016-023-00751-8) |
| Roads, settlements and boundary | The November 9, 2026 OpenStreetMap snapshot provides motorways and trunk, primary and secondary roads, plus mapped cities, towns and national capitals. The national administrative polygon defines the study area. | [Geofabrik Ghana](https://download.geofabrik.de/africa/ghana.html), [OpenStreetMap contributors and ODbL](https://www.openstreetmap.org/copyright) |
| Climate and elevation | Annual mean temperature (BIO1), annual precipitation (BIO12) and SRTM-derived elevation at 2.5 arc-minutes, roughly 4.5 km in Ghana. Climate normals cover 1970-2000. | [WorldClim 2.1](https://www.worldclim.org/data/worldclim21.html) |

The national run used 6,229 major-road line features and 577 urban-centre points after filtering the Geofabrik extract. Multiple line features can belong to the same road. The [data manifest](data/manifest.csv) records input filenames, versions, licences and checksums.

## Methodology

1. **Align the layers on a 1 km grid.** The workflow uses the Ghana National Grid (EPSG:2136) and converts projection units to metres for reporting. Eligible cells fall inside the working Ghana boundary and the cocoa raster coverage mask and have complete environmental predictors. Each cell with mapped cocoa contributes one observation at its centre. The threshold raster is reprojected using mean aggregation, and a cell is marked positive when the resulting value is greater than zero. WorldClim layers are aligned to the same grid by bilinear resampling. The workflow treats NaN within the raster coverage mask as no mapped cocoa, an assumption awaiting provider confirmation.
2. **Measure proximity.** Compute the straight-line distance from each eligible cell centre to the nearest major road and the nearest qualifying urban centre. The primary statistic is median road distance across cocoa-positive cells. Median urban-centre distance is secondary.
3. **Estimate an environmental baseline.** Fit a binomial generalized additive model with smooth terms for mean temperature, log-transformed annual precipitation and elevation. The primary model groups 20 km spatial blocks into five geographic folds using k-means on block coordinates and predicts each fold using a model fitted on the other four. The resulting scores weight how often cells are selected for simulations. Road and settlement distances are excluded from the model.
4. **Simulate and compare.** Generate 999 patterns per baseline, drawing 57,131 distinct cells from the same 240,294 eligible locations on every draw. This holds the overall amount of mapped occupancy fixed so the comparison isolates where positive cells fall. Compare the primary baseline with a 40 km block design, an environmental model fitted to all eligible cells, and uniform sampling over eligible cells. The corrected lower-tail Monte Carlo p-value is (1 + simulated medians at or below the observed median) / 1000. Maps compare observed occupancy with simulated inclusion frequency. An exploratory distance profile compares observed and simulated shares across distance bands with pointwise 95% simulation intervals.

## Results and interpretation

### How much closer is mapped cocoa?

![Observed median distances and distributions from 999 simulated patterns under each of the four baselines](outputs/figures/ghana_simulation_distance_distributions.png)

**Figure 1. Distribution of simulated median distances.** Each coloured curve summarizes 999 pattern-level medians. The dashed black line marks the observed median. The top panel uses distance to major roads, and the bottom uses distance to urban centres. The environmental curves sit closer to the observed line than the uniform curve in both panels, showing how much the baseline choice changes the comparison.

Of the 240,294 eligible 1 km cells, 57,131 contained mapped cocoa (23.8%). This is the prevalence of **positive grid cells** in the study window. A cell contributes equally whether it contains a small cocoa patch or extensive cultivation, so the percentage cannot be read as the share of Ghana's land planted with cocoa.

| Comparison baseline | Major road, simulated median (95% simulation interval), km | Urban centre, simulated median (95% simulation interval), km |
| --- | ---: | ---: |
| **Observed median** | **11.96** | **26.43** |
| Environmental model, 20 km spatial blocks (primary) | 13.25 (13.16-13.34) | 28.58 (28.46-28.70) |
| Environmental model, 40 km spatial blocks | 13.15 (13.06-13.23) | 28.39 (28.28-28.50) |
| Environmental model fitted to all cells | 13.00 (12.92-13.09) | 28.24 (28.13-28.36) |
| Uniform sampling over eligible cells | 19.19 (19.00-19.38) | 39.64 (39.37-39.93) |

Under the primary model, the observed median was **1.30 km shorter for major roads (9.8%)** and **2.15 km shorter for urban centres (7.5%)**, relative to the simulated medians. The direction persisted across the three environmental specifications. The gaps were much larger under uniform sampling, at 7.23 km for roads and 13.21 km for urban centres. Weighting simulated locations by climate and elevation therefore brings the predicted pattern substantially closer to the observed one. The remaining gaps indicate that these environmental predictors alone leave some of the observed proximity unexplained.

For each of the eight comparisons, zero of 999 simulated medians were at or below the observed median, yielding a corrected one-sided Monte Carlo p-value of **0.001**, the smallest value this simulation design can report. Each interval summarizes the middle 95% of simulated medians under its fitted baseline. Uncertainty from cocoa classification, the extent of eligible land, and model fitting is outside those intervals.

### Where along the distance gradient do the differences occur?

![Observed and simulated cocoa occupancy by distance to major roads and urban centres under all four baselines](outputs/figures/ghana_distance_binned_intensity.png)

**Figure 2. Cocoa occupancy across five distance bands.** The dashed black line shows the observed share of positive cells in each band. Each coloured line shows the median simulated share, with a shaded 95% interval for that band. Read across the top row for roads and the bottom row for urban centres. The primary 20 km comparison is in the first column.

Each percentage below is the fraction of **all eligible cells in that distance band** that were cocoa-positive or selected in a simulated pattern. The simulated column gives the median across the 999 primary environmental draws. This denominator matters because some distance bands cover much more land than others.

| Distance band | Road observed | Road simulated | Urban centre observed | Urban centre simulated |
| --- | ---: | ---: | ---: | ---: |
| 0-5 km | 34.4% | 32.3% | 33.7% | 37.1% |
| 5-10 km | 34.4% | 31.3% | 40.0% | 36.7% |
| 10-20 km | 30.8% | 29.0% | 40.3% | 35.9% |
| 20-50 km | 20.4% | 21.1% | 30.5% | 28.9% |
| Over 50 km | 2.1% | 7.6% | 7.9% | 11.5% |

The road profile shows more mapped cocoa than the baseline in every band within 20 km and less in both bands beyond 20 km. The biggest percentage-point deficit is beyond 50 km from a major road. The urban-centre pattern is different. Cocoa occupancy is below the baseline within 5 km and beyond 50 km, and above it from 5 to 50 km. The shorter urban-centre median reflects higher occupancy at intermediate distances and lower occupancy in the most remote band. Each band is an exploratory comparison with a pointwise simulation interval. The [Quarto manuscript](manuscript/manuscript.qmd) shows the profiles and the maps of observed occupancy and simulated inclusion.

### Where does the environmental model reproduce the mapped pattern?

![Maps of observed cocoa occupancy and mean simulated inclusion under the 20 km and 40 km environmental baselines](outputs/figures/ghana_environmental_intensity_comparison.png)

**Figure 3. Mapped cocoa and environmental simulation surfaces.** The top map assigns 1 to each cocoa-positive cell and 0 to each eligible cell without mapped cocoa. The middle and bottom maps show how often each cell was selected across 999 simulations. A value of 0.75 means the cell appeared in about three-quarters of draws. Both environmental baselines capture the broad concentration of cocoa in southern and southwestern Ghana and along a narrower stretch of the eastern border. The simulated surfaces are smoother and extend into areas where the observed map shows gaps.

![Map of observed occupancy minus mean inclusion under the primary 20 km environmental baseline](outputs/figures/ghana_environmental_intensity_residual.png)

**Figure 4. Cell-level differences from the primary environmental baseline.** Red marks cocoa-positive cells selected relatively infrequently in simulation. Blue marks cells without mapped cocoa that the simulation selected relatively often. For example, a cocoa-positive cell selected in 20% of draws has a residual of +0.80, while a cell with no mapped cocoa selected in 80% of draws has a residual of -0.80. These are descriptive differences. A prominent patch can guide follow-up on land cover, mapping quality or omitted predictors, but its colour alone gives no estimate of a road effect.

### How much confidence should we place in the environmental comparison?

| Spatial block size | Pooled held-out AUC | Mean AUC across evaluable folds | Fold AUC range | Folds with no cocoa-positive cells |
| --- | ---: | ---: | ---: | ---: |
| 20 km | 0.934 | 0.815 | 0.463-0.976 | 1 of 5 |
| 40 km | 0.934 | 0.819 | 0.458-0.987 | 1 of 5 |

AUC measures how well the model ranks cocoa-positive cells above other cells. The pooled score looks strong, but individual folds range from below chance-level ranking (0.5) to near-perfect ranking. Both designs have a fold where AUC is undefined because it contains no positives. Some fits also produced warnings about fitted probabilities near zero or one. The proximity result is consistent across the fitted baselines, with substantial uncertainty about how well those baselines represent environmental geography in individual parts of Ghana.

- **Mapped cocoa and eligible land.** Satellite classification errors can affect both the observed pattern and the environmental weights. Treating NaN as absence inside the coverage mask can turn unmapped cells into apparent absences. The raster bounding box overlaps 99.94% of the provisional OSM polygon by area. This measures overlap with the raster rectangle, so coverage within that rectangle remains subject to the NaN assumption. The workflow also lacks an independent land-cover mask for water and other unsuitable areas.
- **Spatial scale and timing.** WorldClim's roughly 4.5 km input resolution limits the environmental detail available on the 1 km analysis grid. The 1970-2000 climate normals describe long-run conditions, while the November 2026 road snapshot may include roads built after the cocoa imagery was acquired.
- **Access and causality.** Straight-line distances omit road connectivity, road condition and travel time. Mapped towns and cities approximate access to population centres. Estimating access to cocoa buyers would require verified purchasing sites. Historical settlement patterns, land use, and infrastructure placement could all contribute to the observed association. The analysis supports a provisional finding about where mapped cocoa occurs relative to roads and settlements. It leaves the causal effect of roads and residual clustering among cocoa locations open.


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

Intermediate spatial data go to data/processed/. Figures, tables and the rendered manuscript go to outputs/. The four PNGs above are outputs of the national reporting script.

| Generated table | Contents |
| --- | --- |
| ghana_simulation_summary.csv | Observed medians, simulated median percentiles, effect ratios and Monte Carlo p-values |
| ghana_distance_binned_intensity_summary.csv | Observed and simulated shares and the number of eligible cells in each distance band |
| ghana_environmental_cv_summary.csv | Pooled and fold-level performance summaries for both block sizes |
| ghana_simulation_qa.csv | Number of draws and unique simulation IDs, fixed sample size and checks for finite road and urban-centre statistics |
