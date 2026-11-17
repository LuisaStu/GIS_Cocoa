# Analysis decisions

## Ghana-wide analysis specification

- Use the Ghana national OSM polygon from the already-downloaded Geofabrik extract. Treat the OSM national boundary vintage as provisional and document its source.
- Use EPSG:2136 (Accra / Ghana National Grid) for the national distance and raster grid. This CRS uses Gold Coast feet, so define a grid cell as 1,000 m converted to Gold Coast feet and convert all reported distances and areas to metric units.
- Use the 0.65-threshold ETH cocoa raster at a 1 km grid-cell-centre occupancy unit. Per project-owner instruction, interpret source NaN as no cocoa only within the source raster footprint. Exclude cells beyond that footprint rather than coding missing coverage as absence; report the spatial coverage limitation.
- Retain only national grid-cell centres inside Ghana with cocoa-source coverage and complete BIO1, BIO12, and elevation values. Use BIO1, log-transformed BIO12, and elevation as GAM predictors, with five-fold spatial-block cross-fitting.
- Use 20 km spatial blocks for the primary cross-fitted environmental null and 40 km blocks as a pre-specified block-size sensitivity. Keep the in-sample environmental GAM and region-uniform fixed-count null as comparisons.
- Define the secondary urban-centre proxy as all Geofabrik `place=city`, `place=town`, and `place=national_capital` points inside the national boundary. It is not a verified cocoa-market inventory. Continue to exclude sparse marketplace POIs as a primary outcome.
- Report median nearest-major-road distance as the primary statistic and median distance to the urban-centre proxy as a secondary statistic. Use 999 fixed-count, no-replacement simulations for each null and the corrected lower-tail Monte Carlo p-value.

## Ghana-wide diagnostics and reporting

- Use the pooled out-of-fold AUC together with the mean valid fold AUC and per-fold range. A fold with no positive observations has undefined AUC and must be reported, not silently imputed. Pooled scores combine predictions from separate fold fits and are not interpreted as a standalone measure of national predictive performance.
- Report the 999-draw median-distance intervals and corrected lower-tail Monte Carlo p-values for both road and urban-centre distances. A p-value of 0.001 is the minimum resolution of this design.
- As a secondary, exploratory functional diagnostic, compare observed positive-cell share and simulated positive-cell shares in fixed distance bands: 0–5, 5–10, 10–20, 20–50, and over 50 km from roads or urban-centre proxies. Report pointwise 2.5th–97.5th percentile simulation intervals; do not describe them as simultaneous/global envelopes or as a pre-registered confirmatory test.
- Compute the national area and raster-bounding-box coverage from the current Geofabrik OSM polygon, but treat the polygon-derived area as a diagnostic only. The boundary source and vintage require independent validation before making definitive national-area or coverage claims.
