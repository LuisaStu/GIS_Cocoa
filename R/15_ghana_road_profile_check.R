# Exploratory whole-profile check using the saved primary baseline simulations.
tables_dir <- file.path("outputs", "tables")
draws <- read.csv(file.path(tables_dir, "ghana_distance_binned_intensity_draws.csv"))
road <- subset(draws, null_model == "spatial_block_20km_crossfit" &
                 statistic == "nearest_major_road")
bands <- c("0-5 km", "5-10 km", "10-20 km", "20-50 km", "50+ km")
profiles <- vapply(bands, function(band) {
  rows <- road[road$distance_bin == band, ]
  rows <- rows[order(rows$simulation_id), ]
  stopifnot(identical(rows$simulation_id, seq_len(999L)),
            length(unique(rows$observed_intensity)) == 1L)
  c(rows$observed_intensity[1], rows$simulated_intensity)
}, numeric(1000L))

# Pool observed and simulated profiles so every pattern gets the same scoring rule.
stopifnot(all(is.finite(profiles)), all(apply(profiles, 2, sd) > 0))
scores <- apply(abs(scale(profiles)), 1, max)
as_extreme <- sum(scores[-1] >= scores[1])
result <- data.frame(
  simulation_count = length(scores) - 1L,
  observed_max_standardized_deviation = scores[1],
  greatest_simulated_deviation = max(scores[-1]),
  simulated_at_least_as_extreme = as_extreme,
  monte_carlo_p = (1 + as_extreme) / length(scores)
)
write.csv(result, file.path(tables_dir, "ghana_road_profile_check.csv"), row.names = FALSE)
print(result)
