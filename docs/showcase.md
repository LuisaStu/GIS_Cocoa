# Ghana cocoa accessibility showcase

The course lists an in-class poster session for Wednesday, December 2, 2026. The Week 12 description uses a conference-table format where visitors circulate and ask about projects. The current Week 14 page gives no poster dimensions, slide count, or talk duration. [Week 14 showcase](https://jjacobs.me/ppol6805/w14/index.html) · [Week 12 format](https://jjacobs.me/ppol6805/w12/slides.html)

## Poster layout

| Panel | What to show |
| --- | --- |
| Question and study area | **Is mapped cocoa in Ghana closer to major roads than climate and elevation predict?** Use the [national overlay](../outputs/figures/ghana_cocoa_access_map.png) and point out the southern and eastern cocoa belts. |
| Spatial comparison | Show how the analysis selects 57,131 of 240,294 eligible 1 km cells in each of 999 simulated maps. Explain that rainfall, temperature and elevation determine the reference weights. |
| Main result | Show the [median distance distributions](../outputs/figures/ghana_simulation_distance_distributions.png). Highlight **11.96 km observed versus 13.25 km simulated**, a difference of **1.30 km** using unrounded values. |
| Where the model misses | Show the [quadrat comparison](../outputs/figures/ghana_quadrat_intensity.png). Explain that the model reproduces the broad cocoa belt while missing local variation. |
| Interpretation | State that the finding is a spatial association under the fitted baseline. Historical roads and verified cocoa purchasing sites would help test actual market access. |

## Short walkthrough

“We wanted to know whether mapped cocoa in Ghana is closer to major roads than we would expect from climate and elevation. We held the number of cocoa-positive grid cells fixed and simulated 999 maps under our primary environmental model. The observed median road distance was about 12 kilometres, versus about 13.3 kilometres in the simulated maps. The result was consistent across all four baselines, but the environmental model performed poorly in the held-out region with most of the mapped cocoa. The next step is to validate the baseline there and use travel times to verified cocoa buyers.”

## Questions to prepare for

- **Why a first-order hypothesis?** We ask how occupied-cell intensity varies with location and road distance. A second-order analysis would ask how cocoa locations relate to one another after accounting for that intensity.
- **Why hold the count fixed?** Equal counts let us compare where cocoa occurs across patterns with the same total occupancy.
- **What does p = 0.001 mean?** None of the 999 primary draws had a median road distance as short as the observed one. The result is conditional on the fitted reference model.
- **Why call this a market-access proxy?** The data measure straight-line proximity to mapped roads and settlements. Actual access depends on road quality, network travel time and the locations of cocoa buyers.

Keep a local copy of the [rendered manuscript](../outputs/index.html) available during the session. Its figures are embedded and can be opened in the browser.

Use the [3D accessibility maps](../outputs/figures/ghana_accessibility_terrain.png) if visitors ask what the distance measures capture. The panels use the same kilometre scale for major roads and urban centres.
