# *Aedes aegypti* surveillance in Reconquista, Argentina

R code from an integrated dengue surveillance study in **Reconquista (Santa Fe, Argentina)**, carried out with the local government during the 2024–2025 season. Two surveillance methods are used: **ovitraps (eggs)** and **gravid *Aedes* traps (adult females)**.

The repository has two parts: a statistical analysis presented at a scientific meeting, and a reporting tool used to communicate results to the municipality.

## 1. Do eggs anticipate adult female abundance?

> Palavecino, C. C., Rodríguez, S. T., Almada, M. S., & López, M. S. (2026). *¿Los huevos anticipan la abundancia de hembras adultas en* Aedes aegypti*? Un análisis temporal mediante modelos aditivos generalizados en Reconquista (Santa Fe).* Poster presented at the II Congreso Argentino de Estadística, Tandil, Argentina, 14–16 October 2026.

`R/temporal_eggs_vs_females_gam.R`

1. **Weekly aggregation** by epidemiological week, counting only active traps with a recorded count, so that numerator and denominator refer to the same sampling units.
2. **Series check:** the script verifies that weeks are consecutive before building lagged variables, because `lag()` shifts rows rather than weeks.
3. **Seasonal GAMs:** negative binomial GAMs on total counts with `offset(log(active traps))`, fitted by REML (`mgcv`), with predictions per active trap and smoothed seasonal peaks.
4. **Lagged GAMs:** female counts modelled as a smooth function of eggs per ovitrap at *t*−1 and *t*−2, plus a seasonal smooth. Both models are fitted on the same weeks so AIC and deviance explained are comparable.
5. **Outputs:** summary table of smooth terms (edf, χ², p, deviance explained), partial-effect plots with confidence bands, and figures exported as PNG (600 dpi) and vector PDF.

## 2. Positivity map for municipal reports

`R/positivity_map_municipal_report.R`

This script produced maps for routine surveillance reports to the **Municipality of Reconquista**, to help local health authorities identify areas with persistent vector activity. Plot labels are in Spanish for that audience.

- Positivity per trap = % of **active** weeks with a positive result (ovitraps and adult traps)
- Attribute join between surveillance records and trap locations, with a check for traps missing coordinates
- Map with neighbourhoods, streets, scale bar and north arrow (`sf`, `ggplot2`, `ggspatial`), distinguishing trap types by shape

## Data

The surveillance data are not public. To run the scripts, place in `data/`:

- `Reconquista_OVI_TRAP_weekly.csv`: weekly records per trap (`;` separator, `dd/mm/yyyy` dates)
- `maps/Barrios.shp`, `maps/Calles.shp`, `maps/puntos_ovi.shp`: neighbourhoods, streets and trap locations

Run the scripts from the repository root. Figures are written to `output/figures/`.

## Main packages

`dplyr`, `tidyr`, `ggplot2`, `mgcv`, `sf`, `ggspatial`

## Author

Cintia C. Palavecino, CONICET / CEVARCAM, Universidad Nacional del Litoral, Argentina
