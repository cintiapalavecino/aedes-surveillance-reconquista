# ============================================================
# Trap positivity map for a municipal surveillance report
# Aedes aegypti surveillance - Reconquista, Santa Fe, Argentina
# Season 2024-2025
#
# Produced for routine reports to the Municipality of Reconquista.
# Plot labels are in Spanish because the audience was local health
# authorities.
# Author: Cintia C. Palavecino
# ============================================================

# 1. Libraries -----------------------------------------------------
library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)

# 2. Surveillance data ----------------------------------------------
Data <- read.delim("data/Reconquista_OVI_TRAP_weekly.csv",
                   header = TRUE, sep = ";", fileEncoding = "latin1") %>%
  filter(!is.na(Anio)) %>%
  rename(Year = Anio, EW = SE)

# Season definition: EW 39 of one year to EW 21 of the next
season_label <- function(year, ew) {
  case_when(
    (year == 2021 & ew >= 39) | (year == 2022 & ew <= 21) ~ "2021 - 2022",
    (year == 2022 & ew >= 39) | (year == 2023 & ew <= 21) ~ "2022 - 2023",
    (year == 2024 & ew >= 39) | (year == 2025 & ew <= 21) ~ "2024 - 2025"
  )
}

selected_season <- "2024 - 2025"

Data <- Data %>%
  mutate(Season = season_label(Year, EW)) %>%
  filter(Season == selected_season)

# 3. Positivity per trap ---------------------------------------------
# A trap-week counts only if the trap was active; positivity is the
# percentage of active weeks in which the trap was positive.
# Ovitraps (O) and adult gravid traps (T) share the same ID column.
positivity <- Data %>%
  mutate(
    # %in% handles numeric or text columns and treats NA as FALSE
    active   = Ovi_activa %in% 1 | Trampa_activa %in% 1,
    positive = Trampa_Ovi_positiva %in% 1
  ) %>%
  group_by(ID_Ovitrampa_trampa, Ovi_Trampa) %>%
  summarise(
    active_weeks   = sum(active, na.rm = TRUE),
    positive_weeks = sum(positive & active, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    positivity_pct = ifelse(active_weeks > 0, 100 * positive_weeks / active_weeks, NA),
    trap_type = recode(Ovi_Trampa, "O" = "Ovitrampa (huevos)", "T" = "Trampa de hembras adultas")
  )

# 4. Spatial layers ---------------------------------------------------
neighbourhoods <- st_read("data/maps/Barrios.shp", quiet = TRUE)
streets        <- st_read("data/maps/Calles.shp", quiet = TRUE) %>%
  st_transform(st_crs(neighbourhoods))
trap_points    <- st_read("data/maps/puntos_ovi.shp", quiet = TRUE) %>%
  st_transform(st_crs(neighbourhoods)) %>%
  rename(ID_Ovitrampa_trampa = id)

# Join by trap ID (attribute join, not a spatial join)
positivity_sf <- trap_points %>%
  mutate(ID_Ovitrampa_trampa = as.character(ID_Ovitrampa_trampa)) %>%
  inner_join(positivity %>% mutate(ID_Ovitrampa_trampa = as.character(ID_Ovitrampa_trampa)),
             by = "ID_Ovitrampa_trampa")

# Report traps without location so they are not silently dropped
missing_ids <- setdiff(positivity$ID_Ovitrampa_trampa, positivity_sf$ID_Ovitrampa_trampa)
if (length(missing_ids) > 0) {
  message("Traps without coordinates: ", paste(missing_ids, collapse = ", "))
}

# 5. Map ----------------------------------------------------------------
bbox <- st_bbox(neighbourhoods)

positivity_map <- ggplot() +
  geom_sf(data = neighbourhoods, fill = "white", color = "grey50") +
  geom_sf(data = streets, color = "grey80", linewidth = 0.3) +
  geom_sf(data = positivity_sf,
          aes(color = positivity_pct, shape = trap_type),
          size = 7, alpha = 0.85) +
  geom_sf_text(data = positivity_sf,
               aes(label = round(positivity_pct)),
               size = 3.5, nudge_y = 60) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
           ylim = c(bbox["ymin"], bbox["ymax"]),
           expand = FALSE) +
  scale_color_gradientn(
    colors = c("#2ECC71", "#F39C12", "#AE123A"),
    limits = c(0, 100),
    name = "% de semanas con\nactividad del vector"
  ) +
  scale_shape_manual(values = c(16, 17), name = "Tipo de trampa") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "tr", style = north_arrow_minimal()) +
  labs(
    title = "Porcentaje de positividad por trampa",
    subtitle = paste("Temporada", selected_season, "| Reconquista, Santa Fe")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

positivity_map

# 6. Save -----------------------------------------------------------------
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/positivity_map.png", positivity_map,
       width = 12, height = 10, dpi = 300, bg = "white")
