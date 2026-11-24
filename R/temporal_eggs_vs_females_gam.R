# ============================================================
# Temporal analysis: ovitraps (eggs) vs. gravid traps (adult females)
# Aedes aegypti surveillance - Reconquista, Santa Fe, Argentina
# Season 2024-2025 (epidemiological week 36/2024 to 20/2025)
# Count models with sampling-effort offsets, fitted by REML
# Author: Cintia C. Palavecino
# ============================================================

# 1. Libraries -----------------------------------------------------
library(dplyr)
library(ggplot2)
library(mgcv)
library(tidyr)

# 2. Data loading and preparation -----------------------------------
data_path <- "data/Reconquista_OVI_TRAP_weekly.csv"

Data <- read.delim(data_path, header = TRUE, sep = ";", fileEncoding = "latin1") %>%
  mutate(Date = as.Date(Fecha_sem_en_terreno, format = "%d/%m/%Y"))

Data_season <- Data %>%
  filter((Anio == 2024 & SE >= 36) | (Anio == 2025 & SE <= 20)) %>%
  mutate(EW = SE, Year = Anio, EW_season = Serie_temporada)

# 3. Weekly aggregation ---------------------------------------------
# Numerator and denominator use exactly the same units (active traps
# WITH a recorded count). This excludes inactive ovitraps that contain
# eggs and active traps without a count.
weekly <- Data_season %>%
  group_by(Year, EW, EW_season) %>%
  summarise(
    Total_eggs    = sum(N_total_ovis[Ovi_activa %in% 1], na.rm = TRUE),
    Active_ovi    = sum(Ovi_activa %in% 1 & !is.na(N_total_ovis)),
    Total_females = sum(Hembras_Ae_aegypti_trampas[Trampa_activa %in% 1], na.rm = TRUE),
    Active_traps  = sum(Trampa_activa %in% 1 & !is.na(Hembras_Ae_aegypti_trampas)),
    .groups = "drop"
  ) %>%
  arrange(EW_season) %>%
  mutate(
    Mean_eggs    = Total_eggs / Active_ovi,
    Mean_females = Total_females / Active_traps,
    EW_label     = paste0(Year, "-EW", sprintf("%02d", EW)),
    Eggs_z       = as.numeric(scale(Mean_eggs)),
    Females_z    = as.numeric(scale(Mean_females))
  )

# 3b. Check the series before building lags --------------------------
# lag() shifts rows, not weeks: if a week were missing, the lag would be
# misaligned. EW_season must therefore be consecutive.
stopifnot(all(diff(weekly$EW_season) == 1))
stopifnot(all(weekly$Active_ovi > 0), all(weekly$Active_traps > 0))

weekly <- weekly %>%
  mutate(
    eggs_lag1 = lag(Mean_eggs, 1),
    eggs_lag2 = lag(Mean_eggs, 2)
  )

# Common style for all figures
fig_theme  <- theme_minimal(base_size = 16) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
stage_cols    <- c(Eggs = "#D55E00", Females = "#56B4E9")
labels_lc  <- c(Eggs = "Eggs",    Females = "Adult females")

check_gam <- function(model) {
  par(mfrow = c(2, 2)); gam.check(model); par(mfrow = c(1, 1))
}

# 4. Standardized observed time series -------------------------------
weekly_long <- weekly %>%
  select(EW, EW_season, Eggs_z, Females_z) %>%
  pivot_longer(c(Eggs_z, Females_z), names_to = "Method", values_to = "z") %>%
  mutate(Method  = recode(Method, Eggs_z = "Eggs", Females_z = "Females"),
         EW_plot = factor(EW, levels = unique(EW)))

Fig_temporal <- ggplot(weekly_long, aes(EW_plot, z, color = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = stage_cols, labels = labels_lc) +
  fig_theme +
  labs(x = "Epidemiological week", y = "Standardized abundance (z-score)", color = "Life stage")
Fig_temporal

# 5. Univariate GAMs: total counts + sampling-effort offset ------------
# log(E[total]) = log(n active traps) + b0 + f(t)  =>  models the mean per
# trap without dividing the counts. Negative binomial for both counts.
gam_eggs <- gam(Total_eggs ~ s(EW_season, k = 10) + offset(log(Active_ovi)),
                family = nb(), method = "REML", data = weekly)

gam_females <- gam(Total_females ~ s(EW_season, k = 10) + offset(log(Active_traps)),
                   family = nb(), method = "REML", data = weekly)

summary(gam_eggs);    check_gam(gam_eggs)
summary(gam_females); check_gam(gam_females)

# Predictions on the "per active trap" scale (offset = log(1))
nd <- weekly %>% mutate(Active_ovi = 1, Active_traps = 1)
weekly <- weekly %>%
  mutate(
    eggs_pred    = as.numeric(predict(gam_eggs,    newdata = nd, type = "response")),
    females_pred = as.numeric(predict(gam_females, newdata = nd, type = "response"))
  )

peak_eggs    <- weekly$EW_season[which.max(weekly$eggs_pred)]
peak_females <- weekly$EW_season[which.max(weekly$females_pred)]
cat("Egg peak (EW):", weekly$EW_label[weekly$EW_season == peak_eggs],
    "| Female peak (EW):", weekly$EW_label[weekly$EW_season == peak_females], "\n")

Fig_gam <- ggplot(weekly, aes(EW_season)) +
  geom_point(aes(y = Eggs_z,    color = "Eggs"),    alpha = 0.6, size = 2) +
  geom_point(aes(y = Females_z, color = "Females"), alpha = 0.6, size = 2) +
  geom_line(aes(y = as.numeric(scale(eggs_pred)),    color = "Eggs"),    linewidth = 1.4) +
  geom_line(aes(y = as.numeric(scale(females_pred)), color = "Females"), linewidth = 1.4) +
  geom_vline(xintercept = peak_eggs,    linetype = "dashed", color = stage_cols["Eggs"]) +
  geom_vline(xintercept = peak_females, linetype = "dashed", color = stage_cols["Females"]) +
  scale_color_manual(values = stage_cols, labels = labels_lc) +
  scale_x_continuous(breaks = weekly$EW_season, labels = weekly$EW) +
  fig_theme +
  labs(x = "Epidemiological week", y = "Standardized smooth", color = "Life stage")
Fig_gam

# 6. Lagged GAMs (eggs at t-1 or t-2 -> females at t) -----------------
# Both models are fitted on the SAME weeks (from the 3rd onwards), so that
# AIC and deviance explained are comparable between lag 1 and lag 2.
lag_data <- weekly %>% filter(!is.na(eggs_lag2))

gam_lag1 <- gam(Total_females ~ s(eggs_lag1, k = 5) + s(EW_season, k = 8) +
                  offset(log(Active_traps)),
                family = nb(), method = "REML", data = lag_data)
gam_lag2 <- gam(Total_females ~ s(eggs_lag2, k = 5) + s(EW_season, k = 8) +
                  offset(log(Active_traps)),
                family = nb(), method = "REML", data = lag_data)

summary(gam_lag1); check_gam(gam_lag1)
summary(gam_lag2); check_gam(gam_lag2)
AIC(gam_lag1, gam_lag2)

# Summary table of smooth terms
gam_table <- function(m, name) {
  s <- summary(m)
  data.frame(Model = name, Term = rownames(s$s.table),
             edf = round(s$s.table[, "edf"], 2),
             Chi2 = round(s$s.table[, "Chi.sq"], 2),
             p = signif(s$s.table[, "p-value"], 3),
             Dev_expl = round(100 * s$dev.expl, 1),
             n = s$n, row.names = NULL)
}
summary_table <- bind_rows(
  gam_table(gam_eggs,    "Eggs ~ s(t)"),
  gam_table(gam_females, "Females ~ s(t)"),
  gam_table(gam_lag1,    "Females ~ s(eggs t-1) + s(t)"),
  gam_table(gam_lag2,    "Females ~ s(eggs t-2) + s(t)")
)
print(summary_table)

# 7. Partial effects (lagged eggs -> females) ---------------------------
partial_effect <- function(model, var, data, n = 100) {
  rng  <- range(data[[var]], na.rm = TRUE)
  grid <- data.frame(x = seq(rng[1], rng[2], length.out = n))
  names(grid) <- var
  grid$EW_season    <- median(data$EW_season)
  grid$Active_traps <- 1
  pr  <- predict(model, newdata = grid, type = "iterms", se.fit = TRUE)
  col <- grep(var, colnames(pr$fit), fixed = TRUE)
  data.frame(x = grid[[var]], fit = pr$fit[, col], se = pr$se.fit[, col])
}

df_lags <- bind_rows(
  partial_effect(gam_lag1, "eggs_lag1", lag_data) %>% mutate(Model = "1-week lag"),
  partial_effect(gam_lag2, "eggs_lag2", lag_data) %>% mutate(Model = "2-week lag")
)

# Rug: each panel shows the observed values of ITS OWN lag
rug_lags <- bind_rows(
  lag_data %>% transmute(x = eggs_lag1, Model = "1-week lag"),
  lag_data %>% transmute(x = eggs_lag2, Model = "2-week lag")
)

# p-values inside each panel
p_lab <- data.frame(
  Model = c("1-week lag", "2-week lag"),
  lab = paste0("p = ", formatC(c(summary(gam_lag1)$s.table[1, 4],
                                  summary(gam_lag2)$s.table[1, 4]),
                                format = "f", digits = 3))
)

lag_cols <- c("1-week lag" = "#D55E00", "2-week lag" = "#009E73")

Fig_lags <- ggplot(df_lags, aes(x, fit, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_rug(data = rug_lags, aes(x = x), inherit.aes = FALSE, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_text(data = p_lab, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = -0.2, vjust = 1.5, size = 5) +
  scale_color_manual(values = lag_cols) +
  scale_fill_manual(values = lag_cols) +
  facet_wrap(~ Model, scales = "free_x") +
  fig_theme +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 14)) +
  labs(x = "Eggs per active ovitrap",
       y = "Partial effect on females\n(linear predictor scale)")
Fig_lags

# 8. Save figures -------------------------------------------------------
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

save_fig <- function(fig, name, width, height) {
  # PNG at 600 dpi
  ggsave(file.path("output/figures", paste0(name, ".png")), fig,
         width = width, height = height, units = "cm", dpi = 600, bg = "white")
  # Vector PDF for print
  ggsave(file.path("output/figures", paste0(name, ".pdf")), fig,
         width = width, height = height, units = "cm", device = cairo_pdf)
}

save_fig(Fig_gam,  "Fig1_gam",  width = 30, height = 14)
save_fig(Fig_lags, "Fig2_lags", width = 30, height = 14)
