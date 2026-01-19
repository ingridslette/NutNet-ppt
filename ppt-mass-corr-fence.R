library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(MuMIn)
library(performance)
library(MASS)
library(broom)
library(broom.mixed)
library(emmeans)
library(ggpubr)
library(purrr)
library(cowplot)
library(ggeffects)

### Loading, viewing, and filtering precipitation and mass data 

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/comb-by-plot-clim-soil-diversity_2024-05-31.csv",
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1_f <- filter(mass, year_trt > 0)
unique(mass1_f$year_trt)

mass1_f <- filter(mass1_f, trt %in% c("Fence", "NPK+Fence")) 
unique(mass1_f$trt)
unique(mass1_f$site_code)

mass1_f <- mass1_f %>%
  mutate(live_mass = case_when(
    !is.na(vascular_live_mass) | !is.na(nonvascular_live_mass) ~ 
      rowSums(across(c(vascular_live_mass, nonvascular_live_mass, standing_dead_mass)), na.rm = TRUE),
    is.na(vascular_live_mass) & is.na(nonvascular_live_mass) ~ 
      unsorted_live_mass
  )
  )

summary(mass1$live_mass)

site_year_counts_f <- mass1_f %>%
  group_by(site_code, trt) %>%
  filter(!is.na(live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years_f <- site_year_counts_f %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

mass2_f <- mass1_f %>%
  filter(site_code %in% sites_with_6_years_f$site_code)

unique(mass2_f$site_code)

## popped over to script "calculate-gs-ppt-pet.R" here, to get growing season ppt for the sites included in mass2
## exported that as csv and now loading it here

ppt_data <- read.csv("/Users/ingridslette/Desktop/NutNet/ppt_pet_annual_gs_only_2025-10-09.csv")

unique(ppt_data$site_code)

ppt_data <- ppt_data %>%
  arrange(site_code, year) %>% 
  group_by(site_code) %>%
  mutate(prev_ppt = lag(ppt)) %>%
  ungroup()

ppt_data <- filter(ppt_data, year >= 1983)
ppt_data <- filter(ppt_data, year < 2025) 
unique(ppt_data$year)

ppt_data <- ppt_data %>%
  group_by(site_code) %>%
  mutate(
    avg_ppt = mean(ppt, na.rm = TRUE),
    avg_pet = mean(pet, na.rm = TRUE),
    sd_ppt = sd(ppt, na.rm = TRUE),
    p10_ppt = quantile(ppt, 0.10, na.rm = TRUE),
    p90_ppt = quantile(ppt, 0.90, na.rm = TRUE),
    p05_ppt = quantile(ppt, 0.05, na.rm = TRUE),
    p95_ppt = quantile(ppt, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

unique(mass2$site_code)
unique(ppt_data$site_code)

mass_ppt_f <- inner_join(mass2_f, ppt_data, by = c("site_code", "year"))

unique(mass_ppt_f$site_code)

mass_ppt_f <- mass_ppt_f %>%
  mutate(log_mass = log10(live_mass),
         log_ppt = log10(ppt))

mass_ppt_f <- mass_ppt_f %>%
  group_by(site_code) %>%
  mutate(min_ppt = min(ppt, na.rm = TRUE),
         max_ppt = max(ppt, na.rm = TRUE)) %>%
  ungroup()


# Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
mass_ppt_f <- mass_ppt_f %>%
  group_by(site_code) %>%
  filter(min_ppt <= (avg_ppt - sd_ppt), max_ppt >= (avg_ppt + sd_ppt)) %>%
  ungroup()

unique(mass_ppt_f$site_code)
unique(mass_ppt_f$trt)


### Main model

mass_ppt_f <- mass_ppt_f[is.finite(mass_ppt_f$log_mass), ]

main_model_f <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year), 
                   data = mass_ppt_f, na.action = na.exclude)

summary(main_model_f)

r2_main_model_f <- performance::r2(main_model_f)


# Model assumptions check 
plot(main_model_f)
resid <- residuals(main_model_f)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid)
plot(fitted(main_model), resid, main = "Residuals vs Fitted")

summary(mass_ppt_f$live_mass)
hist(mass_ppt_f$live_mass)


### Back-transforming data for graphing - allows for non-linear curves on linear scale

# Back transform from log-log scale and graph
fit_model_and_predict_allsites <- function(data) {
  model <- lm(log_mass ~ log_ppt, data = data)
  new_data <- data.frame(log_ppt = seq(min(data$log_ppt, na.rm = TRUE),
                                       max(data$log_ppt, na.rm = TRUE),
                                       length.out = 100))
  preds <- predict(model, newdata = new_data, se.fit = TRUE)
  new_data$predicted_log_mass <- preds$fit
  new_data$se_log_mass <- preds$se.fit
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$mass_lower <- 10^(new_data$predicted_log_mass - 1.96 * new_data$se_log_mass)
  new_data$mass_upper <- 10^(new_data$predicted_log_mass + 1.96 * new_data$se_log_mass)
  new_data$trt <- unique(data$trt)
  return(new_data)
}


predictions_allsites_f <- mass_ppt_f %>%
  group_by(trt) %>%
  group_modify(~ fit_model_and_predict_allsites(.x)) %>%
  ungroup()

ggplot(data = mass_ppt_f, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_line(data = predictions_allsites_f, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m²)") +
  labs(color = "Treatment", shape = "Treatment") +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  theme_bw(base_size = 14)


fit_model_and_predict <- function(data) {
  model <- lm(log_mass ~ log_ppt, data = data)
  p_value <- summary(model)$coefficients["log_ppt", "Pr(>|t|)"]
  
  new_data <- data.frame(
    log_ppt = seq(min(data$log_ppt, na.rm = TRUE),
                  max(data$log_ppt, na.rm = TRUE),
                  length.out = 100)
  )
  new_data$predicted_log_mass <- predict(model, newdata = new_data)
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$p_value <- p_value
  return(new_data)
}

predictions_f <- mass_ppt_f %>%
  group_by(site_code, trt) %>%
  group_modify(~ fit_model_and_predict(.x)) %>%
  ungroup()

predictions_sig_f <- predictions_f %>%
  filter(!is.na(p_value) & p_value < 0.05)

ggplot(mass_ppt_f, aes(x = ppt, y = live_mass, color = trt, shape = trt, fill = trt)) +
  geom_point(alpha = 0.7) +
  geom_line(data = predictions_sig_f,
            aes(x = 10^log_ppt, y = predicted_mass), 
            linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)",
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("Fence" = "#0092E0", "NPK+Fence" = "#ff924c")) +
  scale_fill_manual(values = c("Fence" = "#0092E0", "NPK+Fence" = "#ff924c")) +
  scale_shape_manual(values = c("Fence" = 21, "NPK+Fence" = 24)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))


ggplot(mass_ppt_f, aes(x = ppt, y = live_mass, color = site_code)) +
  geom_line(data = predictions_f, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  geom_line(data = predictions_allsites_f, aes(x = 10^log_ppt, y = predicted_mass), 
            linewidth = 1, color = "black") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14)


### Comparing control vs. NPK R2 - calculate and compare difference at each site

site_codes_f <- unique(mass_ppt_f$site_code)

results_f <- data.frame(site_code = character(), 
                      control_r2 = numeric(), 
                      npk_r2 = numeric(),
                      r2_diff = numeric(),
                      control_slope = numeric(), 
                      npk_slope = numeric(),
                      slope_diff = numeric(),
                      stringsAsFactors = FALSE)

for (site in site_codes_f) {
  site_data_control <- subset(mass_ppt_f, site_code == site & trt == "Fence")
  site_data_npk <- subset(mass_ppt_f, site_code == site & trt == "NPK+Fence")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  results_f <- rbind(results_f, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_diff = npk_r2 - control_r2,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}

paired_t_test_r2_f <- t.test(results_f$control_r2, results_f$npk_r2, paired = TRUE)
paired_t_test_r2_f

