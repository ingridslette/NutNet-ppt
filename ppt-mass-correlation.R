library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(MuMIn)
library(performance)
library(MASS)
library(broom)
library(emmeans)
library(ggpubr)


### Loading, viewing, and filtering precipitation and mass data 

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/comb-by-plot-clim-soil_2025-04-14.csv",
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, year_trt > 0)
unique(mass1$year_trt)
unique(mass1$site_code)

mass1 <- filter(mass1, trt %in% c("Control", "NPK")) 
unique(mass1$trt)
unique(mass1$site_code)

mass1 <- mass1 %>%
  mutate(
    live_mass = if_else(
      if_all(c(vascular_live_mass, nonvascular_live_mass), is.na),
      NA_real_,
      rowSums(across(c(vascular_live_mass, nonvascular_live_mass)), na.rm = TRUE)
    )
  )

site_year_counts <- mass1 %>%
  group_by(site_code, trt) %>%
  filter(!is.na(live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years <- site_year_counts %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

mass2 <- mass1 %>%
  filter(site_code %in% sites_with_6_years$site_code)

unique(mass2$site_code)

## popped over to script "daily-to-gs-ppt.R" here, to get growing season ppt for the sites in that list

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/mswep_ppt_annual_gs_only_2025-04-15.csv")

mswep <- mswep %>%
  arrange(site_code, year) %>% 
  group_by(site_code) %>%
  mutate(prev_ppt = lag(mswep_ppt)) %>%
  ungroup()

mswep <- filter(mswep, year >= 1983)

mswep <- mswep %>%
  group_by(site_code) %>%
  mutate(avg_ppt_site = mean(mswep_ppt, na.rm = TRUE),
         sd_ppt_site = sd(mswep_ppt, na.rm = TRUE),
         mswep_ppt_per = (mswep_ppt / avg_ppt_site) * 100, 
         mswep_ppt_sd = (mswep_ppt - avg_ppt_site) / sd_ppt_site) %>%
  ungroup()


mass_ppt <- left_join(mass2, mswep, by=c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(vascular_live_mass),
         log_mswep_ppt = log10(mswep_ppt))

unique(mass_ppt$trt)


# Filter to keep only sites with a certain ppt range
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  group_by(site_code) %>%
  filter(max(mswep_ppt_sd) > 1 & min(mswep_ppt_sd) < -1) %>%
  ungroup()

unique(mass_ppt_c_npk$site_code)

### Initial graphs
ggplot(data = mass_ppt_c_npk, aes(x = mswep_ppt, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x=year_trt, y=vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Treatment Year") + ylab("Biomass (g/m2)") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt, y= vascular_live_mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation (mm)",
       y = "Biomass (g/m2)",
       color = "Site Code") +
  theme(legend.position = "right")


### Graphing back-transformed data

fit_model_and_predict <- function(data) {
  model <- lm(log_mass ~ log_mswep_ppt, data = data)
  new_data <- data.frame(log_mswep_ppt = seq(min(data$log_mswep_ppt, na.rm = TRUE),
                                             max(data$log_mswep_ppt, na.rm = TRUE),
                                             length.out = 100))
  new_data$predicted_log_mass <- predict(model, newdata = new_data)
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$site_code <- unique(data$site_code)
  new_data$trt <- unique(data$trt)
  return(new_data)
}

predictions <- mass_ppt_c_npk %>%
  group_by(site_code, trt) %>%
  group_modify(~ fit_model_and_predict(.x)) %>%
  ungroup()

ggplot(mass_ppt_c_npk, aes(x = mswep_ppt, y = vascular_live_mass, color = site_code)) +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x = mswep_ppt, y = vascular_live_mass, color = trt)) +
  geom_point() +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

fit_model_and_predict_allsites <- function(data) {
  model <- lm(log_mass ~ log_mswep_ppt, data = data)
  new_data <- data.frame(log_mswep_ppt = seq(min(data$log_mswep_ppt, na.rm = TRUE),
                                             max(data$log_mswep_ppt, na.rm = TRUE),
                                             length.out = 100))
  new_data$predicted_log_mass <- predict(model, newdata = new_data)
  new_data$predicted_mass <- 10^new_data$predicted_log_mass
  new_data$trt <- unique(data$trt)
  return(new_data)
}

predictions_allsites <- mass_ppt_c_npk %>%
  group_by(trt) %>%
  group_modify(~ fit_model_and_predict_allsites(.x)) %>%
  ungroup()

ggplot(mass_ppt_c_npk, aes(x = mswep_ppt, y = vascular_live_mass, color = site_code)) +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            linewidth = 1, color = "black") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14)

ggplot(data = mass_ppt_c_npk,aes(x= mswep_ppt, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + 
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  labs(color = "Treatment", shape = "Treatment") +
  scale_color_manual(values = c("#4267ac", "#ff924c")) +
  theme_bw(14)


### Initial model

c_npk_x_model <- lmer(log_mass ~ log_mswep_ppt * trt + (1 | site_code / year_trt), data = mass_ppt_c_npk)
summary(c_npk_x_model)

# Model assumptions check 
plot(c_npk_x_model)

resid <- residuals(c_npk_x_model)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid)

plot(fitted(c_npk_x_model), resid, main = "Residuals vs Fitted")


### Comparing control vs. NPK R2 - Approach 1: calculate and compare difference at each site

# Get unique site codes
site_codes <- unique(mass_ppt_c_npk$site_code)

results <- data.frame(site_code = character(), 
                      control_r2 = numeric(), 
                      npk_r2 = numeric(), 
                      r2_difference = numeric(),
                      control_slope = numeric(), 
                      npk_slope = numeric(), 
                      slope_difference = numeric(),
                      stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  r2_difference <- control_r2 - npk_r2
  control_slope <- coef(control_model)["log_mswep_ppt"]
  npk_slope <- coef(npk_model)["log_mswep_ppt"]
  slope_difference <- control_slope - npk_slope
  results <- rbind(results, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_difference = r2_difference,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_difference = slope_difference
  ))
}

# t-test on r2 differences
t_test_r2_diff <- t.test(results$r2_difference)
print(t_test_r2_diff)

# paired t-test on r2 differences
paired_t_test_r2 <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
print(paired_t_test_r2)

# sites with higher r2 in NPK plots
sites_higher_r2_npk <- filter(results, r2_difference < 0)

# t-test on slope differences 
t_test_slope_diff <- t.test(results$slope_difference)
print(t_test_slope_diff)

# paired t-test on slope differences 
paired_t_test_slope <- t.test(results$control_slop, results$npk_slope, paired = TRUE)
print(paired_t_test_slope)


### Comparing control vs. NPK R2 - Approach 2: fit separate models for control and NPK data, calculate and compare z scores

model_control <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                      data = subset(mass_ppt_c_npk, trt == "Control"))
model_npk <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                  data = subset(mass_ppt_c_npk, trt == "NPK"))

summary(model_control)
summary(model_npk)

AIC(model_control, model_npk)

r2_control <- performance::r2(model_control)
r2_npk <- performance::r2(model_npk)

conditional_r2_control <- r2_control$R2_conditional
conditional_r2_npk <- r2_npk$R2_conditional

marginal_r2_control <- r2_control$R2_marginal
marginal_r2_npk <- r2_npk$R2_marginal

fixed_effects_control <- fixef(model_control)
fixed_effects_npk <- fixef(model_npk)

slope_control <- fixed_effects_control["log_mswep_ppt"]
slope_npk <- fixed_effects_npk["log_mswep_ppt"]

# Compare R2 values using Fisher's Z transformation
z_control <- 0.5 * log((1 + sqrt(conditional_r2_control)) / (1 - sqrt(conditional_r2_control)))
z_npk <- 0.5 * log((1 + sqrt(conditional_r2_npk)) / (1 - sqrt(conditional_r2_npk)))

n_control <- length(unique(subset(mass_ppt_c_npk, trt == "Control")$site_code))
n_npk <- length(unique(subset(mass_ppt_c_npk, trt == "NPK")$site_code))
se_diff <- sqrt((1 / (n_control - 3)) + (1 / (n_npk - 3)))

# Calculate the Z-score for the difference
z_diff <- (z_control - z_npk) / se_diff

p_value <- 2 * (1 - pnorm(abs(z_diff)))

cat("Z-score for the difference:", z_diff, "\n")
cat("P-value for the difference in conditional R-squared values:", p_value, "\n")


### Comparing control vs. NPK R2 - Approach 3: Bootstrapping to test for difference in R2 between Control and NKP models

r2_diff <- function(data, indices) {
  data_resampled <- data[indices, ]
  model_control <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code), 
                        data = data_resampled[data_resampled$trt == "Control", ])
  model_npk <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code), 
                    data = data_resampled[data_resampled$trt == "NPK", ])
  r2_control <- r.squaredGLMM(model_control)[2]
  r2_npk <- r.squaredGLMM(model_npk)[2]          
  return(r2_control - r2_npk)
}

set.seed(123)

boot_r2 <- boot(data = mass_ppt_c_npk, statistic = r2_diff, R = 1000)

print(boot_r2)

boot_ci <- boot.ci(boot_r2, type = "perc")
print(boot_ci)


## Updating graphs to order by R2 and change color palette 

predictions <- predictions %>%
  left_join(results, by = "site_code")

predictions <- predictions %>%
  mutate(site_code = factor(site_code, levels = unique(site_code[order(r2_difference)])))

pal2 <- c("#800000","#c00000","#ff0000","#ff4040","#ff8080","#a83a01","#e04d01","#f06201","#ff7700","#e0a500",
          "#ffbc00","#ffcd40","#ffde80","#305020","#406a2a","#609f3f","#80d353","#bdda0f","#0b5043","#117864",
          "#1abc9c","#03045e","#125d93","#057dcd","#43b0f1","#96cff1","#503658","#80558c","#af7ab3","#c497b0",
          "#808080")

pal3 <- c("#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080")

ggplot(predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass, colour = site_code)) +
  geom_line() +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14) +
  theme(legend.position = "none")

ggplot(predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass, colour = site_code)) +
  geom_line() +
  scale_color_manual(values = pal2) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14) +
  theme(legend.position = "none")

ggplot(predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass, colour = site_code)) +
  geom_line() +
  scale_color_manual(values = pal3) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14) +
  theme(legend.position = "none")


### Calculating log response ratios

lrr_df <- mass_ppt_c_npk %>%
  group_by(site_code) %>%
  summarize(
    lrr_mass = log(mean(vascular_live_mass[trt == "NPK"], na.rm = TRUE) /
                     mean(vascular_live_mass[trt == "Control"], na.rm = TRUE)),
    lrr_prop_par = log(mean(proportion_par[trt == "NPK"], na.rm = TRUE) /
                         mean(proportion_par[trt == "Control"], na.rm = TRUE))
  )

ggplot(lrr_df, aes(x = lrr_prop_par, y = lrr_mass)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkgrey") +
  labs(x = "LRR proportion par",
       y = "LRR mass") +
  theme_bw()

lrr_par_mass_model <- lm(lrr_mass ~ lrr_prop_par, data = lrr_df)
summary(lrr_par_mass_model)


## incorporating C3/C4 and annual/perennial information from cover data

cover <- read.csv("/Users/ingridslette/Desktop/NutNet/full-cover_2025-01-31.csv",
                 na.strings = c("NULL","NA"))

cover <- cover %>%
  filter(site_code %in% site_codes)

dat1 <- subset(cover, is.na(ps_path) == TRUE)%>%
  separate(Taxon, into = c("Genus", "Species"), remove = FALSE, sep = " ")


dat1$ps_path <- ifelse(dat1$Genus == "BINERTIA" | dat1$Genus == "TIDESTROMIA" | dat1$Genus == "PECTIS" | dat1$Genus == "EUPLOCA" | dat1$Genus == "BULBOSTYLIS" | dat1$Genus == "CYPERUS" | dat1$Genus == "FIMBRISTYLIS" | dat1$Genus == "CHAMAESYCE" | dat1$Genus == "ALLIONIA" | dat1$Genus == "CALLIGONUM" | dat1$Genus == "PORTULACA" | dat1$Genus == "EUPHORBIA", "C4",
                       ifelse(dat1$Genus == "BELAPHARIS" | dat1$Genus == "AERVA" | dat1$Genus == "ALTERNANTHERA" | dat1$Genus == "ATRIPLEX" | dat1$Genus == "SUAEDA" | dat1$Genus == "TECTICORNIA" | dat1$Genus == "FLAVERIA" | dat1$Genus == "POLYCARPOREA" | dat1$Genus == "ELEOCHARS" | dat1$Genus == "RHYNCHOSPORA" | dat1$Genus == "EUPHORBIA" | dat1$Genus == "MOLLUGO" | dat1$Genus == "BOERHAVIA" | dat1$Genus == "BASSIA" | dat1$Family == "Poaceae", NA,
                              "C3"))
unique(dat1$ps_path)

dat1 <- dat1 %>%
  rename(ps_path2 = ps_path)

names(dat1)
names(cover)

cover <- cover %>% 
  left_join(dat1, by = c("year", "site_name", "site_code", "block", "plot", "subplot", "year_trt", "trt", 
                         "Family", "Taxon", "live", "local_provenance", "local_lifeform", "local_lifespan", 
                         "functional_group", "max_cover"))

cover <- cover %>%
  mutate(ps_path2 = if_else(is.na(ps_path2) & !is.na(ps_path), ps_path, ps_path2))

cover_by_site_plot <- cover %>%
  group_by(site_code, plot, trt) %>%
  summarise(
    total_cover = sum(max_cover, na.rm = TRUE),
    c4_cover = if (any(ps_path2 == "C4", na.rm = TRUE)) {
      sum(max_cover[ps_path2 == "C4"], na.rm = TRUE)} else {0},
    c4_proportion = c4_cover / total_cover,
    annual_cover = if (any(local_lifespan == "ANNUAL", na.rm = TRUE)) {
      sum(max_cover[local_lifespan == "ANNUAL"], na.rm = TRUE)} else {0},
    annual_proportion = annual_cover / total_cover,
    .groups = "drop"
  )

cover_by_site_trt <- cover_by_site_plot %>%
  group_by(site_code, trt) %>%
  summarise(
    avg_c4_proportion = mean(c4_proportion, na.rm = TRUE),
    avg_annual_proportion = mean(annual_proportion, na.rm = TRUE),
    .groups = "drop"
  )


### Covariate analysis of mass

mass_ppt_c_npk_edited <- mass_ppt_c_npk %>%
  dplyr::select(site_code, block, plot, continent, country, region, habitat, trt, year, 
                vascular_live_mass, log_mass, mswep_ppt, log_mswep_ppt, prev_ppt, year_trt, 
                proportion_par, avg_ppt_site, richness_vegan, MAT_v2, AI, PET)

unique(mass_ppt_c_npk_edited$site_code)

mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>% 
  left_join(lrr_df, by = "site_code")

mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>% 
  left_join(cover_by_site_trt, by = c("site_code", "trt"))

mass_ppt_c_npk_edited <- na.omit(mass_ppt_c_npk_edited)
unique(mass_ppt_c_npk_edited$site_code)

mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>% 
  mutate(AI2 = avg_ppt_site/PET)

full_model <- lmer(log_mass ~ trt * (log_mswep_ppt + proportion_par + avg_ppt_site + AI2 + MAT_v2 + richness_vegan 
                                     + prev_ppt + lrr_mass + avg_c4_proportion + avg_annual_proportion)
                   + (1 | site_code/year_trt), 
                   data = mass_ppt_c_npk_edited, REML = FALSE, na.action = "na.fail")

summary(full_model)
full_model_table <- dredge(full_model, m.lim=c(NA, 6), fixed = c("c.Control", "c.NPK"))
full_model_avg <- model.avg(get.models(full_model_table, subset = delta < 10))

summary(full_model_avg); sw(full_model_avg)


mass_map_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = avg_ppt_site, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MAP (mm)") + ylab("") +
  theme_bw(12)
mass_map_plot

mass_par_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = proportion_par, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Proportion PAR") + ylab("Bimass (g m-2)") +
  theme_bw(12)
mass_par_plot

mass_rich_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = richness_vegan, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Richness") + ylab("") +
  theme_bw(12)
mass_rich_plot

mass_prev_ppt_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = prev_ppt, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Previous years' precipitation (mm)") + ylab("") +
  theme_bw(12)
mass_prev_ppt_plot

mass_lrr_mass_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = lrr_mass, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw(12)
mass_lrr_mass_plot

mass_covar_figure <- ggarrange(mass_par_plot, mass_rich_plot, mass_lrr_mass_plot, mass_map_plot, mass_prev_ppt_plot, 
                               ncol = 5, common.legend = TRUE, legend = "bottom", align = 'hv')
mass_covar_figure

### Covariate analysis of R2 and slope

results_long <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_mswep_ppt"]
  npk_slope <- coef(npk_model)["log_mswep_ppt"]
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope
  ))
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope
  ))
}

averages <- mass_ppt_c_npk_edited %>%
  group_by(site_code, trt) %>%
  summarise(
    avg_proportion_par = mean(proportion_par, na.rm = TRUE),
    avg_avg_ppt_site = mean(avg_ppt_site, na.rm = TRUE),
    avg_richness = mean(richness_vegan, na.rm = TRUE),
    avg_lrr_mass = mean(lrr_mass, na.rm = TRUE),
    avg_avg_c4_proportion = mean(avg_c4_proportion, na.rm = TRUE),
    avg_avg_annual_proportion = mean(avg_annual_proportion, na.rm = TRUE),
    region = first(region)
  )

results_with_averages <- results_long %>%
  left_join(averages, by = c("site_code", "trt"))

full_r2_model <- lm(r2 ~ trt * (avg_proportion_par + avg_avg_ppt_site + avg_richness + avg_lrr_mass
                                + avg_avg_c4_proportion + avg_avg_annual_proportion), 
                    data = results_with_averages, na.action = "na.fail")
summary(full_r2_model)
model_set <- dredge(full_r2_model)
best_model_r2 <- get.models(model_set, 1)[[1]]
summary(best_model_r2)
r2_best_model_r2 <- performance::r2(best_model_r2)

slope_lrr_model <- lm(slope ~ avg_lrr_mass, data = results_with_averages)
summary(slope_lrr_model)

full_slope_model <- lm(slope ~ trt * (avg_proportion_par + avg_avg_ppt_site + avg_richness + avg_lrr_mass), 
                       data = results_with_averages, na.action = "na.fail")
summary(full_slope_model)
model_set <- dredge(full_slope_model)
best_model_slope <- get.models(model_set, 1)[[1]]
summary(best_model_slope)
r2_best_model_slope <- performance::r2(best_model_slope)

r2_map_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("MAP (mm)") + ylab("R2 of ppt vs. mass") +
  theme_bw(12)
r2_map_plot

r2_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("") +
  theme_bw(12)
r2_par_plot

r2_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Richness") + ylab("") +
  theme_bw(12)
r2_rich_plot

r2_lrr_mass_plot <- ggplot(data = results_with_averages, aes(x = avg_lrr_mass, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw(12)
r2_lrr_mass_plot

r2_covar_figure <- ggarrange(r2_map_plot, r2_par_plot, r2_rich_plot, r2_lrr_mass_plot,
                             ncol = 4, common.legend = TRUE, legend = "bottom", align = 'hv')
r2_covar_figure

slope_map_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("MAP") + ylab("Slope of ppt vs. mass") +
  theme_bw(12)
slope_map_plot

slope_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion Par") + ylab("") +
  theme_bw(12)
slope_par_plot

slope_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Richness") + ylab("") +
  theme_bw(12)
slope_rich_plot

slope_lrr_mass_plot <- ggplot(data = results_with_averages, aes(x = avg_lrr_mass, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw(12)
slope_lrr_mass_plot

slope_covar_figure <- ggarrange(slope_map_plot, slope_par_plot, slope_rich_plot, slope_lrr_mass_plot,
                             ncol = 4, common.legend = TRUE, legend = "bottom", align = 'hv')
slope_covar_figure


## testing for habitat effect

mass_aov <- aov(log_mass ~ habitat * trt, data = mass_ppt_c_npk)
summary(mass_aov)

r2_aov <- aov(r2 ~ habitat * trt, data = results_with_averages)
summary(r2_aov)
r2_emm <- emmeans(r2_aov, ~ trt | habitat)
pairwise_results <- pairs(r2_emm)
print(pairwise_results)

ggplot(results_with_averages, aes(x = trt, y = r2)) +
  geom_point() +
  facet_wrap(~ habitat) +
  theme_bw(14) +
  labs(x = "Habitat",
       y = "R2 of precipitation-mass")

slope_aov <- aov(slope ~ habitat * trt, data = results_with_averages)
summary(slope_aov)
slope_emm <- emmeans(slope_aov, ~ trt | habitat)
pairwise_results_slope <- pairs(slope_emm)
print(pairwise_results_slope)

ggplot(results_with_averages, aes(x = trt, y = slope)) +
  geom_point() +
  facet_wrap(~ habitat) +
  theme_bw(14) +
  labs(x = "Habitat",
       y = "Slope of precipitation-mass")


### Calculating and graphing variance in log_mass

variance <- mass_ppt_c_npk %>%
  group_by(trt) %>%
  summarise(total_variance = var(log_mass, na.rm = TRUE))

variance$conditional_r2 <- c(conditional_r2_control, conditional_r2_npk)
variance$marginal_r2 <- c(marginal_r2_control, marginal_r2_npk)

variance$prop_variance_conditional <- variance$total_variance * variance$conditional_r2
variance$prop_variance_marginal <- variance$total_variance * variance$marginal_r2

variance$unex_variance_conditional <- variance$total_variance - variance$prop_variance_conditional
variance$unex_variance_marginal <- variance$total_variance - variance$prop_variance_marginal

ggplot(variance, aes(x = trt)) +
  geom_bar(aes(y = total_variance), stat = "identity", fill = "darkgrey") +
  geom_bar(aes(y = prop_variance_conditional), stat = "identity", fill = "lightgrey") +
  labs(y = "Variance",
       x = "Treatment",
       fill = "Variance Type") +
  theme_bw(14)


### Additional graphs

# Graphing mean, r2, and slope on absolute scale (not log transformed)
results_graphing <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  control_model <- lm(vascular_live_mass ~ mswep_ppt, data = site_data_control)
  npk_model <- lm(vascular_live_mass ~ mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["mswep_ppt"]
  npk_slope <- coef(npk_model)["mswep_ppt"]
  control_mean <- mean(site_data_control$vascular_live_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$vascular_live_mass, na.rm = TRUE)
  results_graphing <- rbind(results_graphing, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean
  ))
  results_graphing <- rbind(results_graphing, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean
  ))
}

r2_boxplot <- ggplot(results_graphing, aes(x = trt, y = r2, color = trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "", y = "R2 of ppt vs. mass") +
  theme_bw(14)

slope_boxplot <- ggplot(results_graphing, aes(x = trt, y = slope, color = trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "", y = "Slope of ppt vs. mass") +
  theme_bw(14)

mean_boxplot <- ggplot(results_graphing, aes(x = trt, y = mean, color = trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(x = "", y = "Mean mass") +
  theme_bw(14)

boxplots <- ggarrange(mean_boxplot, slope_boxplot, r2_boxplot, 
                      ncol = 3, common.legend = TRUE, legend = "none", align = 'hv')
boxplots


mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>%
  left_join(results_long, by = c("site_code", "trt"))

ggplot(data = mass_ppt_c_npk_edited, aes(x = mswep_ppt, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("R2 of precipitation vs. biomass") +
  theme_bw(14)

ggplot(data = mass_ppt_c_npk_edited, aes(x = mswep_ppt, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Slope of precipitation vs. biomass") +
  theme_bw(14)


## Calculating and graphing effect sizes

mean_model2 <- lmer(log_mass ~ trt + (1 | site_code / year_trt), data = mass_ppt_c_npk)
summary(mean_model)

mean_model <- lmer(mean ~ trt + (1| site_code), data = results_graphing)
summary(mean_model2)

slope_model <- lmer(slope ~ trt + (1| site_code), data = results_graphing)
summary(slope_model)

r2_model <- lmer(r2 ~ trt + (1| site_code), data = results_graphing)
summary(r2_model)


mean_estimate <- 0.1961    
mean_se <- 0.008463        
mean_resid_sd <- 0.2065    
n_mean <- 2415             

slope_estimate <- 0.3570   
slope_se <- 0.1262         
slope_resid_sd <- 0.4968   
n_slope <- 62              

r2_estimate <- -0.006608   
r2_se <- 0.020012          
r2_resid_sd <- 0.07879     
n_r2 <- 62                 

calc_cohen_d <- function(estimate, resid_sd, n) {
  d <- estimate / resid_sd
  n1 <- n / 2
  n2 <- n / 2
  SE_d <- sqrt((n1 + n2) / (n1 * n2) + (d^2) / (2 * (n1 + n2))) 
  lower_CI <- d - 1.96 * SE_d
  upper_CI <- d + 1.96 * SE_d
  return(c(d, lower_CI, upper_CI))
}

mean_results <- calc_cohen_d(mean_estimate, mean_resid_sd, n_mean)
slope_results <- calc_cohen_d(slope_estimate, slope_resid_sd, n_slope)
r2_results <- calc_cohen_d(r2_estimate, r2_resid_sd, n_r2)

cohen_d_df <- data.frame(
  Variable = c("Mean", "Slope", "R²"),
  Cohen_d = c(mean_results[1], slope_results[1], r2_results[1]),
  Lower_CI = c(mean_results[2], slope_results[2], r2_results[2]),
  Upper_CI = c(mean_results[3], slope_results[3], r2_results[3])
)

cohen_d_df$Variable <- factor(cohen_d_df$Variable, levels = c("R²", "Slope", "Mean"))

ggplot(cohen_d_df, aes(x = Cohen_d, y = Variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2) +
  labs(x = "Effect Size (Cohen's d)",
       y = "") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(14) +
  theme(axis.text.y = element_text(size = 14))
  

unique(mass_ppt_c_npk$site_code)
