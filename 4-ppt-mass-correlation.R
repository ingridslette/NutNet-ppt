library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(MuMIn)
library(performance)
library(MASS)
library(broom)
library(purrr)
library(emmeans)
library(ggpubr)

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/mswep_ppt_annual_gs_only.csv")

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

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/comb-by-plot-clim-soil-diversity_2024-05-31.csv",
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, year_trt > 0)
unique(mass1$year_trt)

mass_ppt <- inner_join(mass1, mswep, by=c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(vascular_live_mass),
         log_mswep_ppt = log10(mswep_ppt))

unique(mass_ppt$trt)

mass_ppt_c_npk <- filter(mass_ppt, trt %in% c("Control", "NPK")) 

unique(mass_ppt_c_npk$site_code)
unique(mass_ppt_c_npk$trt)

site_year_counts <- mass_ppt_c_npk %>%
  group_by(site_code, trt) %>%
  filter(!is.na(vascular_live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_6_years <- site_year_counts %>%
  filter(year_count >= 6) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) 

mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(site_code %in% sites_with_6_years$site_code)

unique(mass_ppt_c_npk$site_code)

# Filter to keep only sites with a certain ppt range
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  group_by(site_code) %>%
  filter(max(mswep_ppt_sd) > 1 & min(mswep_ppt_sd) < -1) %>%
  ungroup()

unique(mass_ppt_c_npk$site_code)

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Total live mass") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x=year_trt, y=vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Treatment Year") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt, y= vascular_live_mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation",
       y = "Total live mass",
       color = "Site Code") +
  theme(legend.position = "right")


# Graphing back-transformed data

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
  labs(x = "Total Growing Season Precipitation (mm)", y = "Live Mass") +
  facet_wrap(~ trt) +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x = mswep_ppt, y = vascular_live_mass, color = trt)) +
  geom_point() +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Total Growing Season Precipitation (mm)", y = "Live Mass") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw()

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
  labs(x = "Growing Season Precipitation (mm)", y = "Live Mass") +
  facet_wrap(~ trt) +
  theme_bw()

ggplot(data = mass_ppt_c_npk,aes(x= mswep_ppt, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + 
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("live mass") +
  theme_bw()


### Comparing control vs. NPK R2 - Approach 1: calculate and compare difference at each site

# Get unique site codes
site_codes <- unique(mass_ppt_c_npk$site_code)

# Initialize a dataframe to store results
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

# t-test on the r2_difference values
t_test_r2_diff <- t.test(results$r2_difference)
print(t_test_r2_diff)

# paired t-test on the R2 values for Control and NPK
paired_t_test_result <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
print(paired_t_test_result)

# sites with higher R2 in NPK plots
sites_higher_r2_npk <- filter(results, r2_difference < 0)


### Comparing control vs. NPK R2 - Approach 2: fit separate models for control and NPK data, calculate and compare z scores

model_control <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                      data = subset(mass_ppt_c_npk, trt == "Control"))
model_npk <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                  data = subset(mass_ppt_c_npk, trt == "NPK"))

AIC(model_control, model_npk)

r2_control <- performance::r2(model_control)
r2_npk <- performance::r2(model_npk)

conditional_r2_control <- r2_control$R2_conditional
conditional_r2_npk <- r2_npk$R2_conditional

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

# Get 95% confidence intervals for the R² difference
boot_ci <- boot.ci(boot_r2, type = "perc")
print(boot_ci)


## Updating graphs to order by R2

predictions <- predictions %>%
  left_join(results, by = "site_code")

predictions <- predictions %>%
  mutate(site_code = factor(site_code, levels = unique(site_code[order(r2_difference)])))

pal2 <- c("#800000","#c00000","#ff0000","#ff4040","#ff8080","#a83a01","#e04d01","#f06201","#ff7700","#e0a500",
          "#ffbc00","#ffcd40","#ffde80","#305020","#406a2a","#609f3f","#80d353","#bdda0f","#0b5043","#117864",
          "#1abc9c","#03045e","#125d93","#057dcd","#43b0f1","#96cff1","#503658","#80558c","#af7ab3","#c497b0",
          "#808080")

ggplot(predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass, colour = site_code)) +
  geom_line() +
  scale_color_manual(values = pal2) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            color = "black", linetype = "dashed") +
  labs(x = "Growing Season Precipitation (mm)", y = "Live Mass") +
  facet_wrap(~ trt) +
  theme_bw()


## Calculating log response ratios

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


## Covariate analysis of mass

c_npk_x_model <- lmer(log_mass ~ log_mswep_ppt * trt + (1 | site_code / year_trt), data = mass_ppt_c_npk)
summary(c_npk_x_model)

mass_ppt_c_npk_edited <- mass_ppt_c_npk %>%
  dplyr::select(site_code, block, plot, continent, country, region, habitat, trt, year, 
                vascular_live_mass,log_mass, mswep_ppt,log_mswep_ppt, prev_ppt, year_trt, 
                proportion_par, avg_ppt_site, richness_vegan, litter_mass)

#mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>% group_by(site_code) %>% mutate(PercentSand = if_else(is.na(PercentSand), mean(PercentSand, na.rm = TRUE), PercentSand)) %>% ungroup()

unique(mass_ppt_c_npk_edited$site_code)

mass_ppt_c_npk_edited <- mass_ppt_c_npk_edited %>% 
  left_join(lrr_df, by = "site_code")

mass_ppt_c_npk_edited <- na.omit(mass_ppt_c_npk_edited)
unique(mass_ppt_c_npk_edited$site_code)

full_model_x <- lmer(log_mass ~ trt * (log_mswep_ppt + proportion_par + avg_ppt_site 
                                     + richness_vegan + prev_ppt + lrr_mass) + (1 | site_code/year_trt), 
                     data = mass_ppt_c_npk_edited, REML = FALSE, na.action = "na.fail")
summary(full_model_x)
model_set <- dredge(full_model_x)
best_model <- get.models(model_set, 1)[[1]]
summary(best_model)

#mass_ppt_c <- subset(mass_ppt_c_npk_edited, trt == 'Control')
#mass_ppt_npk <- subset(mass_ppt_c_npk_edited, trt == 'NPK')


mass_map_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = avg_ppt_site, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("") + ylab("Live mass") +
  theme_bw()
mass_map_plot

mass_sand_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = PercentSand, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("") + ylab("") +
  theme_bw()
mass_sand_plot

mass_par_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = proportion_par, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("") + ylab("") +
  theme_bw()
mass_par_plot

mass_rich_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = richness_vegan, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("") + ylab("") +
  theme_bw()
mass_rich_plot

mass_prev_ppt_plot <- ggplot(data = mass_ppt_c_npk_edited, aes(x = prev_ppt, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("") + ylab("") +
  theme_bw()
mass_prev_ppt_plot


## Covariate analysis of R2 and slope

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
    region = first(region),
    country = first(country),
    habitat = first(habitat)
  )

results_with_averages <- results_long %>%
  left_join(averages, by = c("site_code", "trt"))

#results_with_averages <- results_with_averages %>% left_join(lrr_df, by = "site_code")

results_with_averages_edited <- na.omit(results_with_averages)

unique(results_with_averages_edited$site_code)

full_r2_model_x <- lm(r2 ~ trt * (avg_proportion_par + avg_avg_ppt_site + avg_lrr_mass + avg_richness), 
                      data = results_with_averages_edited, na.action = "na.fail")
summary(full_r2_model_x)
model_set <- dredge(full_r2_model_x)
best_model <- get.models(model_set, 1)[[1]]
summary(best_model)

#results_with_averages_c <- subset(results_with_averages, trt == 'Control')
#results_with_averages_npk <- subset(results_with_averages, trt == 'NPK')

full_slope_model_x <- lm(slope ~ trt * (avg_proportion_par + avg_avg_ppt_site + avg_lrr_mass + avg_richness), 
                   data = results_with_averages_edited, na.action = "na.fail")
summary(full_slope_model_x)
model_set <- dredge(full_slope_model_x)
best_model <- get.models(model_set, 1)[[1]]
summary(best_model)


ggplot(data = mass_ppt_c_npk, aes(x = proportion_par, y = vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("Live Mass") +
  theme_bw()

ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = avg_proportion_par, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("MAP") + ylab("Proportion PAR") +
  theme_bw()

ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = avg_PercentSand, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("MAP") + ylab("Percent Sand") +
  theme_bw()

ggplot(results_with_averages, aes(x = lrr_mass, y = slope, color = trt)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Slope of ppt vs. mass",
       y = "LRR mass") +
  theme_bw()

ggplot(results, aes(x = LRR, y = r2_difference)) +
  geom_point() +
  geom_smooth(method = "lm", color = "darkgrey") +
  labs(x = "LRR mass",
       y = "Difference in R2 of ppt vs. mass") +
  theme_bw(14)

r2_map_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("") + ylab("R2 of ppt vs. mass") +
  theme_bw()
r2_map_plot

r2_sand_plot <- ggplot(data = results_with_averages, aes(x = avg_PercentSand, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("") + ylab("") +
  theme_bw()
r2_sand_plot

r2_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("") + ylab("") +
  theme_bw()
r2_par_plot

r2_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("") + ylab("") +
  theme_bw()
r2_rich_plot

slope_map_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_ppt_site, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("MAP") + ylab("Slope of ppt vs. mass") +
  theme_bw()
slope_map_plot

slope_sand_plot <- ggplot(data = results_with_averages, aes(x = avg_PercentSand, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Percent Sand") + ylab("") +
  theme_bw()
slope_sand_plot

slope_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion Par") + ylab("") +
  theme_bw()
slope_par_plot

slope_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Richness") + ylab("") +
  theme_bw()
slope_rich_plot

figure <- ggarrange(mass_map_plot, mass_sand_plot, mass_par_plot, mass_rich_plot,
                    r2_map_plot, r2_sand_plot, r2_par_plot, r2_rich_plot,
                    slope_map_plot, slope_sand_plot, slope_par_plot, slope_rich_plot,
                    ncol = 4, nrow = 3, common.legend = TRUE, legend = "bottom", align = 'hv')
figure

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
  labs(x = "Treatment",
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
  labs(x = "Treatment",
       y = "Slope of precipitation-mass")

## testing for regional effect

mass_aov_reg <- aov(log_mass ~ country * trt, data = mass_ppt_c_npk)
summary(mass_aov_reg)

r2_aov_reg <- aov(r2 ~ country * trt, data = results_with_averages)
summary(r2_aov_reg)
r2_emm_reg <- emmeans(r2_aov_reg, ~ trt | region)
pairwise_results_reg <- pairs(r2_emm_reg)
print(pairwise_results_reg)

ggplot(results_with_averages, aes(x = trt, y = r2)) +
  geom_point() +
  facet_wrap(~ region) +
  theme_bw() +
  labs(x = "Treatment",
       y = "R2 of precipitation-mass")

slope_aov_reg <- aov(slope ~ region * trt, data = results_with_averages)
summary(slope_aov_reg)
slope_emm_reg <- emmeans(slope_aov_reg, ~ trt | region)
pairwise_results_slope_reg <- pairs(slope_emm_reg)
print(pairwise_results_slope_reg)

ggplot(results_with_averages, aes(x = trt, y = slope)) +
  geom_point() +
  facet_wrap(~ region) +
  theme_bw() +
  labs(x = "Treatment",
       y = "Slope of precipitation-mass")





