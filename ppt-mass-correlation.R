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
library(cowplot)


### Loading, viewing, and filtering precipitation and mass data 

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/comb-by-plot-clim-soil-diversity_2024-05-31.csv",
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, year_trt > 0)
unique(mass1$year_trt)

mass1 <- filter(mass1, trt %in% c("Control", "NPK")) 
unique(mass1$trt)

mass1 <- mass1 %>%
  mutate(
    live_mass = case_when(
      !is.na(vascular_live_mass) | !is.na(nonvascular_live_mass) ~ 
        rowSums(across(c(vascular_live_mass, nonvascular_live_mass, standing_dead_mass)), na.rm = TRUE),
      is.na(vascular_live_mass) & is.na(nonvascular_live_mass) ~ 
        unsorted_live_mass
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

## popped over to script "daily-to-gs-ppt.R" here, to get growing season ppt for the sites included in mass2
## exported that as csv and now loading it here

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/mswep_ppt_annual_gs_only_2025-05-07.csv")

unique(mswep$site_code)

mswep <- mswep %>%
  arrange(site_code, year) %>% 
  group_by(site_code) %>%
  mutate(prev_ppt = lag(mswep_ppt)) %>%
  ungroup()

mswep <- filter(mswep, year >= 1983)

mswep <- filter(mswep, year < 2025) 
unique(mswep$year)

mswep <- mswep %>%
  group_by(site_code) %>%
  mutate(avg_ppt = mean(mswep_ppt, na.rm = TRUE),
         sd_ppt = sd(mswep_ppt, na.rm = TRUE)) %>%
  ungroup()

unique(mass2$site_code)
unique(mswep$site_code)

mass_ppt <- inner_join(mass2, mswep, by = c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(live_mass),
         log_mswep_ppt = log10(mswep_ppt))

mass_ppt <- mass_ppt %>%
  group_by(site_code) %>%
  mutate(min_ppt = min(mswep_ppt, na.rm = TRUE),
         max_ppt = max(mswep_ppt, na.rm = TRUE)) %>%
  ungroup()

# Filter to keep only sites with an observed ppt range that spans at least +- 1 sd of long-term avg
mass_ppt <- mass_ppt %>%
  group_by(site_code) %>%
  filter(min_ppt <= (avg_ppt - sd_ppt), max_ppt >= (avg_ppt + sd_ppt)) %>%
  ungroup()

unique(mass_ppt$site_code)


### Initial model

initial_model <- lmer(log_mass ~ log_mswep_ppt * trt + (1 | site_code / block) + (1 | year_trt), data = mass_ppt)
summary(initial_model)

# Model assumptions check 
plot(initial_model)
resid <- residuals(initial_model)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid)
plot(fitted(initial_model), resid, main = "Residuals vs Fitted")


### Initial graphs
ggplot(data = mass_ppt, aes(x = mswep_ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

ggplot(data = mass_ppt, aes(x = mswep_ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

ggplot(data = mass_ppt, aes(x = year_trt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Treatment Year") + ylab("Biomass (g/m2)") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

ggplot(data = mass_ppt, aes(x = mswep_ppt, y = live_mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation (mm)",
       y = "Biomass (g/m2)",
       color = "Site Code") +
  theme(legend.position = "right")


### Graphing back-transformed data, to allow non-linear curves on linear scale

# back transform from log-log scale
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

predictions <- mass_ppt %>%
  group_by(site_code, trt) %>%
  group_modify(~ fit_model_and_predict(.x)) %>%
  ungroup()

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

predictions_allsites <- mass_ppt %>%
  group_by(trt) %>%
  group_modify(~ fit_model_and_predict_allsites(.x)) %>%
  ungroup()

ggplot(data = mass_ppt, aes(x = mswep_ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + 
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  labs(color = "Treatment", shape = "Treatment") +
  scale_color_manual(values = c("#4267ac", "#ff924c")) +
  theme_bw(14)

ggplot(mass_ppt, aes(x = mswep_ppt, y = live_mass, color = trt)) +
  geom_point() +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)", color = "Treatment") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

ggplot(mass_ppt, aes(x = mswep_ppt, y = live_mass, color = site_code)) +
  geom_line(data = predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass), linewidth = 1) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            linewidth = 1, color = "black") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14)

pal2 <- c("#800000","#c00000","#ff0000","#ff4040","#ff8080","#a83a01","#e04d01","#f06201","#ff7700","#e0a500",
          "#ffbc00","#ffcd40","#ffde80","#305020","#406a2a","#609f3f","#80d353","#bdda0f","#0b5043","#117864",
          "#1abc9c","#03045e","#125d93","#057dcd","#43b0f1","#96cff1","#503658","#80558c","#af7ab3","#c497b0",
          "#808080")

pal3 <- c("#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080",
          "#808080")

palette_31 <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#a6cee3", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
  "#b3de69", "#fccde5", "#bc80bd", "#ccebc5", "#d9d9d9", 
  "#ffed6f"
)

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

figure1 <- ggplot(predictions, aes(x = 10^log_mswep_ppt, y = predicted_mass, colour = site_code)) +
  geom_line() +
  scale_color_manual(values = palette_31) +
  geom_line(data = predictions_allsites, aes(x = 10^log_mswep_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14) +
  theme(legend.position = "none")
figure1


### Comparing control vs. NPK R2 - Approach 1: calculate and compare difference at each site

site_codes <- unique(mass_ppt$site_code)

results <- data.frame(site_code = character(), 
                      control_r2 = numeric(), 
                      npk_r2 = numeric(),
                      control_slope = numeric(), 
                      npk_slope = numeric(),
                      stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_mswep_ppt"]
  npk_slope <- coef(npk_model)["log_mswep_ppt"]
  results <- rbind(results, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    control_slope = control_slope,
    npk_slope = npk_slope
  ))
}

paired_t_test_r2 <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
print(paired_t_test_r2)

paired_t_test_slope <- t.test(results$control_slope, results$npk_slope, paired = TRUE)
print(paired_t_test_slope)


### Comparing control vs. NPK R2 - Approach 2: fit separate models for control and NPK data, calculate and compare z scores

model_control <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                      data = subset(mass_ppt, trt == "Control"))
model_npk <- lmer(log_mass ~ log_mswep_ppt + (1 | site_code / year_trt), 
                  data = subset(mass_ppt, trt == "NPK"))

summary(model_control)
summary(model_npk)

AIC(model_control, model_npk)

r2_control <- performance::r2(model_control)
r2_npk <- performance::r2(model_npk)

conditional_r2_control <- r2_control$R2_conditional
conditional_r2_npk <- r2_npk$R2_conditional

marginal_r2_control <- r2_control$R2_marginal
marginal_r2_npk <- r2_npk$R2_marginal

# Compare R2 values using Fisher's Z transformation
z_control <- 0.5 * log((1 + sqrt(marginal_r2_control)) / (1 - sqrt(marginal_r2_control)))
z_npk <- 0.5 * log((1 + sqrt(marginal_r2_npk)) / (1 - sqrt(marginal_r2_npk)))

n_control <- length(unique(subset(mass_ppt, trt == "Control")$site_code))
n_npk <- length(unique(subset(mass_ppt, trt == "NPK")$site_code))
se_diff <- sqrt((1 / (n_control - 3)) + (1 / (n_npk - 3)))

# Calculate the Z-score for the difference
z_diff <- (z_control - z_npk) / se_diff

p_value <- 2 * (1 - pnorm(abs(z_diff)))

cat("Z-score for the difference:", z_diff, "\n")
cat("P-value for the difference in marginal R-squared values:", p_value, "\n")


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

boot_r2 <- boot(data = mass_ppt, statistic = r2_diff, R = 1000)

print(boot_r2)

boot_ci <- boot.ci(boot_r2, type = "perc")
print(boot_ci)


### Covariate analyses

## Incorporating log response ratio of mass to trt
lrr_df <- mass_ppt %>%
  group_by(site_code) %>%
  summarize(
    lrr_mass = log(mean(live_mass[trt == "NPK"], na.rm = TRUE) /
                     mean(live_mass[trt == "Control"], na.rm = TRUE))
  )

## Incorporating C3/C4 and annual/perennial information from cover data
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

mass_ppt_edited <- mass_ppt %>%
  dplyr::select(site_code, block, plot, continent, country, region, habitat, trt, year, 
                live_mass, log_mass, mswep_ppt, log_mswep_ppt, prev_ppt, year_trt, 
                proportion_par, avg_ppt, rich, MAT_v2, AI, PET)

unique(mass_ppt_edited$site_code)

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(lrr_df, by = "site_code")

mass_ppt_edited <- mass_ppt_edited %>% 
  left_join(cover_by_site_trt, by = c("site_code", "trt"))

mass_ppt_edited <- na.omit(mass_ppt_edited)
unique(mass_ppt_edited$site_code)

## covariate model of mass
full_model <- lmer(log_mass ~ trt * (log_mswep_ppt + proportion_par + AI + rich + prev_ppt 
                                     + lrr_mass + avg_c4_proportion + avg_annual_proportion)
                   + (1 | site_code/year_trt) + (1 | site_code/block), 
                   data = mass_ppt_edited, REML = FALSE, na.action = "na.fail")

summary(full_model)
full_model_table <- dredge(full_model, m.lim=c(NA, 6), fixed = c("c.Control", "c.NPK"))
full_model_avg <- model.avg(get.models(full_model_table, subset = delta < 10))
summary(full_model_avg); sw(full_model_avg)

mass_lrr_mass_plot <- ggplot(data = mass_ppt_edited, aes(x = lrr_mass, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_par_plot <- ggplot(data = mass_ppt_edited, aes(x = proportion_par, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Proportion PAR") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_ai_plot <- ggplot(data = mass_ppt_edited, aes(x = AI, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Aridity Index") + ylab("Biomass (g m-2)") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_rich_plot <- ggplot(data = mass_ppt_edited, aes(x = rich, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Richness") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_prev_ppt_plot <- ggplot(data = mass_ppt_edited, aes(x = prev_ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Previous years' precipitation (mm)") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_c4_plot <- ggplot(data = mass_ppt_edited, aes(x = avg_c4_proportion, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Proportion C4 Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

mass_annual_plot <- ggplot(data = mass_ppt_edited, aes(x = avg_annual_proportion, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Proportion Annual Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))


mass_covar_figure <- ggarrange(mass_ai_plot, mass_lrr_mass_plot, mass_par_plot, mass_prev_ppt_plot,
                               mass_rich_plot, mass_c4_plot, mass_annual_plot,
                               ncol = 4, nrow = 2, common.legend = TRUE, legend = "bottom", align = 'hv')
mass_covar_figure

## covariate models of R2 and slope
results_long <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_mswep_ppt"]
  npk_slope <- coef(npk_model)["log_mswep_ppt"]
  control_mean <- mean(site_data_control$live_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$live_mass, na.rm = TRUE)
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean
  ))
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean
  ))
}

averages <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  summarise(
    avg_proportion_par = mean(proportion_par, na.rm = TRUE),
    avg_avg_ppt = mean(avg_ppt, na.rm = TRUE),
    avg_mat = mean(MAT_v2, na.rm = TRUE),
    avg_richness = mean(rich, na.rm = TRUE),
    avg_lrr_mass = mean(lrr_mass, na.rm = TRUE),
    avg_avg_c4_proportion = mean(avg_c4_proportion, na.rm = TRUE),
    avg_avg_annual_proportion = mean(avg_annual_proportion, na.rm = TRUE),
    avg_ai = mean(AI, na.rm = TRUE)
  )

results_with_averages <- results_long %>%
  left_join(averages, by = c("site_code", "trt"))

full_r2_model <- lm(r2 ~ trt * (avg_proportion_par + avg_ai + avg_richness + avg_lrr_mass
                                + avg_avg_c4_proportion + avg_avg_annual_proportion), 
                    data = results_with_averages, na.action = "na.fail")
summary(full_r2_model)

full_slope_model <- lm(slope ~ trt * (avg_proportion_par  + avg_ai + avg_richness + avg_lrr_mass
                                      + avg_avg_c4_proportion + avg_avg_annual_proportion), 
                       data = results_with_averages, na.action = "na.fail")
summary(full_slope_model)

ai_slope_model <- lm(slope ~ trt * avg_ai, data = results_with_averages)
summary(ai_slope_model)

ai_slope_model_quad <- lm(slope ~ trt * poly(avg_ai, 2, raw = TRUE), data = results_with_averages)
summary(ai_slope_model_quad)


r2_lrr_mass_plot <- ggplot(data = results_with_averages, aes(x = avg_lrr_mass, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_ai_plot <- ggplot(data = results_with_averages, aes(x = avg_ai, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Aridity Index") + ylab("R2 of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Richness") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_c4_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_c4_proportion, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion C4 Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_annual_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_annual_proportion, y = r2, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion Annual Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

r2_covar_figure <- ggarrange(r2_ai_plot, r2_lrr_mass_plot, r2_par_plot,  r2_rich_plot, r2_c4_plot, r2_annual_plot,
                             ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", align = 'hv')
r2_covar_figure


slope_lrr_mass_plot <- ggplot(data = results_with_averages, aes(x = avg_lrr_mass, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

slope_par_plot <- ggplot(data = results_with_averages, aes(x = avg_proportion_par, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

slope_ai_plot <- ggplot(data = results_with_averages, aes(x = avg_ai, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Aridity Index") + ylab("slope") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))
slope_ai_plot

slope_ai_plot_quad <- ggplot(data = results_with_averages, aes(x = avg_ai, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, , formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
  xlab("Aridity Index") + ylab("Slope of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))
slope_ai_plot_quad

slope_rich_plot <- ggplot(data = results_with_averages, aes(x = avg_richness, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Richness") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

slope_c4_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_c4_proportion, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion C4 Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

slope_annual_plot <- ggplot(data = results_with_averages, aes(x = avg_avg_annual_proportion, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion Annual Species") + ylab("") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))


slope_covar_figure <- ggarrange(slope_ai_plot_quad, slope_lrr_mass_plot, slope_par_plot, 
                                slope_rich_plot, slope_c4_plot, slope_annual_plot,
                                ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", align = 'hv')
slope_covar_figure


## Looking at site-level trade-offs in magnitude of response to limiting factors

par_lrr_mass_plot <- ggplot(data = results_with_averages, 
                            aes(x = avg_lrr_mass, y = avg_proportion_par, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("Proportion PAR") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))
par_lrr_mass_plot

slope_lrr_mass_plot2 <- ggplot(data = results_with_averages, 
                               aes(x = avg_lrr_mass, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("Slope of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

slope_par_plot2 <- ggplot(data = results_with_averages, 
                         aes(x = avg_proportion_par, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("Slope of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

par_lrr_mass_model <- lm(avg_proportion_par ~ trt * avg_lrr_mass, data = results_with_averages, na.action = "na.fail")
summary(par_lrr_mass_model) # NS

slope_lrr_mass_model <- lm(slope ~ trt * avg_lrr_mass, data = results_with_averages, na.action = "na.fail")
summary(slope_lrr_mass_model) # NS

slope_par_model <- lm(slope ~ trt * avg_proportion_par, data = results_with_averages, na.action = "na.fail")
summary(slope_par_model) # NS

colimitation_figure <- ggarrange(slope_lrr_mass_plot2, slope_par_plot2, par_lrr_mass_plot,
                                 ncol = 1, common.legend = TRUE, legend = "bottom", align = 'hv')
colimitation_figure


### Calculating and graphing effect sizes

results_graphing <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(live_mass ~ mswep_ppt, data = site_data_control)
  npk_model <- lm(live_mass ~ mswep_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["mswep_ppt"]
  npk_slope <- coef(npk_model)["mswep_ppt"]
  control_mean <- mean(site_data_control$live_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$live_mass, na.rm = TRUE)
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

mean_model <- lmer(mean ~ trt + (1| site_code), data = results_graphing)
summary(mean_model)

slope_model <- lmer(slope ~ trt + (1| site_code), data = results_graphing)
summary(slope_model)

r2_model <- lmer(r2 ~ trt + (1| site_code), data = results_graphing)
summary(r2_model)


mean_estimate <- 191.32
mean_se <- 22.75
mean_resid_sd <- 89.55
n_mean <- 62

slope_estimate <- 0.3336
slope_se <- 0.1338
slope_resid_sd <- 0.5267
n_slope <- 62   

r2_estimate <- -0.01054
r2_se <- 0.02046
r2_resid_sd <- 0.080554
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


## Calculating and testing trt effect on RUE 

mass_ppt_edited <- mass_ppt_edited %>% 
  mutate(rue = live_mass/mswep_ppt)

rue_trt_model <- lmer(rue ~ trt + (1 | site_code / block) + (1 | year_trt), data = mass_ppt_edited)
summary(rue_trt_model)

ggplot(mass_ppt_edited, aes(x = trt, y = rue)) +
  geom_boxplot() +
  labs(x= "Treatment", y= "Rain Use Efficiency")
  theme_bw(14)

mass_ppt_edited_subset <- mass_ppt_edited %>%
  filter(rue < 2)

rue_trt_model_subset <- lmer(rue ~ trt + (1 | site_code / block) + (1 | year_trt), data = mass_ppt_edited_subset)
summary(rue_trt_model_subset)

ggplot(mass_ppt_edited_subset, aes(x = trt, y = rue)) +
  geom_boxplot() +
  labs(x= "Treatment", y= "Rain Use Efficiency")
theme_bw(14)


## Analyzing only data from the driest year at each site
### testing for convergence of RUE a la Huxman and Smith 2004

driest_year_mass_ppt <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  slice_min(mswep_ppt) %>%
  ungroup()

driest_year_plot <- ggplot(data = driest_year_mass_ppt, 
                           aes(x = mswep_ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  labs(x = "Precipitation", y = "Biomass", color = "Treatment", shape = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))
driest_year_plot

initial_model_driest <- lmer(log_mass ~ log_mswep_ppt * trt + (1 | site_code / block) + (1 | year_trt), data = driest_year_mass_ppt)
summary(initial_model_driest)


