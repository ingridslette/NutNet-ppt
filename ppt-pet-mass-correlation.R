
## Digging into PPT-PET model in more detail

### Comparing control vs. NPK R2 - calculate and compare difference at each site

results_pet <- data.frame(site_code = character(),
                          control_r2 = numeric(), 
                          npk_r2 = numeric(),
                          r2_diff = numeric(),
                          control_slope = numeric(), 
                          npk_slope = numeric(),
                          slope_diff = numeric(),
                          stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(live_mass ~ ppt_pet, data = site_data_control)
  npk_model <- lm(live_mass ~ ppt_pet, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["ppt_pet"]
  npk_slope <- coef(npk_model)["ppt_pet"]
  results_pet <- rbind(results_pet, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_diff = npk_r2 - control_r2,
    control_slope = control_slope,
    npk_slope = npk_slope,
    slope_diff = npk_slope - control_slope
  ))
}

paired_t_test_r2 <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
print(paired_t_test_r2)


## Covariate model of mass

full_model_pet <- lmer(live_mass ~ trt * (ppt_pet + proportion_par + AI + rich + prev_ppt + lrr_mass
                                     + MAT_v2 + MAP_v2 + avg_c4_proportion + avg_annual_proportion)
                   + (1 | site_code/year_trt) + (1 | site_code/block), 
                   data = mass_ppt_edited, REML = FALSE, na.action = "na.fail")

summary(full_model_pet)
full_model_table_pet <- dredge(full_model_pet, m.lim = c(NA, 6), fixed = c("c.Control", "c.NPK"))
full_model_avg_pet <- model.avg(get.models(full_model_table_pet, subset = delta < 10))
summary(full_model_avg_pet); sw(full_model_avg_pet)


## Covariate models of R2 and slope

results_long_pet <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(live_mass ~ ppt_pet, data = site_data_control)
  npk_model <- lm(live_mass ~ ppt_pet, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["ppt_pet"]
  npk_slope <- coef(npk_model)["ppt_pet"]
  control_mean <- mean(site_data_control$log_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$log_mass, na.rm = TRUE)
  results_long_pet <- rbind(results_long_pet, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean
  ))
  results_long_pet <- rbind(results_long_pet, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean
  ))
}

results_with_averages_pet <- results_long_pet %>%
  left_join(averages, by = c("site_code", "trt"))

full_r2_model_pet <- lm(r2 ~ trt * (avg_proportion_par + avg_ai + avg_richness + avg_lrr_mass
                                + avg_avg_c4_proportion + avg_avg_annual_proportion), 
                    data = results_with_averages_pet, na.action = "na.fail")

summary(full_r2_model)


full_slope_model_pet <- lm(slope ~ trt * (avg_proportion_par  + avg_ai + avg_richness + avg_lrr_mass
                                      + avg_avg_c4_proportion + avg_avg_annual_proportion), 
                       data = results_with_averages_pet, na.action = "na.fail")

summary(full_slope_model_pet)


## Main graph

fig2_pet <- ggplot(mass_ppt_edited, aes(x = ppt_pet, y = live_mass)) +
  geom_smooth(method = lm, aes(color = site_code), se = F) +
  geom_smooth(method = lm, color = "black") +
  labs(x = "Annual GSP-PET(mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw()

fig2_pet

## Analyzing only data from the driest year at each site
### testing for convergence of RUE a la Huxman and Smith 2004

driest_year_mass_ppt_pet <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  slice_min(ppt_pet) %>%
  ungroup()

main_model_driest_pet <- lmer(live_mass ~ ppt_pet * trt + (1 | site_code / block) + (1 | year_trt), 
                          data = driest_year_mass_ppt_pet)

summary(main_model_driest_pet)

driest_year_plot_pet <- ggplot(data = driest_year_mass_ppt_pet, 
                           aes(x = ppt_pet, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", alpha = 0.3) +
  labs(x = "Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

driest_year_plot_pet


main_model_p05_pet <- lmer(live_mass ~ ppt_pet * trt + (1 | site_code / block) + (1 | year_trt), 
                       data = subset(mass_ppt_edited, ppt<p05_ppt))

summary(main_model_p05_pet)

main_model_p10_pet <- lmer(live_mass ~ ppt_pet * trt + (1 | site_code / block) + (1 | year_trt), 
                       data = subset(mass_ppt_edited, ppt<p10_ppt))

summary(main_model_p10_pet)

### Calculating and graphing effect sizes

mean_model_pet <- lmer(mean ~ trt + (1| site_code), data = results_long_pet)
summary(mean_model_pet)

slope_model_pet <- lmer(slope ~ trt + (1| site_code), data = results_long_pet)
summary(slope_model_pet)

r2_model_pet <- lmer(r2 ~ trt + (1| site_code), data = results_long_pet)
summary(r2_model_pet)

mean_estimate <- 0.17159
mean_se <- 0.01737
mean_resid_sd <- 0.07268
n_mean <- 70

slope_estimate <- 0.2718
slope_se <- 0.1041
slope_resid_sd <- 0.4354
n_slope <- 70   

r2_estimate <- -0.002222
r2_se <- 0.020551
r2_resid_sd <- 0.08597
n_r2 <- 70

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
  Variable = c("Biomass", "Sensitivity", "R²"),
  Cohen_d = c(mean_results[1], slope_results[1], r2_results[1]),
  Lower_CI = c(mean_results[2], slope_results[2], r2_results[2]),
  Upper_CI = c(mean_results[3], slope_results[3], r2_results[3])
)

cohen_d_df$Variable <- factor(cohen_d_df$Variable, levels = c("R²", "Sensitivity", "Biomass"))

es_fig <- ggplot(cohen_d_df, aes(x = Cohen_d, y = Variable, color = Variable)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), width = 0.2) +
  labs(x = "Effect Size",
       y = "") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#6F6F6F", "#ff924c", "#ff924c")) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 14), legend.position = "none")

es_fig

