
# "just playing around" code

fig2_4 <- ggarrange(
  es_fig,
  fig2_control + rremove("xlab"),
  fig2_npk + rremove("ylab") + rremove("xlab") +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ),
  slope_ai_plot_quad,
  legend = "bottom",
  common.legend = TRUE,
  ncol = 4,
  align = 'hv'
)

fig2_4 <- annotate_figure(
  fig2_4,
  bottom = text_grob("Growing Season Precipitation (mm)", size = 14)
)

fig2_4

fig2_both_ef <- ggarrange(es_fig, fig2_both,
                          widths = c(0.5, 1))

fig2_both_ef

my_palette <- colorRampPalette(c("#a50026","#d73027","#f46d43","#fdae61",
                                 "#CAF0F8","#90E0EF","#00A5D0","#0077B6",
                                 "#023E8A","#032174","#030455"))


my_palette <- colorRampPalette(c("#533C88","#7251B5","#B185DB","#D8C0E7",
                                 "#99E2B4","#78C6A3","#56AB91","#358F80",
                                 "#14746F","#116460"))

ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2)) +
  geom_line(aes(group = site_code)) +
  scale_color_gradientn(
    colors = my_palette(100),
    limits = c(-0.25, 0.4),
    name = "R²") +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 0.75) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

fig2 <- ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) +
  scale_color_gradientn(
    colors = my_palette2(50),
    limits = c(-0.25, 0.4),
    name = "Δ R²") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt, ncol = 1) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")
fig2

fig2_inset <- ggplot(data = mass_ppt, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
  geom_ribbon(data = predictions_allsites, 
              aes(x = 10^log_ppt, ymin = mass_lower, ymax = mass_upper, fill = trt),
              inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = predictions_allsites, 
            aes(x = 10^log_ppt, y = predicted_mass),
            linewidth = 1) +
  labs(x = "GSP (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  scale_color_manual(values = c("#4267ac", "#ff924c")) +
  scale_fill_manual(values = c("#4267ac", "#ff924c")) +
  theme_bw(16) +
  theme(legend.title = element_blank(), legend.position = "bottom")
fig2_inset

fig2_with_inset <- ggdraw() + 
  draw_plot(fig2) +
  draw_plot(fig2_inset, x = 0.09, y = 0.52, width = 0.22, height = 0.38)
fig2_with_inset


summary_stats <- mass_ppt %>%
  group_by(trt) %>%
  summarise(
    total_variance = var(log_mass, na.rm = TRUE),
    r2 = summary(lm(log_mass ~ log_ppt, data = cur_data()))$r.squared,
    exp_variance = total_variance * r2
  )

ggplot(summary_stats, aes(x = trt, fill = trt)) +
  geom_bar(aes(y = total_variance), stat = "identity", alpha = 0.3) +
  geom_bar(aes(y = exp_variance), stat = "identity") +
  theme_bw()



### Calculating and graphing variance in log_mass
results_long <- data.frame(site_code = character(), 
                           trt = character(), 
                           r2 = numeric(), 
                           slope = numeric(),
                           mean = numeric(),
                           variance = numeric(),
                           stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt, site_code == site & trt == "NPK")
  control_model <- lm(log_mass ~ log_ppt, data = site_data_control)
  npk_model <- lm(log_mass ~ log_ppt, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  control_slope <- coef(control_model)["log_ppt"]
  npk_slope <- coef(npk_model)["log_ppt"]
  control_mean <- mean(site_data_control$log_mass, na.rm = TRUE)
  npk_mean <- mean(site_data_npk$log_mass, na.rm = TRUE)
  control_variance <- var(site_data_control$log_mass, na.rm = TRUE)
  npk_variance <- var(site_data_npk$log_mass, na.rm = TRUE)
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "Control",
    r2 = control_r2,
    slope = control_slope,
    mean = control_mean,
    variance = control_variance
  ))
  results_long <- rbind(results_long, data.frame(
    site_code = site,
    trt = "NPK",
    r2 = npk_r2,
    slope = npk_slope,
    mean = npk_mean,
    variance = npk_variance
  ))
}

variance_summary <- results_long %>%
  group_by(trt) %>%
  summarise(
    mean_variance = mean(variance, na.rm = TRUE),
    se_variance = sd(variance, na.rm = TRUE) / sqrt(n()),
    mean_r2 = mean(r2, na.rm = TRUE),
    se_r2 = sd(r2, na.rm = TRUE) / sqrt(n()),
  )


ggplot(variance_summary, aes(x = trt), fill = trt) +
  geom_bar(aes(y = mean_variance), stat = "identity") +
  geom_bar(aes(y = mean_r2), stat = "identity", alpha = 0.3) +
  labs(y = "Variance",
       x = "Treatment") +
  theme_bw()

variance <- mass_ppt_edited %>%
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
  geom_bar(aes(y = unex_variance_marginal), stat = "identity", fill = "lightgrey") +
  labs(y = "Variance",
       x = "Treatment",
       fill = "Variance Type") +
  theme_bw()


# Graphing mean, r2, and slope on absolute scale (not log transformed)
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


mass1 %>%
  filter(is.na(vascular_live_mass) & is.na(nonvascular_live_mass) & !is.na(unsorted_live_mass)) %>%
  count(site_code)

mass1 %>%
  filter(!is.na(standing_dead_mass)) %>%
  count(site_code)


fig1_inset <- ggplot(data = results_with_averages, aes(x = avg_ai, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, , formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
  labs(x = "Aridity Index", y = "Sensitivity", color = "Treatment", shape = "Treatment") +
  theme_bw(12) +
  theme(legend.title = element_blank(), legend.position = c(0.77, 0.85), 
        legend.background = element_rect(fill = alpha("white", 0))) +
  scale_color_manual(values = c("#BFBFBF", "#666666"))

figure1_with_inset <- ggdraw() + 
  draw_plot(figure1) +
  draw_plot(fig1_inset, x = 0.08, y = 0.55, width = 0.22, height = 0.38)
figure1_with_inset


### Initial graphs
ggplot(data = mass_ppt, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  theme_bw() +
  scale_color_manual(values = c("#4267ac", "#ff924c"))

ggplot(data = mass_ppt, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
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

ggplot(data = mass_ppt, aes(x = ppt, y = live_mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation (mm)",
       y = "Biomass (g/m2)",
       color = "Site Code") +
  theme(legend.position = "right")


ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) +  # group applied here, where site_code exists
  scale_color_gradient2(low = "#a50026", mid = '#ffffbf', high = "#313695", 
                        midpoint = 0, 
                        limits = c(-0.1, 0.1),
                        oob = scales::squish) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", color = "Δ R²") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")


ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) + 
  scale_color_gradientn(
    colours = c("#a50026",
      "#d73027",
      "#f46d43",
      "#fdae61",
      "grey",
      "#d1e5f0",
      "#92c5de",
      "#4393c3",
      "#2166ac",
      "#053061"),
    values = scales::rescale(c(-0.25, -0.15, -0.03, 0, 0.03, 0.1, 0.15, 0.2, 0.3)),
    limits = c(-0.25, 0.4),
    oob = scales::squish
  ) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)", color = "Δ R²") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")


my_palette <- colorRampPalette(c("#a50026",
                                 "#d73027",
                                 "#f46d43",
                                 "#fdae61",
                                 "#CAF0F8",
                                 "#90E0EF",
                                 "#00A5D0",
                                 "#0077B6",
                                 "#023E8A",
                                 "#032174",
                                 "#030455"))

ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) +
  scale_color_gradientn(
    colors = my_palette(1000),
    limits = c(-0.25, 0.4),
    oob = scales::squish,
    name = "Δ R²"
  ) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")


ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) +
  scale_color_gradientn(
    colors = my_palette(1000),
    limits = range(predictions$r2_diff),
    name = "Δ R²"
  ) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")


ggplot(predictions, aes(x = 10^log_ppt, y = predicted_mass, colour = r2_diff)) +
  geom_line(aes(group = site_code)) +
  scale_fill_distiller(palette = "Spectral", 
                       aesthetics = "colour",
                       direction = 1,
                       limits = c(-0.1, 0.2),
                       oob = scales::squish) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            color = "black", linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m²)") +
  facet_wrap(~ trt) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")


## Looking at site-level trade-offs in magnitude of response to limiting factors

par_lrr_mass_plot <- ggplot(data = results_with_averages, 
                            aes(x = avg_lrr_mass, y = avg_proportion_par, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("Proportion PAR") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c"))
par_lrr_mass_plot

slope_lrr_mass_plot2 <- ggplot(data = results_with_averages, 
                               aes(x = avg_lrr_mass, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Log Response Ratio of Mass") + ylab("Slope of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c"))

slope_par_plot2 <- ggplot(data = results_with_averages, 
                          aes(x = avg_proportion_par, y = slope, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm, se = FALSE) +
  xlab("Proportion PAR") + ylab("Slope of ppt vs. mass") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c"))

par_lrr_mass_model <- lm(avg_proportion_par ~ trt * avg_lrr_mass, data = results_with_averages, na.action = "na.fail")
summary(par_lrr_mass_model) # NS

slope_lrr_mass_model <- lm(slope ~ trt * avg_lrr_mass, data = results_with_averages, na.action = "na.fail")
summary(slope_lrr_mass_model) # NS

slope_par_model <- lm(slope ~ trt * avg_proportion_par, data = results_with_averages, na.action = "na.fail")
summary(slope_par_model) # NS

colimitation_figure <- ggarrange(slope_lrr_mass_plot2, slope_par_plot2, par_lrr_mass_plot,
                                 ncol = 1, common.legend = TRUE, legend = "bottom", align = 'hv')
colimitation_figure

ai_lrr_model <- lm(avg_lrr_mass ~ avg_ai, data = results_with_averages)
summary(ai_lrr_model)

ai_lrr_plot <- ggplot(data = results_with_averages, aes(x = avg_ai, y = avg_lrr_mass)) +
  geom_point(color = "darkgrey") + geom_smooth(method = lm, se = FALSE, color = "black") +
  labs(x = "Aridity Index",  y = "LRR mass") +
  theme_bw() 
ai_lrr_plot



## Code loop to model mass vs. each predictor individually and summarize results in a table

library(lme4)
library(purrr)
library(glue)
library(broom.mixed)
library(MuMIn)
library(dplyr)
library(stringr)

predictors <- c(
  "log_ppt",
  "proportion_par",
  "AI",
  "rich",
  "prev_ppt",
  "lrr_mass",
  "avg_c4_proportion",
  "avg_annual_proportion"
)

base_formula <- "log_mass ~ trt * {predictor} + (1 | site_code/year_trt) + (1 | site_code/block)"

models <- map(
  predictors,
  ~ lmer(as.formula(glue(base_formula, predictor = .x)), data = mass_ppt_edited)
)

names(models) <- paste0(predictors, "_mass_model")

model_summaries <- map_dfr(
  names(models),
  function(name) {
    model <- models[[name]]
    tidy_df <- broom.mixed::tidy(model, effects = "fixed") %>% mutate(model = name)
    r2 <- suppressWarnings(MuMIn::r.squaredGLMM(model))
    tidy_df %>%
      mutate(Marginal_R2 = r2[1], Conditional_R2 = r2[2])
  }
)

summary_table <- model_summaries %>%
  dplyr::select(model, term, estimate, std.error, statistic, p.value, Marginal_R2, Conditional_R2)

summary_table



## Analyzing only data from the driest year and from below 10th ppt percentile at each site

driest_year_mass_ppt <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  slice_min(ppt) %>%
  ungroup()

main_model_driest <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year_trt), 
                          data = driest_year_mass_ppt)

summary(main_model_driest)

driest_year_plot <- ggplot(data = driest_year_mass_ppt, 
                           aes(x = ppt, y = live_mass)) +
  geom_point(aes(color = trt, fill = trt, shape = trt), alpha = 0.7) + 
  geom_smooth(method = "lm", color = "#6F6F6F", alpha = 0.3) +
  labs(x = "Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

driest_year_plot


main_model_p10 <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year_trt), 
                       data = subset(mass_ppt_edited, ppt<p10_ppt))

summary(main_model_p10)

p10_plot <- ggplot(data = subset(mass_ppt_edited, ppt<p10_ppt), 
                   aes(x = ppt, y = live_mass)) +
  geom_point(aes(color = trt, fill = trt, shape = trt), alpha = 0.7) + 
  geom_smooth(method = "lm", color = "#6F6F6F", alpha = 0.3) +
  labs(x = "Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

p10_plot


## Analyzing only data from the wettest year and from above the 90th ppt percentile at each site

wettest_year_mass_ppt <- mass_ppt_edited %>%
  group_by(site_code, trt) %>%
  slice_max(ppt) %>%
  ungroup()

main_model_wettest <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year_trt), 
                           data = wettest_year_mass_ppt)

summary(main_model_wettest)

wettest_year_plot <- ggplot(data = wettest_year_mass_ppt, 
                            aes(x = ppt, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", alpha = 0.2) +
  labs(x = "Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

wettest_year_plot


main_model_p90 <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year_trt), 
                       data = subset(mass_ppt_edited, ppt>p90_ppt))

summary(main_model_p90)

wettest_year_plot_p90 <- ggplot(data = subset(mass_ppt_edited, ppt>p90_ppt), 
                                aes(x = ppt, y = live_mass, color = trt, fill = trt, shape = trt)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm", alpha = 0.2) +
  labs(x = "Precipitation (mm)", y = "Biomass (g/m²)", 
       color = "Treatment", shape = "Treatment", fill = "Treatment") +
  theme_bw() +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

wettest_year_plot_p90


## old code for back-transforming and graphing ppt-biomass faceted by site

fit_model_and_predict <- function(data) {
  model <- lm(log_mass ~ log_ppt, data = data)
  new_data <- data.frame(log_ppt = seq(min(data$log_ppt, na.rm = TRUE),
                                       max(data$log_ppt, na.rm = TRUE),
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

ggplot(mass_ppt, aes(x = ppt, y = live_mass, color = trt)) +
  geom_point(alpha = 0.7) +
  geom_line(data = predictions, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Annual Growing Season Precipitation (mm)", y = "Biomass (g/m²)", color = "Treatment") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")



