
my_palette <- colorRampPalette(c("#a50026","#d73027","#f46d43","#fdae61",
                                 "#CAF0F8","#90E0EF","#00A5D0","#0077B6",
                                 "#023E8A","#032174","#030455"))



### Calculating and graphing variance in log_mass

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
  geom_bar(aes(y = prop_variance_conditional), stat = "identity", fill = "lightgrey") +
  labs(y = "Variance",
       x = "Treatment",
       fill = "Variance Type") +
  theme_bw(14)


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


