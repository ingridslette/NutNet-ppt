
### Doing everything again using only non-extreme ppt values
### Picking up after line 106 of main script (right before running main model)

# Filter to keep only ppt values within non-extreme range
mass_ppt_nominal <- mass_ppt %>%
  filter(ppt > p05_ppt, ppt < p95_ppt)

View(mass_ppt_nominal)
unique(mass_ppt_nominal$site_code)

### Main model

main_model_nominal <- lmer(log_mass ~ log_ppt * trt + (1 | site_code / block) + (1 | year_trt), data = mass_ppt_nominal)
summary(main_model_nominal)

# Model assumptions check 
plot(main_model_nominal)
resid <- residuals(main_model_nominal)
hist(resid, breaks = 30, main = "Histogram of Residuals")
qqnorm(resid)
qqline(resid)
plot(fitted(main_model_nominal), resid, main = "Residuals vs Fitted")

### Models per site

# Split data by site_code
site_models <- mass_ppt_nominal %>%
  group_by(site_code) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$site_code)))

# Fit models for each site
fit_site_models <- map(site_models, ~ {
  lmer(log_mass ~ log_ppt * trt + (1 | year_trt), data = .x)
})

# Extract tidy fixed effect results for each model
site_results <- map2_df(
  fit_site_models,
  names(fit_site_models),
  ~ tidy(.x, effects = "fixed") %>%
    mutate(site_code = .y)
)

View(site_results)

# Extract R² values for each model
site_r2 <- map2_df(
  fit_site_models,
  names(fit_site_models),
  ~ {
    r2_vals <- suppressWarnings(MuMIn::r.squaredGLMM(.x))
    tibble(
      site_code = .y,
      R2_marginal = r2_vals[1, "R2m"],
      R2_conditional = r2_vals[1, "R2c"]
    )
  }
)

View(site_r2)

# Join results
site_results_r2 <- site_results %>%
  left_join(site_r2, by = "site_code")

View(site_results_r2)

# Reorder columns for clarity
site_results_r2 <- site_results_r2 %>%
  dplyr::select(site_code, term, estimate, std.error, statistic, df, p.value, R2_marginal, R2_conditional)

View(site_results_r2)

# Pivot site results wider
site_results_wide <- site_results %>%
  dplyr::select(site_code, term, estimate, std.error, statistic, df, p.value) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, std.error, statistic, df, p.value),
    names_sep = "_"
  )

View(site_results_wide)

# Join site R2 values
site_results_r2_wide <- site_results_wide %>%
  left_join(site_r2, by = "site_code") %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  arrange(site_code)

View(site_results_r2_wide)


### Models per site-trt

# Fit models for each site-trt
site_trt_models <- mass_ppt_nominal %>%
  group_by(site_code, trt) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(unique(.x$site_code), unique(.x$trt), sep = "_"))) %>%
  map(~ lm(log_mass ~ log_ppt, data = .x))

# Extract tidy fixed effect results for model coefficients
site_trt_results <- map2_df(
  site_trt_models,
  names(site_trt_models),
  ~ tidy(.x) %>%
    rename(
      p.value_coef = p.value,
      t.value = statistic
    ) %>%
    mutate(site_trt = .y) %>%
  separate(site_trt, into = c("site_code", "trt"), sep = "_")
)

View(site_trt_results)

# Extract model-level R² and p-value
site_trt_r2 <- map2_df(
  site_trt_models,
  names(site_trt_models),
  ~ glance(.x) %>%
    dplyr::select(r.squared, adj.r.squared, p.value, df = df.residual) %>%
    rename(
      p.value_model = p.value
    ) %>%
    mutate(site_trt = .y) %>%
    separate(site_trt, into = c("site_code", "trt"), sep = "_")
)

View(site_trt_r2)

# Join coefficient and R2 values
site_trt_results_r2 <- site_trt_results %>%
  left_join(site_trt_r2, by = c("site_code", "trt"))

View(site_trt_results_r2)

# Reorder columns for clarity
site_trt_results_r2 <- site_trt_results_r2 %>%
  dplyr::select(site_code, trt, term, estimate, std.error, t.value, p.value_coef, df,
         r.squared, adj.r.squared, p.value_model)

View(site_trt_results_r2)

# Pivot site-trt results wider
site_trt_results_wide <- site_trt_results_r2 %>%
  dplyr::select(site_code, trt, term, estimate, std.error, t.value, p.value_coef) %>%
  pivot_wider(
    names_from = term,
    values_from = c(estimate, std.error, t.value, p.value_coef),
    names_sep = "_"
  )

View(site_trt_results_wide)

# Join R2 values
site_trt_results_r2_wide <- site_trt_results_wide %>%
  left_join(site_trt_r2, by = c("site_code", "trt")) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>% 
  arrange(site_code, trt)

View(site_trt_results_r2_wide)


## Back-transforming data for graphing - allowing for non-linear curves on linear scale

# back transform from log-log scale
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

predictions <- mass_ppt_nominal %>%
  group_by(site_code, trt) %>%
  group_modify(~ fit_model_and_predict(.x)) %>%
  ungroup()

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

predictions_allsites <- mass_ppt_nominal %>%
  group_by(trt) %>%
  group_modify(~ fit_model_and_predict_allsites(.x)) %>%
  ungroup()

ggplot(data = mass_ppt_nominal, aes(x = ppt, y = live_mass, color = trt, shape = trt)) +
  geom_point() + 
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  xlab("Growing Season Precipitation (mm)") + ylab("Biomass (g/m2)") +
  labs(color = "Treatment", shape = "Treatment") +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  theme_bw(14)

ggplot(mass_ppt_nominal, aes(x = ppt, y = live_mass, color = trt)) +
  geom_point() +
  geom_line(data = predictions, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)", color = "Treatment") +
  facet_wrap(~ site_code, scales = "free") +
  theme_bw(14) +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  theme(legend.position = "bottom")

ggplot(mass_ppt_nominal, aes(x = ppt, y = live_mass, color = site_code)) +
  geom_line(data = predictions, aes(x = 10^log_ppt, y = predicted_mass), linewidth = 1) +
  geom_line(data = predictions_allsites, aes(x = 10^log_ppt, y = predicted_mass), 
            linewidth = 1, color = "black") +
  labs(x = "Growing Season Precipitation (mm)", y = "Biomass (g/m2)") +
  facet_wrap(~ trt) +
  theme_bw(14)





