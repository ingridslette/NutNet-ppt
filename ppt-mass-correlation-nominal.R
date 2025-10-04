
### Doing everything again using only non-extreme ppt values
### Picking up after line 106 of main script (right before running main model)

# Filter to keep only ppt values within non-extreme range
mass_ppt_nominal <- mass_ppt %>%
  filter(ppt > p10_ppt, ppt < p90_ppt)

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



library(ggplot2)

ggplot(mass_ppt_nominal, aes(x = ppt, y = live_mass, color = site_code)) +
  # points
  geom_point(alpha = 0.7, size = 2) +
  
  # site-specific regression lines (colored)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, alpha = 0.7) +
  
  # overall treatment-level regression line (black)
  geom_smooth(aes(group = trt), method = "lm", se = FALSE,
              color = "black", linewidth = 1.5) +
  
  facet_wrap(~ trt, scales = "free") +
  theme_bw() +
  labs(
    x = "Precipitation (ppt)",
    y = "Live Mass",
    color = "Site",
    title = "Precipitation vs. Live Mass by Treatment and Site"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "gray70"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )




library(ggplot2)
library(lme4)
library(dplyr)
library(ggeffects)

# 1. Overall fixed-effects predictions (log scale → back-transform)
overall_pred <- ggpredict(main_model, terms = c("log_ppt [all]", "trt"), type = "fixed") %>%
  as.data.frame() %>%
  mutate(
    ppt = 10^x,                      # back-transform x-axis
    live_mass = 10^predicted,        # back-transform predictions
    conf.low_bt = 10^conf.low,       # back-transform lower CI
    conf.high_bt = 10^conf.high      # back-transform upper CI
  )

# 2. Site-specific predictions (log scale → back-transform)
site_pred <- ggpredict(main_model, terms = c("log_ppt [all]", "trt", "site_code"), type = "random") %>%
  as.data.frame() %>%
  mutate(
    ppt = 10^x,
    live_mass = 10^predicted
  )

# 3. Plot
ggplot(mass_ppt_nominal, aes(x = ppt, y = live_mass, color = site_code)) +
  geom_point(alpha = 0.6, size = 2) +
  
  geom_line(data = site_pred,
            aes(x = ppt, y = live_mass, color = group),
            linewidth = 0.7, alpha = 0.7) +
  
  geom_ribbon(
    data = overall_pred,
    aes(x = ppt, ymin = conf.low_bt, ymax = conf.high_bt, fill = group),
    inherit.aes = FALSE, alpha = 0.2, color = NA
  ) +
  
  geom_line(
    data = overall_pred,
    aes(x = ppt, y = live_mass),
    color = "black", linewidth = 1.5
  ) +
  
  facet_wrap(~ group, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Precipitation (ppt)",
    y = "Live Mass",
    color = "Site",
    fill = "Treatment") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )


