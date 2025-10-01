
# Doing everything again using only non-extreme ppt values
# Picking up after line 106 of main script (right before running main model)

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

# Join results
final_results_site <- site_results %>%
  left_join(site_r2, by = "site_code")

# final_results now contains:
# site_code | term | estimate | std.error | statistic | df | p.value | R2_marginal | R2_conditional



# Fit lm for each site × trt
site_trt_models <- mass_ppt_nominal %>%
  group_by(site_code, trt) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(unique(.x$site_code), unique(.x$trt), sep = "_"))) %>%
  map(~ lm(log_mass ~ log_ppt, data = .x))

# Extract coefficient-level estimates
site_trt_estimates <- map2_df(
  site_trt_models,
  names(site_trt_models),
  ~ tidy(.x) %>%
    rename(
      p.value_coef = p.value,
      t.value = statistic
    ) %>%
    mutate(site_trt = .y)
)

# Extract model-level R² and p-value
site_trt_r2 <- map2_df(
  site_trt_models,
  names(site_trt_models),
  ~ glance(.x) %>%
    dplyr::select(r.squared, adj.r.squared, p.value, df = df.residual) %>%
    rename(
      p.value_model = p.value
    ) %>%
    mutate(site_trt = .y)
)

# Split site_trt into separate columns
site_trt_estimates <- site_trt_estimates %>%
  separate(site_trt, into = c("site_code", "trt"), sep = "_")

site_trt_r2 <- site_trt_r2 %>%
  separate(site_trt, into = c("site_code", "trt"), sep = "_")

# Join coefficient and model-level results
final_results_site_trt <- site_trt_estimates %>%
  left_join(site_trt_r2, by = c("site_code", "trt"))

# final_results now contains:
# site_code | trt | term | estimate | std.error | t.value | p.value_coef | r.squared | adj.r.squared | p.value_model | df


# Keep model-level metrics only for slope term
final_results_clean <- final_results %>%
  mutate(
    r.squared = ifelse(term != "log_ppt", NA, r.squared),
    adj.r.squared = ifelse(term != "log_ppt", NA, adj.r.squared),
    p.value_model = ifelse(term != "log_ppt", NA, p.value_model)
  )

# Optional: reorder columns for clarity
final_results_clean <- final_results_clean %>%
  dplyr::select(site_code, trt, term, estimate, std.error, t.value, p.value_coef, df,
         r.squared, adj.r.squared, p.value_model)


