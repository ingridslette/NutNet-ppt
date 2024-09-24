library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(MuMIn)
library(performance)

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/mswep_ppt_annual_gs_only.csv")

mswep <- filter(mswep, year >= 1983)

mswep <- mswep %>%
  group_by(site_code) %>%
  mutate(avg_ppt_site = mean(mswep_ppt, na.rm = TRUE),
         sd_ppt_site = sd(mswep_ppt, na.rm = TRUE),
         mswep_ppt_per = (mswep_ppt / avg_ppt_site) * 100, 
         mswep_ppt_sd = (mswep_ppt - avg_ppt_site) / sd_ppt_site) %>%
  ungroup()

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/comb-by-plot-2024-09-24.csv",
                 na.strings = c("NULL","NA"))

unique(mass$site_code)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, year_trt > 0)
unique(mass1$year_trt)

mass_ppt <- inner_join(mass1, mswep, by=c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt_c_npk <- filter(mass_ppt, trt %in% c("Control", "NPK")) 

unique(mass_ppt_c_npk$site_code)
unique(mass_ppt_c_npk$trt)

site_year_counts <- mass_ppt_c_npk %>%
  group_by(site_code, trt) %>%
  filter(!is.na(vascular_live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_8_years <- site_year_counts %>%
  filter(year_count >= 8) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 2) %>% 
  select(site_code)

mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(site_code %in% sites_with_8_years$site_code)

unique(mass_ppt_c_npk$site_code)

# Filter to keep only sites with a certain ppt range
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  group_by(site_code) %>%
  filter(max(mswep_ppt_sd) > 1 & min(mswep_ppt_sd) < -1) %>%
  ungroup()

unique(mass_ppt_c_npk$site_code)

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt_per, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt_per, y= vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x=year_trt, y=vascular_live_mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Treatment Year") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_c_npk, !is.na(vascular_live_mass)), aes(x= mswep_ppt_per, y= vascular_live_mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation",
       y = "Total live mass",
       color = "Site Code") +
  theme(legend.position = "right")

c_npk_x_model <- lmer(vascular_live_mass ~ mswep_ppt_per * trt + (1 | site_code) + (1 | year_trt), data = mass_ppt_c_npk)
summary(c_npk_x_model)


### Approach 1: calculate and compare difference in R2 between Control and NKP at each site

# Get unique site codes
site_codes <- unique(mass_ppt_c_npk$site_code)

# Initialize a dataframe to store results
results <- data.frame(site_code = character(), 
                      control_r2 = numeric(), 
                      npk_r2 = numeric(), 
                      r2_difference = numeric(), 
                      stringsAsFactors = FALSE)

for (site in site_codes) {
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  control_model <- lm(vascular_live_mass ~ mswep_ppt_per, data = site_data_control)
  npk_model <- lm(vascular_live_mass ~ mswep_ppt_per, data = site_data_npk)
  control_r2 <- summary(control_model)$r.squared
  npk_r2 <- summary(npk_model)$r.squared
  r2_difference <- control_r2 - npk_r2
  results <- rbind(results, data.frame(
    site_code = site,
    control_r2 = control_r2,
    npk_r2 = npk_r2,
    r2_difference = r2_difference
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


### Approach 2: fit separate models for control and NPK data, calculate and compare z scores

model_control <- lmer(vascular_live_mass ~ mswep_ppt_per + (1 | site_code) + (1 | year_trt), 
                      data = subset(mass_ppt_c_npk, trt == "Control"))
model_npk <- lmer(vascular_live_mass ~ mswep_ppt_per + (1 | site_code) + (1 | year_trt), 
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


### Approach 3: Bootstrapping to test for difference in R2 between Control and NKP models

r2_diff <- function(data, indices) {
  data_resampled <- data[indices, ]
  model_control <- lmer(vascular_live_mass ~ mswep_ppt_per + (1 | site_code), 
                        data = data_resampled[data_resampled$trt == "Control", ])
  model_npk <- lmer(vascular_live_mass ~ mswep_ppt_per + (1 | site_code), 
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


### Investigating the effect of light limitation on the change in R2
## Using approach 2 from above: fit separate models for different trts, calculate and compare z scores

model_control <- lmer(vascular_live_mass ~ mswep_ppt * proportion_par + (1 | site_code) + (1 | year_trt), 
                      data = subset(mass_ppt_c_npk, trt == "Control"))
model_npk <- lmer(vascular_live_mass ~ mswep_ppt * proportion_par + (1 | site_code) + (1 | year_trt), 
                  data = subset(mass_ppt_c_npk, trt == "NPK"))

summary(model_control)
summary(model_npk)

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

