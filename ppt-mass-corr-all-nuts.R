library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

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
unique(mass_ppt$trt)

mass_ppt_nuts <- filter(mass_ppt, trt %in% c("K", "NP", "Control", "NK", "PK", "NPK", "P", "N")) 
unique(mass_ppt_nuts$trt)

desired_order <- c("Control", "K", "P", "N", "PK", "NK", "NP", "NPK")
mass_ppt_nuts$trt <- factor(mass_ppt_nuts$trt, levels = desired_order)

# Create the live_mass column by summing vascular_live_mass and nonvascular_live_mass
mass_ppt_nuts <- mass_ppt_nuts %>%
  mutate(live_mass = ifelse(is.na(vascular_live_mass) & is.na(nonvascular_live_mass),
                            NA,
                            rowSums(cbind(vascular_live_mass, nonvascular_live_mass), na.rm = TRUE)))


mass_ppt_nuts <- mass_ppt_nuts %>%
  mutate(log_live_mass = log10(live_mass),
         log_mswep_ppt = log10(mswep_ppt))

site_year_counts <- mass_ppt_nuts %>%
  group_by(site_code, trt) %>%
  filter(!is.na(live_mass)) %>% 
  summarise(year_count = n_distinct(year), .groups = 'drop')

sites_with_8_years <- site_year_counts %>%
  filter(year_count >= 8) %>%
  group_by(site_code) %>% 
  filter(n_distinct(trt) == 8) %>% 
  select(site_code)

unique(sites_with_8_years$site_code)

mass_ppt_nuts <- mass_ppt_nuts %>%
  filter(site_code %in% sites_with_8_years$site_code)

unique(mass_ppt_nuts$site_code)

# Filter to keep only sites where ppt range spans at least +- 1 standard deviation
mass_ppt_nuts <- mass_ppt_nuts %>%
  group_by(site_code) %>%
  filter(max(mswep_ppt_sd) > 1 & min(mswep_ppt_sd) < -1) %>%
  ungroup()

unique(mass_ppt_nuts$site_code)

mass_ppt_nuts <- filter(mass_ppt_nuts, !(site_code == 'koffler.ca' & year == 2022 & plot == 2))
mass_ppt_nuts <- filter(mass_ppt_nuts, !(site_code == 'smith.us' & year == 2016 & plot == 1))
mass_ppt_nuts <- filter(mass_ppt_nuts, !(site_code == 'marc.ar' & year == 2017 & plot == 24))
mass_ppt_nuts <- filter(mass_ppt_nuts, !(site_code == 'bldr.us' & year == 2012 & plot == 5))
mass_ppt_nuts <- filter(mass_ppt_nuts, !(site_code == 'bnch.us' & year == 2021 & plot == 2))

ggplot(data = subset(mass_ppt_nuts, !is.na(live_mass)), aes(x= mswep_ppt, y= live_mass, color = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(data = subset(mass_ppt_nuts, !is.na(live_mass)), aes(x= mswep_ppt, y= live_mass, color = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  theme_bw()

# Calculate R² for each treatment at each site
calculate_r2 <- function(df) {
  model <- lm(live_mass ~ mswep_ppt, data = df)
  summary(model)$r.squared
}

# Apply R² calculation to each combination of trt and site_code
r2_results <- mass_ppt_nuts %>%
  group_by(site_code, trt) %>%
  summarise(r2 = calculate_r2(cur_data())) %>%
  ungroup()

View(r2_results)

# Create a function to compare non-control trts against Control
compare_to_control <- function(df) {
  control_r2 <- df %>% filter(trt == "Control") %>% pull(r2)
  df %>%
    filter(trt != "Control") %>%
    mutate(diff_r2 = r2 - control_r2)
}

# Apply the comparison for each site
comparison_results <- r2_results %>%
  group_by(site_code) %>%
  do(compare_to_control(.)) %>%
  ungroup()

View(comparison_results)

# Calculate average and standard deviation of R2 for each treatment
average_r2_by_trt <- r2_results %>%
  group_by(trt) %>%
  summarise(avg_r2 = mean(r2, na.rm = TRUE),
            sd_r2 = sd(r2))

View(average_r2_by_trt)

# Fit a linear mixed-effects model
lme_model <- lmer(r2 ~ trt + (1 | site_code), data = r2_results)

# Perform pairwise comparisons of treatments
emmeans_trt <- emmeans(lme_model, pairwise ~ trt)
emmeans_trt$contrasts


