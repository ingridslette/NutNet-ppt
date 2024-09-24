library(tidyverse)
library(lme4)

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

site_year_counts <- mass_ppt_nuts %>%
  group_by(site_code, trt) %>%
  filter(!is.na(vascular_live_mass)) %>% 
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

# Calculate R² for each treatment at each site
calculate_r2 <- function(df) {
  model <- lm(vascular_live_mass ~ mswep_ppt, data = df)
  summary(model)$r.squared
}

# Apply R² calculation to each combination of trt and site_code
r2_results <- mass_ppt_nuts %>%
  group_by(site_code, trt) %>%
  summarise(r2 = calculate_r2(cur_data())) %>%
  ungroup()

View(r2_results)





