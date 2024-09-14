library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(broom.mixed)
library(boot)
library(MuMIn)
library(purrr)
library(tidyr)

mass <- read.csv("/Users/ingridslette/Desktop/full-biomass_2024-05-31.csv")

unique(mass$site_code)
unique(mass$live)
unique(mass$trt)
unique(mass$year_trt)

mass1 <- filter(mass, live == 1)
mass1 <- filter(mass1, trt %in% c("Control", "K", "P", "N", "PK", "NK", "NP", "NPK"))
mass1 <- filter(mass1, year_trt > 0)

unique(mass1$live)
unique(mass1$trt)
unique(mass1$year_trt)

unique(mass1$category)
mass2 <- filter(mass1, category != 'WOODY')
unique(mass2$category)

site_year_counts <- mass2 %>%
  group_by(site_code) %>%
  summarise(year_count = n_distinct(year))

sites_with_10_years <- site_year_counts %>%
  filter(year_count >= 10) %>%
  select(site_code)

mass3 <- mass2 %>%
  filter(site_code %in% sites_with_10_years$site_code)

total_mass <- aggregate(
  mass ~ year + year_trt + trt + site_name + site_code + block + plot + subplot, 
  data = mass3, sum)

mass_ppt <- inner_join(total_mass, mswep_cru, by=c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(mass),
         log_mswep_ppt = log10(mswep_ppt),
         log_cru_ppt = log10(cru_ppt))

desired_order <- c("Control", "K", "P", "N", "PK", "NK", "NP", "NPK")
mass_ppt$trt <- factor(mass_ppt$trt, levels = desired_order)

mass_ppt_c_npk <- filter(mass_ppt, trt %in% c("Control", "NPK")) 

ggplot(mass_ppt_c_npk, aes(x=mswep_ppt, y=mass, color = trt, shape = trt, label = site_code)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  geom_text(aes(label=ifelse(mass>2500, as.character(site_code), '')), hjust=-0.1, vjust=0.1) +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x=mswep_ppt, y=mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x = mswep_ppt, y = mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation",
       y = "Total live mass",
       color = "Site Code") +
  theme(legend.position = "right")

c_npk_x_model <- lm(mass ~ mswep_ppt * trt, data = mass_ppt_c_npk)
summary(c_npk_x_model)

#mass_ppt_c_npk_filtered <- mass_ppt_c_npk %>% filter(log_mass != -Inf)

#c_npk_x_model2 <- lm(log_mass ~ log_mswep_ppt * trt, data = mass_ppt_c_npk_filtered)
#summary(c_npk_x_model2)

control_data <- subset(mass_ppt_c_npk, trt == "Control")
npk_data <- subset(mass_ppt_c_npk, trt == "NPK")

# Define a function to calculate R^2 for bootstrapping
calc_r2 <- function(data, indices) {
  # Resample the data
  d <- data[indices, ]
  # Fit the linear model
  model <- lm(mass ~ mswep_ppt, data = d)
  # Return the R^2 value
  return(summary(model)$r.squared)
}
# Set number of bootstrap iterations
n_boot <- 1000

# Perform bootstrap for "Control" group
control_boot <- boot(data = control_data, statistic = calc_r2, R = n_boot)
control_r2_bootstrap <- control_boot$t
# Perform bootstrap for "NPK" group
npk_boot <- boot(data = npk_data, statistic = calc_r2, R = n_boot)
npk_r2_bootstrap <- npk_boot$t

# Print the means and standard deviations of the bootstrapped R^2 values
cat("Mean R^2 for Control (bootstrapped):", mean(control_r2_bootstrap), "\n")
cat("Mean R^2 for NPK (bootstrapped):", mean(npk_r2_bootstrap), "\n")
cat("Standard deviation of R^2 for Control (bootstrapped):", sd(control_r2_bootstrap), "\n")
cat("Standard deviation of R^2 for NPK (bootstrapped):", sd(npk_r2_bootstrap), "\n")

# Perform t-test to compare the two distributions of R^2 values
r2_t_test <- t.test(control_r2_bootstrap, npk_r2_bootstrap)
print(r2_t_test)

# Get unique site codes
site_codes <- unique(mass_ppt_c_npk$site_code)
# Initialize a dataframe to store results
results <- data.frame(site_code = character(), 
                      control_r2 = numeric(), 
                      npk_r2 = numeric(), 
                      r2_difference = numeric(), 
                      stringsAsFactors = FALSE)

# Loop over each site_code
for (site in site_codes) {
  # Subset data for the site
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  # Only proceed if both Control and NPK have sufficient data for the site
  if (nrow(site_data_control) > 1 & nrow(site_data_npk) > 1) {
    # Fit the linear models for Control and NPK
    control_model <- lm(mass ~ mswep_ppt, data = site_data_control)
    npk_model <- lm(mass ~ mswep_ppt, data = site_data_npk)
    # Extract the R^2 values
    control_r2 <- summary(control_model)$r.squared
    npk_r2 <- summary(npk_model)$r.squared
    # Calculate the difference in R^2
    r2_difference <- control_r2 - npk_r2
    # Store the results in the dataframe
    results <- rbind(results, data.frame(
      site_code = site,
      control_r2 = control_r2,
      npk_r2 = npk_r2,
      r2_difference = r2_difference
    ))
  }
}

# Perform a one-sample t-test on the r2_difference values
t_test_r2_diff <- t.test(results$r2_difference)
print(t_test_r2_diff)

