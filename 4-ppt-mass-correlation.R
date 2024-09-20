library(tidyverse)
library(lme4)
library(lmerTest)
library(boot)
library(MuMIn)
library(performance)

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/mswep_ppt_annual_gs_only.csv")

mass <- read.csv("/Users/ingridslette/Desktop/NutNet/full-biomass-2024-09-17.csv")

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

sites_with_8_years <- site_year_counts %>%
  filter(year_count >= 8) %>%
  select(site_code)

mass3 <- mass2 %>%
  filter(site_code %in% sites_with_8_years$site_code)

total_mass <- aggregate(
  mass ~ year + year_trt + trt + site_name + site_code + block + plot + subplot, 
  data = mass3, sum)

#total_mass_plots_avg <- aggregate(mass ~ year + year_trt + trt + site_name + site_code, data = mass3, mean)

mass_ppt <- inner_join(total_mass, mswep, by=c("site_code", "year"))

unique(mass_ppt$site_code)

mass_ppt <- mass_ppt %>%
  mutate(log_mass = log10(mass),
         log_mswep_ppt = log10(mswep_ppt))

desired_order <- c("Control", "K", "P", "N", "PK", "NK", "NP", "NPK")
mass_ppt$trt <- factor(mass_ppt$trt, levels = desired_order)

mass_ppt_c_npk <- filter(mass_ppt, trt %in% c("Control", "NPK")) 

str(mass_ppt_c_npk)

# remove high biomass outlier
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(!(year == 2021 & site_code == "ukul.za" & plot == 23))

mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(!(year == 2015 & site_code == "trel.us" & plot == 24))

# remove 2013 comp.pt - mass data is incorrect, I'm looking into it...
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(!(year == 2013 & site_code == "comp.pt"))

# remove 2022 ukul.za - data looks incorrect, I'm looking into it...
mass_ppt_c_npk <- mass_ppt_c_npk %>%
  filter(!(year == 2022 & site_code == "ukul.za"))

ggplot(mass_ppt_c_npk, aes(x= mswep_ppt, y= mass, color = trt, shape = trt, 
                           label = site_code
                           )) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  geom_text(aes(label=ifelse(mass>2500, as.character(site_code), '')), hjust=-0.1, vjust=0.1) +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x= mswep_ppt, y= mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("MSWEP Growing Season Precipitation (mm)") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x=year_trt, y=mass, color = trt, shape = trt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("Treatment Year") + ylab("Total live mass") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(mass_ppt_c_npk, aes(x= mswep_ppt, y= mass)) +
  geom_smooth(aes(group = site_code, color = site_code), method = "lm", se = FALSE) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  facet_wrap(~ trt, nrow = 2) +
  theme_bw() +
  labs(x = "Growing Season Precipitation",
       y = "Total live mass",
       color = "Site Code") +
  theme(legend.position = "right")

c_npk_x_model <- lmer(mass ~ mswep_ppt * trt + (1 | site_code) + (1 | year_trt), data = mass_ppt_c_npk)
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
  # Subset data for the site
  site_data_control <- subset(mass_ppt_c_npk, site_code == site & trt == "Control")
  site_data_npk <- subset(mass_ppt_c_npk, site_code == site & trt == "NPK")
  # Only proceed if both Control and NPK have sufficient data for the site
  if (nrow(site_data_control) > 1 & nrow(site_data_npk) > 1) {
    # Fit the linear models for Control and NPK
    control_model <- lm(mass ~ mswep_ppt, data = site_data_control)
    npk_model <- lm(mass ~ mswep_ppt, data = site_data_npk)
    # Extract the R2 values
    control_r2 <- summary(control_model)$r.squared
    npk_r2 <- summary(npk_model)$r.squared
    # Calculate the difference in R2
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

# Perform a t-test on the r2_difference values
t_test_r2_diff <- t.test(results$r2_difference)
print(t_test_r2_diff)

# Perform a paired t-test on the R^2 values for Control and NPK
paired_t_test_result <- t.test(results$control_r2, results$npk_r2, paired = TRUE)
print(paired_t_test_result)


### Approach 2: fit separate models for control and NPK data, calculate and compare z scores

# Fit linear mixed-effects models for each treatment level
model_control <- lmer(mass ~ mswep_ppt + (1 | site_code) + (1 | year_trt), data = subset(mass_ppt_c_npk, trt == "Control"))
model_npk <- lmer(mass ~ mswep_ppt + (1 | site_code) + (1 | year_trt), data = subset(mass_ppt_c_npk, trt == "NPK"))

summary(model_control)
summary(model_npk)

# Compare AIC of the two models
AIC(model_control, model_npk)

# Use the performance package to calculate R2
r2_control <- performance::r2(model_control)
r2_npk <- performance::r2(model_npk)

# Extract R2 values
conditional_r2_control <- r2_control$R2_conditional
conditional_r2_npk <- r2_npk$R2_conditional

# Compare R2 values using Fisher's Z transformation
z_control <- 0.5 * log((1 + sqrt(conditional_r2_control)) / (1 - sqrt(conditional_r2_control)))
z_npk <- 0.5 * log((1 + sqrt(conditional_r2_npk)) / (1 - sqrt(conditional_r2_npk)))

# Calculate the standard error
n_control <- length(unique(subset(mass_ppt_c_npk, trt == "Control")$site_code))
n_npk <- length(unique(subset(mass_ppt_c_npk, trt == "NPK")$site_code))
se_diff <- sqrt((1 / (n_control - 3)) + (1 / (n_npk - 3)))

# Calculate the Z-score for the difference
z_diff <- (z_control - z_npk) / se_diff

# Calculate the p-value
p_value <- 2 * (1 - pnorm(abs(z_diff)))

# Print the results
cat("Z-score for the difference:", z_diff, "\n")
cat("P-value for the difference in conditional R-squared values:", p_value, "\n")


### Approach 3: Bootstrapping to test for difference in R2 between Control and NKP models

# Function to calculate the difference in marginal R²
r2_diff <- function(data, indices) {
  data_resampled <- data[indices, ]
  model_control <- lmer(mass ~ mswep_ppt + (1 | site_code) + (1 | year_trt), 
                        data = data_resampled[data_resampled$trt == "Control", ])
  model_npk <- lmer(mass ~ mswep_ppt + (1 | site_code) + (1 | year_trt), 
                    data = data_resampled[data_resampled$trt == "NPK", ])
  r2_control <- r.squaredGLMM(model_control)[2]
  r2_npk <- r.squaredGLMM(model_npk)[2]          
  return(r2_control - r2_npk)
}

# Set seed for reproducibility
set.seed(123)

# Perform bootstrapping with 1000 resamples
boot_r2 <- boot(data = mass_ppt_c_npk, statistic = r2_diff, R = 1000)

print(boot_r2)

# Get 95% confidence intervals for the R² difference
boot_ci <- boot.ci(boot_r2, type = "perc")
print(boot_ci)
