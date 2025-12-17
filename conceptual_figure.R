library(ggplot2)
library(dplyr)
library(tidyr)

x <- seq(0, 10, length.out = 25)

panel_specs <- tibble(
  panel = factor(1:4),
  
  # Control (constant across panels)
  intercept_control = 2.0,
  slope_control     = 0.5,
  
  # Fertilized (intercept fixed, slope varies)
  slope_fert = c(0.5, 1.0, 0.5, 1.0),
  
  # Fertilized-only SE scaling
  se_fert = c(1, 1, 0.5, 0.5)
)

# Constants
intercept_fert <- 4.0
se_control     <- 1


plot_data <- panel_specs %>%
  crossing(
    x = x,
    Treatment = c("Control", "Fertilized")
  ) %>%
  mutate(
    intercept = if_else(
      Treatment == "Control",
      intercept_control,
      intercept_fert
    ),
    
    slope = if_else(
      Treatment == "Control",
      slope_control,
      slope_fert
    ),
    
    y = intercept + slope * x,
    
    se = if_else(
      Treatment == "Control",
      se_control * x / 10,
      se_fert * x / 10
    ),
    
    ymin = y - se,
    ymax = y + se
  )


fig1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Treatment), alpha = 0.3) +
  geom_line(aes(color = Treatment), linewidth = 1.2) +
  facet_wrap(~ panel, nrow = 1, ncol = 4) +
  theme_bw(14) +
  scale_color_manual(values = c("#0092E0", "#ff924c")) +
  scale_fill_manual(values = c("#0092E0", "#ff924c")) +
  labs(x = "Precipitation", y = "Biomass") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

fig1

fig1 +
  geom_text(
    data = tibble(panel = factor(1:4), label = letters[1:4]),
    aes(x = 0.5, y = 14, label = label),
    inherit.aes = FALSE
  )
