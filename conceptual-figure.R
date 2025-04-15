library(ggplot2)

# Create a data frame for plotting
x <- seq(0, 10, length.out = 25)
y1 <- 2 + 0.5 * x  
y1.2 <- 2 + 0.5 * x
y1.3 <- 2 + 0.5 * x
y1.4 <- 2 + 0.5 * x
y2 <- 3.5 + 1 * x  # Steeper slope
y3 <- 3.8 + 0.5 * x  
y4 <- 3.5 + 1 * x  # Steeper slope

# Standard errors
se1 <- 1.2 * x / 10
se1.2 <- 1.2 * x / 10
se1.3 <- 1.2 * x / 10
se1.4 <- 1.2 * x / 10
se2 <- 1.2 * x / 10   
se3 <- 0.6 * x / 10   # Smaller standard error
se4 <- 0.6 * x / 10   # Smaller standard error

data_main_line <- data.frame(
  x = rep(x, 4),
  y = c(y1, y1.2, y1.3, y1.4),
  ymin = c(y1 - se1, y1.2 - se1.2, y1.3 - se1.3, y1.4 - se1.4),
  ymax = c(y1 + se1, y1.2 + se1.2, y1.3 + se1.3, y1.4 + se1.4),
  panel = factor(rep(1:4, each = length(x))),
  Treatment = "Control"
)

data_unique_line <- data.frame(
  x = rep(x, 3),
  y = c(y2, y3, y4),
  ymin = c(y2 - se2, y3 - se3, y4 - se4),
  ymax = c(y2 + se2, y3 + se3, y4 + se4),
  panel = factor(rep(2:4, each = length(x))),
  Treatment = "Fertilized"
)

final_data <- rbind(data_main_line, data_unique_line)

fig1 <- ggplot(final_data, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Treatment), alpha = 0.3, show.legend = FALSE) +
  geom_line(aes(color = Treatment), linewidth = 1.2) +
  facet_wrap(.~panel, nrow = 2, ncol = 2) +
  theme_bw(14) +
  scale_color_manual(values = c("#4267ac", "#ff924c")) +
  scale_fill_manual(values = c("#4267ac", "#ff924c")) +
  labs(x = "Precipitation", y = "Biomass") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank()
  )

fig1

fig1_labels <- data.frame(panel = c(1,2,3,4), label = c("a", "b", "c", "d"))

fig1 + geom_text(x = 0.5, y = 14, aes(label = label), data = fig1_labels)


