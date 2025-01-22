# Load necessary library
library(ggplot2)

# Create a data frame for plotting
x <- seq(0, 10, length.out = 25)
y1 <- 2 + 0.5 * x  # First line
y1.2 <- 2 + 0.5 * x
y1.3 <- 2 + 0.5 * x
y1.4 <- 2 + 0.5 * x
y2 <- 3 + 1 * x  # Second line, steeper slope
y3 <- 3 + 0.5 * x  # Third line, same slope, higher intercept
y4 <- 3 + 1 * x  # Fourth line, steeper slope

# Standard errors
se1 <- 1 * x / 10   # Standard error for line 1
se1.2 <- 1 * x / 10
se1.3 <- 1 * x / 10
se1.4 <- 1 * x / 10
se2 <- 1 * x / 10   # Same as se1
se3 <- 0.5 * x / 10   # Smaller standard error for line 3
se4 <- 0.5 * x / 10   # Smaller standard error for line 4

data_main_line <- data.frame(
  x = rep(x, 4),
  y = c(y1, y1.2, y1.3, y1.4),
  ymin = c(y1 - se1, y1.2 - se1.2, y1.3 - se1.3, y1.4 - se1.4),
  ymax = c(y1 + se1, y1.2 + se1.2, y1.3 + se1.3, y1.4 + se1.4),
  panel = factor(rep(1:4, each = length(x))),
  line_type = "Main Line"
)

data_unique_line <- data.frame(
  x = rep(x, 3),
  y = c(y2, y3, y4),
  ymin = c(y2 - se2, y3 - se3, y4 - se4),
  ymax = c(y2 + se2, y3 + se3, y4 + se4),
  panel = factor(rep(2:4, each = length(x))),
  line_type = "Unique Line"
)

final_data <- rbind(data_main_line, data_unique_line)

fig1 <- ggplot(final_data, aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
  geom_ribbon(aes(fill = line_type), alpha = 0.2, show.legend = FALSE) +
  geom_line(aes(color = line_type), size = 1.2, show.legend = FALSE) +
  facet_wrap(~panel, nrow = 2, ncol = 2) +
  theme_bw() +
  scale_color_manual(values = c("#1982c4", "#af7ab3")) +
  scale_fill_manual(values = c("#1982c4", "#af7ab3")) +
  theme(
    strip.text = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

fig1
