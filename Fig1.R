# Load libraries
library(ggplot2)
library(dplyr)

# Function to generate Gaussian phenology curve
phenology_curve <- function(DOY, peak, spread, amplitude = 1) {
  amplitude * exp(-((DOY - peak)^2) / (2 * spread^2))
}

# Function to compute overlap
compute_overlap <- function(resource, consumer) {
  pmin(resource, consumer)
}

# Set DOY range
DOY <- 1:365

# Base resource parameters
resource_peak <- 150
resource_spread <- 20
resource <- phenology_curve(DOY, resource_peak, resource_spread)

# Parameter grid for consumer
timing_shifts <- c(-20, 0, 20)         # earlier, aligned, later
consumer_spreads <- c(10, 20, 40)      # narrow, normal, wide

# Create all combinations
params <- expand.grid(timing_shift = timing_shifts,
                      consumer_spread = consumer_spreads)

# Generate long-format data for plotting
plot_data <- data.frame()

for (i in 1:nrow(params)) {
  shift <- params$timing_shift[i]
  spread <- params$consumer_spread[i]

  consumer <- phenology_curve(DOY, resource_peak + shift, spread)
  overlap <- compute_overlap(resource, consumer)

  df <- data.frame(
    DOY = rep(DOY, 3),
    value = c(resource, consumer, overlap),
    type = rep(c("Resource", "Consumer", "Overlap"), each = length(DOY)),
    timing_shift = shift,
    consumer_spread = spread
  )

  plot_data <- rbind(plot_data, df)
}

# Convert factors for facet labels with descriptive names
plot_data$timing_shift <- factor(plot_data$timing_shift,
                                 levels = c(-20, 0, 20),
                                 labels = c("earlier", "no change", "later"))

plot_data$consumer_spread <- factor(plot_data$consumer_spread,
                                    levels = c(10, 20, 40),
                                    labels = c("narrower", "no change", "broader"))

# Plot phenology curves with overlap shaded
ggplot(plot_data, aes(x = DOY, y = value, color = type, fill = type)) +
  geom_line(size = 1) +
  geom_area(data = subset(plot_data, type == "Overlap"), alpha = 0.4, fill = "grey") +
  scale_color_manual(values = c("Resource" = "navyblue", "Consumer" = "gold")) +
  facet_grid(consumer_spread ~ timing_shift) +
  labs(x = "Day of Year (DOY)", y = "Phenology Intensity") +
  theme_classic() +
  theme(legend.position = "bottom")
