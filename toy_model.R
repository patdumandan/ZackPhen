library(dplyr)
library(ggplot2)
library(purrr)
library(broom)

# Simulate sample data
set.seed(42)
years <- rep(2000:2009, each = 50)
doy <- rep(seq(100, 300, length.out = 50), 10)
intensity <- -0.01 * (doy - 200)^2 + 10 + rnorm(500, 0, 0.5)

df <- data.frame(year = years, doy = doy, intensity = intensity)

# Function to fit quadratic model and extract peak DOY
fit_quadratic_peak <- function(data) {
  model <- lm(intensity ~ doy + I(doy^2), data = data)
  coefs <- coef(model)
  b <- coefs["doy"]
  c <- coefs["I(doy^2)"]
  peak_doy <- -b / (2 * c)
  rsq <- summary(model)$r.squared
  return(data.frame(peak_doy = peak_doy, r_squared = rsq))
}

# Apply to each year
results <- df %>%
  group_by(year) %>%
  group_modify(~fit_quadratic_peak(.x))

# Plot peak DOY over time
ggplot(results, aes(x = year, y = peak_doy)) +
  geom_line() +
  geom_point() +
  labs(title = "Estimated Peak Phenological Timing Over Years",
       x = "Year",
       y = "Peak Day of Year (DOY)") +
  theme_minimal()
