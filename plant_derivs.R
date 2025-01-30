# Load necessary libraries
library(lme4)
library(ggplot2)
library(dplyr)


# Get the predicted values
predicted_values <- predict(s1, type = "response", re.form = NA)

# Extract rand effects for Plot and Year
random_effects_plot <- ranef(s1)$Plot
random_effects_year <- ranef(s1)$year

# data frame with DOYs, years, plots, and the predicted values
prediction_data <- dry_datA %>%
  cbind(predicted = predicted_values)

# function to calculate the first derivative (numerical approximation)
first_derivative <- function(predicted_values, DOYs) {
  # Numerical derivative using finite difference method
  diff_pred <- diff(predicted_values)  # First differences
  diff_DOYs <- diff(DOYs)  # Differences in DOYs
  first_deriv <- diff_pred / diff_DOYs
  return(first_deriv)
}

# function to calculate the second derivative
second_derivative <- function(first_derivatives, DOYs) {
  # Numerical derivative of the first derivative using finite differences
  diff_first_deriv <- diff(first_derivatives)  # Second differences
  diff_DOYs <- diff(DOYs)[-1]  # Differences in DOYs for second derivative
  second_deriv <- diff_first_deriv / diff_DOYs
  return(second_deriv)
}

# Add random effects for Plot and Year to the data frame
prediction_data <- prediction_data %>%
  left_join(data.frame(Plot = unique(dry_datA$Plot), random_effect_plot = random_effects_plot), by = "Plot") %>%
  left_join(data.frame(year = unique(dry_datA$year), random_effect_year = random_effects_year), by = "year")

# data frame to hold the second derivatives for each Plot and Year combination

##blah not working right
derivatives_data <- prediction_data %>%
  arrange(Plot, year, DOY) %>%
  mutate(first_deriv = first_derivative(predicted, DOYs),
         second_deriv = second_derivative(first_deriv, DOYs))

# Plot
ggplot(derivatives_data, aes(x = DOYs, y = second_deriv, color = interaction(Plot, year))) +
  geom_line() +
  labs(x = "Day of Year (DOYs)", y = "Second Derivative", color = "Plot-Year Combination") +
  theme_minimal() +
  theme(legend.position = "bottom")
