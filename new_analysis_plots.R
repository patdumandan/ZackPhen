get_yearly_peaks <- function(model, doy_var = "DOYs", doy_sq_var = "DOYsqs", year_var = "yearc", raw_doy = NULL) {
  coefs <- fixef(model)
  year_levels <- levels(model@frame[[year_var]])
  base_year <- year_levels[1]

  b1 <- coefs[doy_var]
  b2 <- coefs[doy_sq_var]

  peak_std <- numeric(length(year_levels))
  names(peak_std) <- year_levels

  for (y in year_levels) {
    # Construct interaction coefficient names
    b1_y_name <- paste0(doy_var, ":", year_var, y)
    b2_y_name <- paste0(doy_sq_var, ":", year_var, y)

    # Fallback to 0 if interaction term does not exist (i.e., base year)
    b1_y <- if (b1_y_name %in% names(coefs)) coefs[b1_y_name] else 0
    b2_y <- if (b2_y_name %in% names(coefs)) coefs[b2_y_name] else 0

    # Compute year-specific linear and quadratic terms
    bt <- b1 + b1_y
    ct <- b2 + b2_y

    # Avoid division by zero
    if (!is.na(bt) && !is.na(ct) && ct != 0) {
      peak_std[y] <- -bt / (2 * ct)
    } else {
      peak_std[y] <- NA
    }
  }

  # Convert from standardized DOY back to actual DOY
  if (!is.null(raw_doy)) {
    peak_doy <- peak_std * sd(raw_doy) + mean(raw_doy)
    return(data.frame(year = as.numeric(year_levels), peak_doy = peak_doy, peak_std = peak_std))
  } else {
    return(data.frame(year = as.numeric(year_levels), peak_std = peak_std))
  }
}

peaks_df <- get_yearly_peaks(model = s1, raw_doy = dry_datA$DOY)

# Compute mean peak DOY
mean_peak <- mean(peaks_df$peak_doy, na.rm = TRUE)

# Fit linear model
lm_peak <- lm(peak_doy ~ year, data = peaks_df)
slope <- coef(lm_peak)[2]
pval <- summary(lm_peak)$coefficients[2, 4]

# Color by above/below average
cols <- ifelse(peaks_df$peak_doy >= mean_peak, "darkred", "blue")

# Base plot
plot(peaks_df$year, peaks_df$peak_doy, type = "n",
     xlab = "Year", ylab = "Estimated Peak DOY",
     main = "Year-Specific Peak Flowering Dates (with Trend)")

# Add points colored by above/below average
points(peaks_df$year, peaks_df$peak_doy, col = cols, pch = 19)

# Add connecting line
lines(peaks_df$year, peaks_df$peak_doy, col = "gray40", lwd = 1)

# Add horizontal mean line
abline(h = mean_peak, col = "black", lty = 3)

# Add linear trend line
abline(lm_peak, col = "black", lwd = 2)

# Determine central text position
text_x <- mean(range(peaks_df$year, na.rm = TRUE))
text_y <- mean(range(peaks_df$peak_doy, na.rm = TRUE))

# Add annotation text at center
text(text_x, text_y,
     labels = paste0("Trend: ", round(slope, 3), " DOY/year\np = ", signif(pval, 3)),
     col = "black", cex = 0.95, font = 2)

# Add legend
legend("topleft", legend = c("Above avg", "Below avg", "Trend", "Mean"),
       col = c("darkred", "blue", "black", "black"),
       pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 3), lwd = c(NA, NA, 2, 1),
       bty = "n")


plot_obs_vs_fit_raw_actual_DOY <- function(model, data, years_to_plot,
                                           doy_var_std = "DOYs", doy_raw_var = "DOY") {
  library(dplyr)
  library(ggplot2)

  data_filtered <- data %>%
    filter(yearc %in% years_to_plot)

  # Add squared standardized DOY if missing
  if (!"DOYsqs" %in% names(data_filtered)) {
    data_filtered$DOYsqs <- data_filtered[[doy_var_std]]^2
  }

  # Compute observed flowering proportion
  data_filtered <- data_filtered %>%
    mutate(obs_flowering = tot_flwr / (tot_flwr + tot_NF))

  # Predict fitted values from model
  data_filtered$predicted <- predict(model, newdata = data_filtered, type = "response", re.form = NA)

  # Plot using actual DOY
  ggplot(data_filtered, aes_string(x = doy_raw_var)) +
    geom_point(aes(y = obs_flowering, color = Plot), alpha = 0.6, size = 2) +
    geom_smooth(aes(y = predicted), color = "black", method = "loess", se = FALSE, size = 1) +
    facet_wrap(~yearc, scales = "free_y") +
    labs(title = "Observed vs. Fitted Flowering by Year (Colored by Plot)",
         x = "Actual DOY",
         y = "Proportion Flowering",
         color = "Plot") +
    theme_classic() +
    theme(legend.position = "right")
}

plot_obs_vs_fit_raw_actual_DOY(model = s1, data = dry_datA, years_to_plot = levels(dry_datA$yearc))

extract_shape_terms <- function(model, doy_var = "DOYs", doy_sq_var = "DOYsqs",
                                year_var = "yearc") {
  coefs <- fixef(model)
  year_levels <- levels(model@frame[[year_var]])
  base_year <- year_levels[1]

  b1_main <- coefs[doy_var]
  b2_main <- coefs[doy_sq_var]

  # Initialize data frame
  shape_df <- data.frame(
    year = as.numeric(year_levels),
    bt = NA,
    ct = NA
  )

  for (y in year_levels) {
    # Interaction term names
    b1_y_name <- paste0(doy_var, ":", year_var, y)
    b2_y_name <- paste0(doy_sq_var, ":", year_var, y)

    # If interaction not present (base year), default to 0
    b1_y <- if (b1_y_name %in% names(coefs)) coefs[b1_y_name] else 0
    b2_y <- if (b2_y_name %in% names(coefs)) coefs[b2_y_name] else 0

    # Year-specific slope and curvature
    bt <- b1_main + b1_y
    ct <- b2_main + b2_y

    shape_df[shape_df$year == as.numeric(y), "bt"] <- bt
    shape_df[shape_df$year == as.numeric(y), "ct"] <- ct
  }

  # Add ratio and peak location
  shape_df <- shape_df %>%
    mutate(bt_ct_ratio = bt / ct,
           peak_std = -bt / (2 * ct))

  return(shape_df)
}

shape_df <- extract_shape_terms(model = s1)

flat_threshold <- -0.05   # near-zero = flat curve
positive_ct <- shape_df$ct > 0
flat_ct <- shape_df$ct > flat_threshold

# Plot with annotations
# Fit linear model to curvature trend
ct_lm <- lm(ct ~ year, data = shape_df)
ct_slope <- coef(ct_lm)[2]
ct_pval <- summary(ct_lm)$coefficients[2, 4]

# Base plot
plot(shape_df$year, shape_df$ct, type = "l", col = "red", lwd = 2,
     ylab = "Curvature (ct)", xlab = "Year", main = "Phenological Shapes Over Time")

# Highlight years 1997 and 2021
highlight_years <- c(2018)
points(shape_df$year[shape_df$year %in% highlight_years],
     shape_df$ct[shape_df$year %in% highlight_years],
    col = "darkgreen", pch = 19, cex = 1.5)
text(shape_df$year[shape_df$year %in% highlight_years],
     shape_df$ct[shape_df$year %in% highlight_years],
     labels = highlight_years, pos = 3, cex = 0.9)

# Highlight years with ct > 0
positive_ct <- shape_df$ct > 0
points(shape_df$year[positive_ct], shape_df$ct[positive_ct],
       col = "purple", pch = 4, cex = 1.5)
text(shape_df$year[positive_ct], shape_df$ct[positive_ct],
     labels = shape_df$year[positive_ct], pos = 3, col = "purple", cex = 0.75)

# Horizontal reference lines
abline(h = 0, lty = 2, col = "black")       # ct = 0
abline(h = -0.05, lty = 3, col = "gray40")  # flat-ish threshold

# Linear trend line
abline(ct_lm, col = "black", lwd = 2)

# Add slope and p-value
text_x <- min(shape_df$year) + 1
text_y <- max(shape_df$ct, na.rm = TRUE) * 0.9
text(text_x, text_y,
     labels = paste0("Trend: ", round(ct_slope, 4),
                     " /yr\np = ", signif(ct_pval, 3)),
     col = "black", adj = 0)

# Legend
legend("bottomleft", legend = c("ct", "Trend line",  "ct > 0"),
       col = c("red", "blue", "purple"),
       pch = c(NA, NA, 19), lty = c(1, 1, NA), lwd = 2, bty = "n")

# Base plot
plot(peaks_df$year, peaks_df$peak_doy, type = "n",
     xlab = "Year", ylab = "Estimated Peak DOY",
     main = "Year-Specific Peak Flowering Dates (with Trend)")

# Add points colored by above/below average
points(peaks_df$year, peaks_df$peak_doy, col = cols, pch = 19)

# Add connecting line
lines(peaks_df$year, peaks_df$peak_doy, col = "gray40", lwd = 1)

# Add horizontal mean line
abline(h = mean_peak, col = "black", lty = 3)

# Add linear trend line
abline(lm_peak, col = "black", lwd = 2)

# Highlight the 2018 point
highlight_index <- which(peaks_df$year == 2018)
if (length(highlight_index) == 1) {
  points(peaks_df$year[highlight_index], peaks_df$peak_doy[highlight_index],
         col = "red", pch = 21, bg = "yellow", cex = 1.5, lwd = 2)
}

# Determine central text position
text_x <- mean(range(peaks_df$year, na.rm = TRUE))
text_y <- mean(range(peaks_df$peak_doy, na.rm = TRUE))

# Add annotation text at center
text(text_x, text_y,
     labels = paste0("Trend: ", round(slope, 3), " DOY/year\np = ", signif(pval, 3)),
     col = "black", cex = 0.95, font = 2)

