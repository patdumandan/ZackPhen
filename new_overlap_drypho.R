#overlap computation hierarchical slopes model
shared_years <- intersect(unique(dry_datA$year),
                          unique(pho_datA$Year))

N_shared = length(shared_years)

pl_shared_id = match(shared_years,sort(unique(dry_datA$year)))

ar_shared_id = match(shared_years,sort(unique(pho_datA$Year)))

#extract posterior draws
draws_pho=pho_mod$draws(variables=c("alpha", "beta_DOYs", "beta_DOYsqs"), format="df")
draws_dry=dry_mod2$draws(variables=c("alpha", "beta_DOYs", "beta_DOYsqs"), format="df")

alpha_pl        <- as_draws_matrix(draws_dry[,"alpha"])
beta_DOYs_pl    <- as_draws_matrix(draws_dry[, grepl("^beta_DOYs\\[", colnames(draws_dry))])
beta_DOYsqs_pl  <- as_draws_matrix(draws_dry[, grepl("^beta_DOYsqs\\[", colnames(draws_dry))])

alpha_ar        <- as_draws_matrix(draws_pho[,"alpha"])
beta_DOYs_ar    <- as_draws_matrix(draws_pho[, grepl("^beta_DOYs\\[", colnames(draws_pho))])
beta_DOYsqs_ar  <- as_draws_matrix(draws_pho[, grepl("^beta_DOYsqs\\[", colnames(draws_pho))])

# Standardized DOY sequences (plant + arthropod)
DOY_seq <- 150:270

doy_std_pl <- (DOY_seq - dryas_data$DOY_mean) / dryas_data$DOY_sd
doy_sq_std_pl <- doy_std_pl^2

doy_std_ar <- (DOY_seq - pho_data$DOY_mean) / pho_data$DOY_sd
doy_sq_std_ar <- doy_std_ar^2

n_draws <- nrow(alpha_pl)
Nyr_pl  <- ncol(beta_DOYs_pl)
Nyr_ar  <- ncol(beta_DOYs_ar)
ND      <- length(150:270)

#generate fitted curves
generate_fitted_curve <- function(draw, yr, alpha, beta_DOYs, beta_DOYsqs,
                                  z, z2, link) {

  eta <- alpha[draw] +
    beta_DOYs[draw, yr]   * z +
    beta_DOYsqs[draw, yr] * z2

  if (link == "logit") {
    out <- plogis(eta)       # probability scale
  } else if (link == "log") {
    out <- exp(eta)          # mean abundance
  } else {
    stop("Link must be 'log' or 'logit'")
  }

  return(out)
}

plant_fits <- array(NA, c(n_draws, Nyr_pl, ND))
arth_fits  <- array(NA, c(n_draws, Nyr_ar, ND))

for (draw in 1:n_draws) {
  for (y in 1:Nyr_pl) {
    plant_fits[draw, y, ] <- generate_fitted_curve(
      draw, y, alpha_pl, beta_DOYs_pl, beta_DOYsqs_pl,
      doy_std_pl, doy_sq_std_pl,
      link="logit")

    #remove for plotting with raw data
    # plant_fits[draw, y, ] <- plant_fits[draw, y, ] /
    #   sum(plant_fits[draw, y, ])
  }

  for (y in 1:Nyr_ar) {
    arth_fits[draw, y, ] <- generate_fitted_curve(
      draw, y, alpha_ar, beta_DOYs_ar, beta_DOYsqs_ar,
      doy_std_ar, doy_sq_std_ar,
      link="log")

    # arth_fits[draw, y, ] <- arth_fits[draw, y, ] /
    #   sum(arth_fits[draw, y, ])
  }
}

drypho_overlap <- matrix(NA, n_draws, N_shared)

for (s in 1:N_shared) {
  yp <- pl_shared_id[s]
  ya <- ar_shared_id[s]

  for (draw in 1:n_draws) {
    drypho_overlap[draw, s] <- sum(
      pmin(plant_fits[draw, yp, ], arth_fits[draw, ya, ])
    )
  }
}

summary_overlap=apply(drypho_overlap, 2, quantile, probs = c(0.025, 0.5, 0.975))%>%
  as.data.frame()%>%
  mutate(stat = c("lower", "median", "upper")) %>%
  pivot_longer(
    cols = -stat,
    names_to = "year",
    values_to = "overlap"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = overlap
  )

summary_overlap$year=shared_years

ggplot(summary_overlap, aes(x = year, y = median)) +
  geom_point(size = 3, color = "steelblue4") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, color = "steelblue4") +
  theme_classic() +
  labs(
    x = "Year",
    y = "Temporal Overlap",
    title = "Posterior Temporal Overlap (Dryas-Phoridae)")

#trend of overlap
head(drypho_overlap)
colnames(drypho_overlap)=shared_years

# vector of slopes
slopes <- numeric(n_draws)

for(d in 1:n_draws) {
  df <- data.frame(
    year = shared_years,
    overlap = drypho_overlap[d, ]
  )

  fit <- lm(overlap ~ year, data = df)
  slopes[d] <- coef(fit)["year"]
  slopes_per_decade=slopes*10
}

quantile(slopes_per_decade, probs = c(0.025, 0.5, 0.975))

#plot posteriors
drypho_df=as.data.frame(drypho_overlap)
drypho_df$draw <- 1:8000

drypho_long=drypho_df%>%
  pivot_longer(cols = -draw, names_to = "year", values_to = "overlap")

# Plot individual posterior draw trends
ggplot(drypho_long, aes(x = year, y = overlap, group = draw)) +
  geom_line(alpha = 0.05, color = "steelblue") +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "red", size = 1.2) +
  theme_classic() +
  labs(x = "Year", y = "Temporal overlap",
       title = "Posterior draws of temporal overlap")
