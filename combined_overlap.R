library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)

drypho=list(
DOY_mean_pl=mean(plant_datA$DOY),
DOY_sd_pl=sd(plant_datA$DOY),
Nyr_pl=length(unique(dry_datA$year)),
Nplant=nrow(dry_datA),
Nplant_plots=length(unique(dry_datA$plot_id)),
plant_plot_id=dry_datA$plot_id,
plant_year_id=as.integer(factor(dry_datA$year)),
tot_F=dry_datA$tot_F,
tot_NF=dry_datA$tot_NF,
plant_DOY=dry_datA$DOYs,
plant_DOYsqs=dry_datA$DOYsqs,
DOY_mean_ar=mean(arth_df$DOY),
DOY_sd_ar=sd(arth_df$DOY),
Nyr_ar=length(unique(pho_datA$Year)),
Narth=nrow(pho_datA),
Narth_plots=length(unique(pho_datA$plot_id)),
arth_plot_id=pho_datA$plot_id,
arth_year_id=as.integer(factor(pho_datA$Year)),
arth_y=pho_datA$TotalCatch1,
arth_effort=pho_datA$TrapDays,
arth_DOY=pho_datA$DOYs,
arth_DOYsqs=pho_datA$DOYsqs,
ND=121)

#later for overlap computation:
shared_years <- intersect(unique(dry_datA$year),
                          unique(pho_datA$Year))

drypho$N_shared = length(shared_years)

drypho$pl_shared_id = match(shared_years,
                            sort(unique(dry_datA$year)))

drypho$ar_shared_id = match(shared_years,
                            sort(unique(pho_datA$Year)))


#fit model
combined_mod=cmdstan_model("combined_plantarth.stan")

drypho_mod <- combined_mod$sample(
  data = drypho,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000,
  iter_warmup = 500)

#extract posterior draws
draws_drypho=drypho_mod$draws(variables=c("alpha_pl", "beta_DOYs_pl", "beta_DOYsqs_pl",
                                           "alpha_ar", "beta_DOYs_ar", "beta_DOYsqs_ar"),
                               format="df")

alpha_pl        <- as_draws_matrix(draws_drypho[,"alpha_pl"])
beta_DOYs_pl    <- as_draws_matrix(draws_drypho[,grepl("beta_DOYs_pl", colnames(draws_drypho))])
beta_DOYsqs_pl  <- as_draws_matrix(draws_drypho[,grepl("beta_DOYsqs_pl", colnames(draws_drypho))])

alpha_ar        <- as_draws_matrix(draws_drypho[,"alpha_ar"])
beta_DOYs_ar    <- as_draws_matrix(draws_drypho[,grepl("beta_DOYs_ar", colnames(draws_drypho))])
beta_DOYsqs_ar  <- as_draws_matrix(draws_drypho[,grepl("beta_DOYsqs_ar", colnames(draws_drypho))])

# Standardized DOY sequences (plant + arthropod)
DOY_seq <- 150:270

doy_std_pl <- (DOY_seq - drypho$DOY_mean_pl) / drypho$DOY_sd_pl
doy_sq_std_pl <- doy_std_pl^2

doy_std_ar <- (DOY_seq - drypho$DOY_mean_ar) / drypho$DOY_sd_ar
doy_sq_std_ar <- doy_std_ar^2

n_draws <- nrow(alpha_pl)
Nyr_pl  <- ncol(beta_DOYs_pl)
Nyr_ar  <- ncol(beta_DOYs_ar)
ND      <- length(150:270)

generate_fitted_curve <- function(draw, yr, alpha, beta_DOYs, beta_DOYsqs, z, z2) {
  eta <- alpha[draw] +
    beta_DOYs[draw, yr]   * z + #z=standardized DOY seq, z2=quad for DOY seq
    beta_DOYsqs[draw, yr] * z2
  plogis(eta)
}

plant_fits <- array(NA, c(n_draws, Nyr_pl, ND))
arth_fits  <- array(NA, c(n_draws, Nyr_ar, ND))

for (draw in 1:n_draws) {
  for (y in 1:Nyr_pl) {
    plant_fits[draw, y, ] <- generate_fitted_curve(
      draw, y, alpha_pl, beta_DOYs_pl, beta_DOYsqs_pl,
      doy_std_pl, doy_sq_std_pl)

    s <- sum(plant_fits[draw, y, ]) #normalize values to get PDF
    if (s > 0) plant_fits[draw, y, ] <- plant_fits[draw, y, ] / s
  }

  for (y in 1:Nyr_ar) {
    arth_fits[draw, y, ] <- generate_fitted_curve(
      draw, y, alpha_ar, beta_DOYs_ar, beta_DOYsqs_ar,
      doy_std_ar, doy_sq_std_ar)

    s <- sum(arth_fits[draw, y, ])
    if (s > 0) arth_fits[draw, y, ] <- arth_fits[draw, y, ] / s
  }
}

#calculate overlap

N_shared <- length(drypho$pl_shared_id)
pl_shared_id=drypho$pl_shared_id
ar_shared_id=drypho$ar_shared_id
n_draws <- nrow(alpha_pl)

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
