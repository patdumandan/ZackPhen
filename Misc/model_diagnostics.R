#packages####
library(bayesplot)
library(posterior)

#Dryas####

mu_draws= dry_mod6b1$draws("mu_year")
mcmc_areas(mu_draws, regex = "mu_year")

width_draws= dry_mod6b1$draws("width")
mcmc_areas(width_draws, regex = "width")

drymod6b1_draws <- dry_mod6b1$draws()%>%posterior::as_draws_array()

varnames <- variables(drymod_draws)
mu_params <- grep("^mu\\[", varnames, value = TRUE)
width_params <- grep("^width\\[", varnames, value = TRUE)

mcmc_trace(drymod_draws, pars = mu_params)
mcmc_trace(drymod_draws, pars = width_params)
mcmc_trace(drymod_draws, pars = c("sigma_mu","sigma_width", "sigma_year"))

# Rhat values
ess_vals <- posterior::ess_bulk(drymod6b1_draws)
n_draws  <- posterior::ndraws(drymod6b1_draws)

neff_ratio_vals <- ess_vals / n_draws

p1= mcmc_rhat_hist(dry_mod6b1$summary()$rhat)
p2= mcmc_rhat(dry_mod6b1$summary()$rhat)
p3= mcmc_neff_hist(neff_ratio_vals, bins=10)

ggarrange(p1, p2, p3, ncol = 3)

# Regular, visual, posterior predictive checks
color_scheme_set("blue")
yrep = dry_mod6b1$draws(format="draws_matrix", variables="y_pred")
yrep <- yrep[(nrow(yrep)-99):nrow(yrep), ]  # replicate data, take only the last 100 samples to speed up drawing
ggarrange(ppc_dens_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_ecdf_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
          nrow = 2, ncol=2)
ppc_hist(y = dry_datA$tot_F, yrep[1:11,])

mu         <- dry_mod6b1$draws("mu", format="draws_matrix")
width      <- dry_mod6b1$draws("width", format="draws_matrix")
alpha_year <- dry_mod6b1$draws("alpha_year", format="draws_matrix")
u_plot     <- dry_mod6b1$draws("u_plot", format="draws_matrix")
u_plot_mu  <- dry_mod6b1$draws("u_plot_mu", format="draws_matrix")

#per year####
DOY = seq(-2, 2, length = 100)
invLogit <- function(x) exp(x) / (1 + exp(x))

plot.ids <- sort(unique(dryas_data$plot_id))
P <- length(plot.ids)
cols <- rainbow(P)

par(mfrow = c(3,4))

eta_post <- numeric(length(DOY)) #posterior preds

for (y in 1:ncol(mu)) {            # years
  for (pl in 1:P) {                # plots
    p_idx <- plot.ids[pl]
    eta_hat <- numeric(length(DOY))

    for (nd in 1:length(DOY)) {
      # compute eta for each draw, then average

      eta_post= alpha_year[, y]+u_plot[, p_idx] -
          (DOY[nd]-(mu[, y] + u_plot_mu[, p_idx]))^2/(width[, y]^2)

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    if (pl == 1) {
      plot(DOY, eta_hat, type="l", col=cols[pl],
           ylim=c(0,1), main=paste("Year", y))
    } else {
      lines(DOY, eta_hat, col=cols[pl])
    }

    ind <- dryas_data$plot_id == p_idx & dryas_data$year_id == y
    points(dryas_data$DOYs[ind],
           dryas_data$tot_F[ind] / (dryas_data$tot_F[ind] + dryas_data$tot_NF[ind]),
           col=cols[pl])
  }
}

#per plot####
DOY <- seq(-2, 2, length = 100)
invLogit <- function(x) exp(x) / (1 + exp(x))

year.ids <- sort(unique(dryas_data$year_id))
YR <- length(year.ids)
cols <- rainbow(YR)

plot.ids <- sort(unique(dryas_data$plot_id))
P <- length(plot.ids)

par(mfrow = c(3, 2))

for (pl in plot.ids) {                 # LOOP OVER PLOTS

  plot(NULL, xlim = range(DOY), ylim = c(0, 1),
       xlab = "DOY", ylab = "Proportion F",
       main = paste("Plot", pl))

  for (y in 1:YR) {                    # LOOP OVER YEARS

    yr_idx <- year.ids[y]
    eta_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {

      eta_post <- alpha_year[, y] +
        u_plot[, pl] -
        (DOY[nd] - (mu[, y] + u_plot_mu[, pl]))^2 /
        (width[, y]^2)

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    # posterior mean line for this year
    lines(DOY, eta_hat, col = cols[y], lwd = 2)

    # observed data
    ind <- dryas_data$plot_id == pl &
      dryas_data$year_id == yr_idx

    points(dryas_data$DOYs[ind],
           dryas_data$tot_F[ind] /
             (dryas_data$tot_F[ind] + dryas_data$tot_NF[ind]),
           col = cols[y], pch = 16)
  }
}

#diagnostics for plant_phen7
drymod7_draws <- dry_mod7$draws()%>%posterior::as_draws_array()
ess_vals <- posterior::ess_bulk(drymod7_draws)   # numeric vector
n_draws  <- posterior::ndraws(drymod7_draws)              # total draws per chain × chains

neff_ratio_vals <- ess_vals / n_draws

p1= mcmc_rhat_hist(dry_mod7$summary()$rhat)
p2= mcmc_rhat(dry_mod7$summary()$rhat)
p3= mcmc_neff_hist(neff_ratio_vals, bins=10)

ggarrange(p1, p2, p3, ncol = 3)

color_scheme_set("blue")
yrep = dry_mod7$draws(format="draws_matrix", variables="y_pred")
yrep <- yrep[(nrow(yrep)-99):nrow(yrep), ]  # replicate data, take only the last 100 samples to speed up drawing
ggarrange(ppc_dens_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_ecdf_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
          nrow = 2, ncol=2)
mu              <- dry_mod7$draws("mu", format = "draws_matrix")
alpha_year      <- dry_mod7$draws("alpha_year", format = "draws_matrix")
u_plot          <- dry_mod7$draws("u_plot", format = "draws_matrix")
u_plot_mu       <- dry_mod7$draws("u_plot_mu", format = "draws_matrix")

width_L_year    <- dry_mod7$draws("width_L_year", format = "draws_matrix")
width_R_year    <- dry_mod7$draws("width_R_year", format = "draws_matrix")

width_L_plot    <- dry_mod7$draws("width_L_plot", format = "draws_matrix")
width_R_plot    <- dry_mod7$draws("width_R_plot", format = "draws_matrix")

DOY <- seq(-2, 2, length = 150)
invLogit <- plogis

plot.ids <- sort(unique(dryas_data$plot_id))
P <- length(plot.ids)

cols <- rainbow(P)

par(mfrow = c(3, 4))
for (y in 1:ncol(mu)) {          # loop over years
  for (pl in 1:P) {              # loop over plots

    p_idx <- plot.ids[pl]
    eta_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {

      # distance from peak (per posterior draw)
      d <- DOY[nd] - (mu[, y] + u_plot_mu[, p_idx])

      # year × plot widths
      width_L_np <- width_L_year[, y] * width_L_plot[, p_idx]
      width_R_np <- width_R_year[, y] * width_R_plot[, p_idx]

      penalty <- ifelse(
        d < 0,
        d^2 / (2 * width_L_np^2),
        d^2 / (2 * width_R_np^2)
      )

      eta_post <-
        alpha_year[, y] +
        u_plot[, p_idx] -
        penalty

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    # draw curves
    if (pl == 1) {
      plot(
        DOY, eta_hat, type = "l",
        ylim = c(0, 1),
        col = cols[pl],
        main = paste("Year", y),
        xlab = "DOY", ylab = "Probability"
      )
    } else {
      lines(DOY, eta_hat, col = cols[pl])
    }

    # observed data
    ind <- dryas_data$plot_id == p_idx &
      dryas_data$year_id == y

    points(
      dryas_data$DOYs[ind],
      dryas_data$tot_F[ind] /
        (dryas_data$tot_F[ind] + dryas_data$tot_NF[ind]),
      col = cols[pl],
      pch = 16
    )
  }
}

#per plot#####
DOY <- seq(-2, 2, length = 150)
invLogit <- plogis

year.ids <- sort(unique(dryas_data$year_id))
YR <- length(year.ids)

plot.ids <- sort(unique(dryas_data$plot_id))
P <- length(plot.ids)

cols <- rainbow(YR)

par(mfrow = c(3, 2))

for (pl in plot.ids) {   # LOOP OVER PLOTS

  plot(
    NULL,
    xlim = range(DOY),
    ylim = c(0, 1),
    xlab = "DOY",
    ylab = "Proportion F",
    main = paste("Plot", pl)
  )

  for (y in seq_along(year.ids)) {  # LOOP OVER YEARS

    yr_idx <- year.ids[y]
    eta_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {

      # distance from peak (vector over posterior draws)
      d <- DOY[nd] - (mu[, y] + u_plot_mu[, pl])

      # year × plot widths
      width_L_np <- width_L_year[, y] * width_L_plot[, pl]
      width_R_np <- width_R_year[, y] * width_R_plot[, pl]

      penalty <- ifelse(
        d < 0,
        d^2 / (2 * width_L_np^2),
        d^2 / (2 * width_R_np^2)
      )

      eta_post <-
        alpha_year[, y] +
        u_plot[, pl] -
        penalty

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    # posterior mean curve
    lines(DOY, eta_hat, col = cols[y], lwd = 2)

    # observed data
    ind <- dryas_data$plot_id == pl &
      dryas_data$year_id == yr_idx

    points(
      dryas_data$DOYs[ind],
      dryas_data$tot_F[ind] /
        (dryas_data$tot_F[ind] + dryas_data$tot_NF[ind]),
      col = cols[y],
      pch = 16
    )
  }
}

#saxifraga####

saxmod6b_draws <- sax_mod6b$draws()%>%posterior::as_draws_array()

ess_vals <- posterior::ess_bulk(saxmod6b_draws)
n_draws  <- posterior::ndraws(saxmod6b_draws)

neff_ratio_vals <- ess_vals / n_draws

p1= mcmc_rhat_hist(sax_mod6b$summary()$rhat)
p2= mcmc_rhat(sax_mod6b$summary()$rhat)
p3= mcmc_neff_hist(neff_ratio_vals, bins=10)

ggarrange(p1, p2, p3, ncol = 3)

# Regular, visual, posterior predictive checks
color_scheme_set("blue")
yrep = sax_mod6b$draws(format="draws_matrix", variables="y_pred")
yrep <- yrep[(nrow(yrep)-99):nrow(yrep), ]  # replicate data, take only the last 100 samples to speed up drawing
ggarrange(ppc_dens_overlay(y = sax_datA$tot_F,
                           yrep = yrep),
          ppc_ecdf_overlay(y = sax_datA$tot_F,
                           yrep = yrep),
          ppc_pit_ecdf(y = sax_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
          ppc_pit_ecdf(y = sax_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
          nrow = 2, ncol=2)


mu_draws= sax_mod6b$draws("mu")
mcmc_areas(mu_draws, regex = "mu")

width_draws= sax_mod6b$draws("width")
mcmc_areas(width_draws, regex = "width")

mu         <- sax_mod6b$draws("mu", format="draws_matrix")
width      <- sax_mod6b$draws("width", format="draws_matrix")
alpha_year <- sax_mod6b$draws("alpha_year", format="draws_matrix")
u_plot     <- sax_mod6b$draws("u_plot", format="draws_matrix")
u_plot_mu  <- sax_mod6b$draws("u_plot_mu", format="draws_matrix")

#per year####
DOY = seq(-2, 2, length = 100)
invLogit <- function(x) exp(x) / (1 + exp(x))

plot.ids <- sort(unique(saxifraga_data$plot_id))
P <- length(plot.ids)
cols <- rainbow(P)

par(mfrow = c(3,4))

eta_post <- numeric(length(DOY)) #posterior preds

for (y in 1:ncol(mu)) {            # years
  for (pl in 1:P) {                # plots
    p_idx <- plot.ids[pl]
    eta_hat <- numeric(length(DOY))

    for (nd in 1:length(DOY)) {
      # compute eta for each draw, then average

      eta_post= alpha_year[, y]+u_plot[, p_idx] -
        (DOY[nd]-(mu[, y] + u_plot_mu[, p_idx]))^2/(width[, y]^2)

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    if (pl == 1) {
      plot(DOY, eta_hat, type="l", col=cols[pl],
           ylim=c(0,1), main=paste("Year", y))
    } else {
      lines(DOY, eta_hat, col=cols[pl])
    }

    ind <- saxifraga_data$plot_id == p_idx & saxifraga_data$year_id == y
    points(saxifraga_data$DOYs[ind],
           saxifraga_data$tot_F[ind] / (saxifraga_data$tot_F[ind] + saxifraga_data$tot_NF[ind]),
           col=cols[pl])
  }
}

#per plot####
DOY <- seq(-2, 2, length = 100)
invLogit <- function(x) exp(x) / (1 + exp(x))

year.ids <- sort(unique(saxifraga_data$year_id))
YR <- length(year.ids)
cols <- rainbow(YR)

plot.ids <- sort(unique(saxifraga_data$plot_id))
P <- length(plot.ids)

par(mfrow = c(3, 1))

for (pl in plot.ids) {                 # LOOP OVER PLOTS

  plot(NULL, xlim = range(DOY), ylim = c(0, 1),
       xlab = "DOY", ylab = "Proportion F",
       main = paste("Plot", pl))

  for (y in 1:YR) {                    # LOOP OVER YEARS

    yr_idx <- year.ids[y]
    eta_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {

      eta_post <- alpha_year[, y] +
        u_plot[, pl] -
        (DOY[nd] - (mu[, y] + u_plot_mu[, pl]))^2 /
        (width[, y]^2)

      eta_hat[nd] <- mean(invLogit(eta_post))
    }

    # posterior mean line for this year
    lines(DOY, eta_hat, col = cols[y], lwd = 2)

    # observed data
    ind <- saxifraga_data$plot_id == pl &
      saxifraga_data$year_id == yr_idx

    points(saxifraga_data$DOYs[ind],
           saxifraga_data$tot_F[ind] /
             (saxifraga_data$tot_F[ind] + saxifraga_data$tot_NF[ind]),
           col = cols[y], pch = 16)
  }
}

#muscid####
mu         <- mus_mod$draws("mu", format="draws_matrix")
width      <- mus_mod$draws("width", format="draws_matrix")
alpha_year <- mus_mod$draws("alpha_year", format="draws_matrix")
u_plot     <- mus_mod$draws("u_plot", format="draws_matrix")
u_plot_mu  <- mus_mod$draws("u_plot_mu", format="draws_matrix")

DOY <- seq(-2, 2, length = 100)
obs_days <- mus_data$obs_days

plot.ids <- sort(unique(mus_data$plot_id))
P <- length(plot.ids)
cols <- rainbow(P)

par(mfrow = c(3, 4))

for (y in 1:ncol(mu)) {

  ymax <- 0

  for (pl in 1:P) {
    p_idx <- plot.ids[pl]

    mu_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {
      eta_post <-
        exp(
          alpha_year[, y] +
            u_plot[, p_idx] -
            (DOY[nd] - (mu[, y] + u_plot_mu[, p_idx]))^2 /
            (width[, y]^2) +
            log(obs_days)
        )

      mu_hat[nd] <- mean(eta_post)
    }

    ymax <- max(ymax, mu_hat, na.rm = TRUE)

    # observed values for this plot–year
    ind <- mus_data$plot_id == p_idx &
      mus_data$year_id == y
    ymax <- max(ymax, mus_data$y[ind], na.rm = TRUE)
  }

  ylim <- c(0, ymax * 1.1)

  first_plot <- TRUE

  for (pl in 1:P) {

    p_idx <- plot.ids[pl]
    mu_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {
      eta_post <-
        exp(
          alpha_year[, y] +
            u_plot[, p_idx] -
            (DOY[nd] - (mu[, y] + u_plot_mu[, p_idx]))^2 /
            (width[, y]^2) +
            log(obs_days)
        )

      mu_hat[nd] <- mean(eta_post)
    }

    if (first_plot) {
      plot(
        DOY, mu_hat,
        type = "l",
        col = cols[pl],
        lwd = 2,
        ylim = ylim,
        main = paste("Year", y),
        xlab = "DOY",
        ylab = "Expected count"
      )
      first_plot <- FALSE
    } else {
      lines(DOY, mu_hat, col = cols[pl], lwd = 2)
    }

    # observed data
    ind <- mus_data$plot_id == p_idx &
      mus_data$year_id == y

    points(
      mus_data$DOYs[ind],
      mus_data$y[ind],
      col = cols[pl],
      pch = 16
    )
  }
}

#per plot####
mu         <- mus_mod$draws("mu", format="draws_matrix")
width      <- mus_mod$draws("width", format="draws_matrix")
alpha_year <- mus_mod$draws("alpha_year", format="draws_matrix")
u_plot     <- mus_mod$draws("u_plot", format="draws_matrix")
u_plot_mu  <- mus_mod$draws("u_plot_mu", format="draws_matrix")

DOY <- seq(-2, 2, length = 100)
obs_days <- mus_data$obs_days

year.ids <- sort(unique(mus_data$year_id))
Y <- length(year.ids)

plot.ids <- sort(unique(mus_data$plot_id))
P <- length(plot.ids)

cols <- rainbow(Y)

par(mfrow = c(2, 2))

for (pl in plot.ids) {
  ymax <- 0

  for (y in 1:Y) {

    mu_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {
      eta_post <-
        exp(
          alpha_year[, y] +
            u_plot[, pl] -
            (DOY[nd] - (mu[, y] + u_plot_mu[, pl]))^2 /
            (width[, y]^2) +
            log(obs_days)
        )

      mu_hat[nd] <- mean(eta_post)
    }

    ymax <- max(ymax, mu_hat, na.rm = TRUE)

    ind <- mus_data$plot_id == pl &
      mus_data$year_id == y
    ymax <- max(ymax, mus_data$y[ind], na.rm = TRUE)
  }

  ylim <- c(0, ymax * 1.1)

  first_year <- TRUE

  for (y in 1:Y) {

    mu_hat <- numeric(length(DOY))

    for (nd in seq_along(DOY)) {
      eta_post <-
        exp(
          alpha_year[, y] +
            u_plot[, pl] -
            (DOY[nd] - (mu[, y] + u_plot_mu[, pl]))^2 /
            (width[, y]^2) +
            log(obs_days)
        )

      mu_hat[nd] <- mean(eta_post)
    }

    if (first_year) {
      plot(
        DOY, mu_hat,
        type = "l",
        col = cols[y],
        lwd = 2,
        ylim = ylim,
        main = paste("Plot", pl),
        xlab = "DOY",
        ylab = "Expected count"
      )
      first_year <- FALSE
    } else {
      lines(DOY, mu_hat, col = cols[y], lwd = 2)
    }

    # observed data for this plot–year
    ind <- mus_data$plot_id == pl &
      mus_data$year_id == y

    points(
      mus_data$DOYs[ind],
      mus_data$y[ind],
      col = cols[y],
      pch = 16
    )
  }
}
