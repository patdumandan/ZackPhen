#packages####
library(bayesplot)
library(posterior)

DOY <- seq(-2, 2, length = 100)

#Dryas####

drymod_draws <- dry_mod4b$draws()%>%posterior::as_draws_array()

varnames <- variables(drymod_draws)
mu_params <- grep("^mu\\[", varnames, value = TRUE)
width_params <- grep("^width\\[", varnames, value = TRUE)

mcmc_trace(drymod_draws, pars = mu_params)
mcmc_trace(drymod_draws, pars = width_params)
mcmc_trace(drymod_draws, pars = c("sigma_mu","sigma_width", "sigma_year"))

mu_draws= dry_mod4b$draws("mu_year")
mcmc_areas(mu_draws, regex = "mu_year")

width_draws= dry_mod4b$draws("width")
mcmc_areas(width_draws, regex = "width")

# Rhat values
rhat_val=dry_mod4b$summary()$rhat
rhat_vals <- posterior::rhat(drymod_draws)
ess_vals <- posterior::ess_bulk(drymod_draws)   # numeric vector
n_draws  <- posterior::ndraws(drymod_draws)              # total draws per chain × chains

neff_ratio_vals <- ess_vals / n_draws

p1 <- mcmc_rhat_hist(rhat_val, bins=10)
p2 <- mcmc_rhat(rhat_val)
p3 <- mcmc_neff_hist(neff_ratio_vals)

ggarrange(p1, p2, p3, nrow = 3)

# Regular, visual, posterior predictive checks
color_scheme_set("red")
yrep = dry_mod4b$draws(format="draws_matrix", variables="y_pred")
yrep <- yrep[(nrow(yrep)-99):nrow(yrep), ]  # replicate data, take only the last 100 samples to speed up drawing
ggarrange(ppc_dens_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_ecdf_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
          nrow = 2, ncol=2)
ppc_hist(y = dry_datA$tot_F, yrep[1:11,])

mu         <- dry_mod4b$draws("mu_year", format="draws_matrix")
width      <- dry_mod4b$draws("width", format="draws_matrix")
alpha_year <- dry_mod4b$draws("alpha_year", format="draws_matrix")
u_plot     <- dry_mod4b$draws("u_plot", format="draws_matrix")
u_plot_mu     <- dry_mod4b$draws("u_plot_mu", format="draws_matrix")

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
          (DOY[nd]-(mu[, y] + u_plot_mu[, p_idx]))^2/(2 * width[, y]^2)

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
