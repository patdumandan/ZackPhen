#packages####
library(bayesplot)
library(posterior)

DOY <- seq(-2, 2, length = 100)

#Dryas####

drymod_draws <- dry_mod6c$draws()%>%posterior::as_draws_array()

varnames <- variables(drymod_draws)
mu_params <- grep("^mu\\[", varnames, value = TRUE)
width_params <- grep("^width\\[", varnames, value = TRUE)

mcmc_trace(drymod_draws, pars = mu_params)
mcmc_trace(drymod_draws, pars = width_params)
mcmc_trace(drymod_draws, pars = c("sigma_mu","sigma_width", "sigma_year"))

plot(drymod_draws, pars="mu")     # Peak dates
plot(drymod_draws, pars="width")  # width parameter of the phenology

# Rhat values
rhat_vals <- posterior::rhat(drymod_draws)
ess_vals <- posterior::ess_bulk(drymod_draws)   # numeric vector
n_draws  <- posterior::ndraws(drymod_draws)              # total draws per chain × chains

neff_ratio_vals <- ess_vals / n_draws

p1 <- mcmc_rhat_hist(rhat_vals)
p2 <- mcmc_rhat(rhat_vals)
p3 <- mcmc_neff_hist(neff_ratio_vals)

ggarrange(p1, p2, p3, nrow = 3)

# Regular, visual, posterior predictive checks
color_scheme_set("red")
yrep = dry_mod$draws(format="draws_matrix", variables="y_pred")
yrep <- yrep[(nrow(yrep)-99):nrow(yrep), ]  # replicate data, take only the last 100 samples to speed up drawing
ggarrange(ppc_dens_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_ecdf_overlay(y = dry_datA$tot_F,
                           yrep = yrep),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = FALSE),
          ppc_pit_ecdf(y = dry_datA$tot_F, yrep = yrep, prob = 0.99, plot_diff = TRUE),
          nrow = 2, ncol=2)
ppc_hist(y = dry_datA$tot_F, yrep[1:11,])

mu         <- dry_mod$draws("mu", format="draws_matrix")
width      <- dry_mod$draws("width", format="draws_matrix")
alpha_year <- dry_mod$draws("alpha_year", format="draws_matrix")
u_plot     <- dry_mod$draws("u_plot", format="draws_matrix")

invLogit <- function(x) exp(x) / (1 + exp(x))

plot.id <- sort(unique(dryas_data$plot_id))
P <- length(plot.id)
cols <- rainbow(P)

par(mfrow = c(3,4))

eta <- numeric(length(DOY))

for (y in 1:ncol(mu)) {
  for (pl in 1:P) {
    for (nd in 1:length(DOY)) {
      eta[nd] <- mean(invLogit(alpha_year[, y] -(DOY[nd] - mu[, y])^2 / width[, y]^2 +
            u_plot[, plot.id[pl]]))
    }
    if (pl == 1) {
      plot(DOY, eta, type="l", col=cols[pl],
           ylim=c(0,1),
           main=paste("Year", y))
    } else {
      lines(DOY, eta, col=cols[pl])
    }
    ind <- dryas_data$plot_id == plot.id[pl] & dryas_data$year_id == y
    points(dryas_data$DOYs[ind],
           dryas_data$tot_F[ind]/(dryas_data$tot_F[ind] + dryas_data$tot_NF[ind]),
           col=cols[pl])
  }
}
