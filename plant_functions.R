make_stan_data <- function(df, global_df) {
  years <- unique(df$year)
  list(N= nrow(df),
       tot_F= df$tot_F,
       tot_NF= df$tot_NF,
       Nyr= length(years),
       year_id= as.integer(factor(df$year)),
       DOYs= df$DOYs,
       year= df$year,
       year_vec = (years-mean(years))/sd(years),
       Nplots  = length(unique(df$plot_id)),
       plot_id = df$plot_id,
       DOY_sd  = sd(global_df$DOY),
       DOY_mean = mean(global_df$DOY))
}

fit_species_model <- function(species_name, data, model, out_dir = "models") {

  message("Fit model for: ", species_name)

  sp_df <- data %>% filter(species == species_name)

  stan_data <- make_stan_data(sp_df, data)

  fit <- model$sample(
    data = stan_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 2000,
    iter_warmup = 500,
    output_dir = "cmdstan_csv")

  if (!dir.exists(out_dir)) dir.create(out_dir)

  saveRDS(fit,file = file.path(out_dir, paste0("phenology_", species_name, ".rds")))

  return(fit)
}

extract_diagnostics <- function(fit) {

  draws_array <- fit$draws() %>% posterior::as_draws_array()

  ess_vals  <- posterior::ess_bulk(draws_array)
  n_draws  <- posterior::ndraws(draws_array)
  neff_ratio_vals <- ess_vals / n_draws

  list(
    draws_array = draws_array,
    rhat = fit$summary()$rhat,
    neff_ratio = neff_ratio_vals
  )
}

plot_mcmc_diagnostics <- function(diagnost_list) {

  p1 <- mcmc_rhat_hist(diagnost_list$rhat)
  p2 <- mcmc_rhat(diagnost_list$rhat)
  p3 <- mcmc_neff_hist(diagnost_list$neff_ratio)

  ggarrange(p1, p2, p3, ncol = 3)
}

pars=c("mu")

chain_check <- function(fit, pars) {
  draws <-fit$draws()%>%posterior::as_draws_array()

  varnames <- variables(draws)
  mu_params <- grep("^mu\\[", varnames, value = TRUE)
  width_params <- grep("^width\\[", varnames, value = TRUE)

  mcmc_trace(draws, pars = mu_params)
  mcmc_trace(draws, pars = width_params)
  mcmc_trace(draws, pars = c("sigma_mu","sigma_width", "sigma_year"))
}

plot_ppc <- function(fit, y_obs, yrep_var = "y_pred", ndraws = 100) {

  color_scheme_set("blue")

  yrep <- fit$draws(format = "draws_matrix",variables = yrep_var)

  yrep <- yrep[(nrow(yrep) - ndraws + 1):nrow(yrep), ]

  ggarrange(
    ppc_dens_overlay(y = y_obs, yrep = yrep),
    ppc_ecdf_overlay(y = y_obs, yrep = yrep),
    ppc_pit_ecdf(y = y_obs, yrep = yrep, prob = 0.99, plot_diff = FALSE),
    ppc_pit_ecdf(y = y_obs, yrep = yrep, prob = 0.99, plot_diff = TRUE),
    nrow = 2, ncol = 2
  )
}

run_model_diagnostics <- function(
    fit,
    species_name,
    data,
    pars,
    out_dir = "diagnostics"
) {

  message("Running diagnostics for: ", species_name)

  # Ensure central diagnostics directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  sp_df <- data %>% filter(species == species_name)

  diag <- extract_diagnostics(fit)

  # ---- Rhat / Neff ----
  pdf(
    file.path(out_dir, paste0(species_name, "_mcmc_diagnostics.pdf")),
    width = 12, height = 4
  )
  print(plot_mcmc_diagnostics(diag))
  dev.off()

  # ---- Caterpillar plots ----
  pdf(
    file.path(out_dir, paste0(species_name, "_caterpillar_plots.pdf")),
    width = 8, height = 6
  )
  print(chain_check(fit, pars))
  dev.off()

  # ---- Posterior predictive checks ----
  pdf(
    file.path(out_dir, paste0(species_name, "_ppc_plots.pdf")),
    width = 10, height = 8
  )
  print(plot_ppc(fit, y_obs = sp_df$tot_F))
  dev.off()

  invisible(TRUE)
}


plot_mu_est <- function(fit) {

  mu_draws= fit$draws("mu")
  mcmc_areas(mu_draws, regex = "mu")
}

plot_width_est <- function(fit) {

  width_draws= fit$draws("width")
  mcmc_areas(width_draws, regex = "width")
}



###moving window###
rolling_mod=function(split) {

  analysis_set=analysis(split)
  Plot=analysis_set[, "Plot"]

  fit_model= lme4::glmer(analysis_set[, "DOY"]~analysis_set[, "Year"]+
                           (1|Plot),data=analysis_set)
  return(summary(fit_model))

}
rolling_plot=function(split, slope, intercept) {

  analysis_set=analysis(split)
  slope=slope
  intercept=intercept

  plot_model=ggplot2::ggplot(data=analysis_set,
                             aes(x=as.integer(analysis_set[, "Year"]), y=analysis_set[, "DOY"]))+
    geom_point()+
    geom_abline(slope=slope, intercept=intercept)+theme_classic()+
    ylab("DOY")+ xlab("Year")+
    stat_smooth(method="lm")

  print(plot_model)
}

get_year=function(split) {

  analysis_set=analysis(split)

  start_yr=as.numeric(analysis_set[,"Year"][1])

}

get_dat=function(split) {

  analysis_set=analysis(split)

  yr=as.vector(analysis_set[,"Year"])
  doys=as.vector(analysis_set[,"DOY"])

  return(cbind(yr, doys))
}

generate_fitted_curves <- function(n_years, alpha_mean, beta_doy, beta_doy_sq,
                                   doy_std, doy_sq_std, doy_seq, year_levels) {
  lapply(1:n_years, function(i) {
    eta <- alpha_mean + beta_doy[i] * doy_std + beta_doy_sq[i] * doy_sq_std
    p <- plogis(eta)

    data.frame(
      DOY = doy_seq,
      prob = p,
      year = rep(year_levels[i], length(doy_seq))
    )
  })
}

get_dat
