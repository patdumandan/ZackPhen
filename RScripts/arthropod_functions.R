make_arth_data <- function(df, global_df) {
  years <- unique(df$year)
  list(N= nrow(df),
       y=df$TotalCatch1,
       obs_days=df$TrapDays,
       Nyr = length(unique(df$Year)),
       year_id= as.integer(factor(df$year)),
       DOYs= df$DOYs,
       year= df$Year,
       year_vec = (year-mean(year))/sd(year),
       Nplots  = length(unique(df$plot_id)),
       plot_id = df$plot_id,
       DOY_sd  = sd(global_df$DOY),
       DOY_mean = mean(global_df$DOY))
}

fit_arth_model <- function(
    species_name,
    data,
    model,
    out_dir = "models",
    csv_dir = "cmdstan_csv") {

  message("Fit model for: ", species_name)

  sp_df <- data %>% filter(species == species_name)

  arth_data <- make_arth_data(sp_df, data)

  # Ensure output directories exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  species_csv_dir <- file.path(csv_dir, species_name)
  if (!dir.exists(species_csv_dir)) dir.create(species_csv_dir, recursive = TRUE)

  fit <- model$sample(
    data = arth_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 2000,
    iter_warmup = 500,
    output_dir = species_csv_dir,
    save_warmup = TRUE)

  saveRDS(fit, file = file.path(out_dir, paste0("phenology_", species_name, ".rds")))

  return(fit)
}

extract_diagnostics <- function(fit) {

  draws_array <- fit$draws() %>% posterior::as_draws_array()

  ess_vals  <- posterior::ess_bulk(draws_array)
  n_draws  <- posterior::ndraws(draws_array)
  neff_ratio_vals <- ess_vals / n_draws

  list(draws_array = draws_array,
       rhat = fit$summary()$rhat,
       neff_ratio = neff_ratio_vals)
}

plot_mcmc_diagnostics <- function(diagnost_list) {

  p1 <- mcmc_rhat_hist(diagnost_list$rhat)
  p2 <- mcmc_rhat(diagnost_list$rhat)
  p3 <- mcmc_neff_hist(diagnost_list$neff_ratio)

  ggarrange(p1, p2, p3, ncol = 3)
}

chain_check <- function(fit, pars) {
  draws <-fit$draws()%>%posterior::as_draws_array()

  # varnames <- variables(draws)
  # mu_params <- grep("^mu\\[", varnames, value = TRUE)
  # width_params <- grep("^width\\[", varnames, value = TRUE)

  #mcmc_trace(draws, pars = mu_params)
  #mcmc_trace(draws, pars = width_params)
  mcmc_trace(draws, pars = c("sigma_mu","sigma_width", "sigma_year", "phi"))
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
    nrow = 2, ncol = 2)
}

run_model_diagnostics <- function(species_name, data, pars,
                                  model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\",
                                  out_dir = "diagnostics") {

  message("Running diagnostics for: ", species_name)

  #open file
  model_file <- file.path(model_dir,paste0("phenology_", species_name, ".rds"))

  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }

  #read model
  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  sp_df <- data %>% dplyr::filter(species == species_name)

  diag <- extract_diagnostics(fit)

  #Rhat / Neff
  pdf(file.path(out_dir, paste0(species_name, "_mcmc_diagnostics.pdf")),width = 12, height = 4)
  print(plot_mcmc_diagnostics(diag))
  dev.off()

  # chains
  pdf(file.path(out_dir, paste0(species_name, "_caterpillar_plots.pdf")),width = 8, height = 6)
  print(chain_check(fit, pars))
  dev.off()

  # ppc
  pdf(file.path(out_dir, paste0(species_name, "_ppc_plots.pdf")),width = 10, height = 8)
  print(plot_ppc(fit, y_obs = sp_df$tot_F))
  dev.off()

  invisible(TRUE)
}

plot_params= function(species_name, data,
                      model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\",
                      out_dir = "predictions") {

  message("Plotting parameters of: ", species_name)

  #open file
  model_file <- file.path(model_dir,paste0("phenology_", species_name, ".rds"))

  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }

  #read model
  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  sp_df <- data %>% dplyr::filter(species == species_name)

  mu_draws= fit$draws("mu")
  width_draws= fit$draws("width")

  #peak and width
  pdf(file.path(out_dir, paste0(species_name, "_params.pdf")),width = 12, height = 4)
  m1=mcmc_areas(mu_draws, regex = "mu")+labs(title=paste0("Peak of ", species_name))
  m2=mcmc_areas(width_draws, regex = "width")+ labs(title=paste0("Width of ", species_name))

  m12=ggarrange(m1, m2, ncol = 2)
  print(m12)

  dev.off()
}

plot_arth_preds <- function(species_name, data,
                            model_dir = "C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\",
                            out_dir = "predictions") {

  message("Plotting predictions of: ", species_name)

  # Open model
  model_file <- file.path(model_dir, paste0("phenology_", species_name, ".rds"))
  if (!file.exists(model_file)) stop("Model file not found: ", model_file)
  fit <- readRDS(model_file)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Subset data
  sp_df <- dplyr::filter(data, species == species_name)

  mu         <- fit$draws("mu", format="draws_matrix")
  width      <- fit$draws("width", format="draws_matrix")
  alpha_year <- fit$draws("alpha_year", format="draws_matrix")
  u_plot     <- fit$draws("u_plot", format="draws_matrix")
  u_plot_mu  <- fit$draws("u_plot_mu", format="draws_matrix")
  obs_days   <- data$TrapDays

  # Setup
  DOY <- seq(-2, 2, length.out=100)

  plot_ids = sort(unique(sp_df$plot_id))
  year_ids = sort(as.integer(factor(sp_df$year)))
  P <- length(plot_ids)
  cols <- rainbow(P)

  # PDF
  pdf(file.path(out_dir, paste0(species_name, "_preds.pdf")))
  par(mfrow=c(4, 4))

  plot_ids <- sort(unique(sp_df$plot_id))
  year_ids <- sort(unique(sp_df$year))  # actual year numbers
  P <- length(plot_ids)
  cols <- rainbow(P)

  for (y_idx in seq_along(year_ids)) {
    y <- year_ids[y_idx]  # actual year

    for (pl in seq_along(plot_ids)) {
      p_id <- plot_ids[pl]
      pl_idx <- which(plot_ids == p_id)  # correct column index

      eta_hat <- numeric(length(DOY))

      for (nd in seq_along(DOY)) {
        eta_post <- alpha_year[, y_idx] + u_plot[, pl_idx] -
          ((DOY[nd] - (mu[, y_idx] + u_plot_mu[, pl_idx]))^2 / (width[, y_idx]^2)) +
          log(obs_days[nd])

        eta_hat[nd] <- mean(exp(eta_post))
      }

      if (pl == 1) {
        plot(DOY, eta_hat, type="l", col=cols[pl],
             main=paste(species_name, "Year", y),
             xlab="Scaled DOY", ylab="Abundance")
      } else {
        lines(DOY, eta_hat, col=cols[pl])
      }

      # raw points
      ind <- sp_df$plot_id == p_id & sp_df$Year == y
      points(sp_df$y[ind],
             col=cols[pl], pch=16)
    }
  }

  dev.off()
}

