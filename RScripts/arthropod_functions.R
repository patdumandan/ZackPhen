make_sliding_windows=function(dat, species_name) {

  dat=dat%>%filter(HoyeTaxon== species_name)

  step_size=1
  year_var=dat$Year
  nyrs=length(unique(year_var))
  years_all=sort(unique(year_var))

  min_win=5
  max_win=length(years_all)

  slide_data <- list()

  for (nyr in seq(min_win, max_win, by = step_size)) {
    for (start in seq(1, length(years_all) - nyr + 1)) {

      yrs <- years_all[start:(start + nyr - 1)]

      slide_data [[length(slide_data) + 1]] <-
        dat%>%filter(Year %in% yrs)
    }}
  return(slide_data)
}
make_arth_data <- function(df, global_df) {
  years <- unique(df$Year)
  list(N= nrow(df),
       y=df$TotalCatch1,
       obs_days=df$TrapDays,
       Nyr = length(unique(df$Year)),
       year_id= as.integer(factor(df$Year)),
       DOYs= df$DOYs,
       year= df$Year,
       year_vec = (years-mean(years))/sd(years),
       plot_id = as.integer(factor(df$plot_id)),
       Nplots  = length(unique(df$plot_id)),
       DOY_sd  = sd(global_df$DOY),
       DOY_mean = mean(global_df$DOY))
}

#use for sliding window model fitting
fit_arth_mod=function(data, model,
                      out_dir = "arth_models") {

  species_name=unique(data$HoyeTaxon)
  sp_df= data%>%filter(HoyeTaxon == species_name)
  yrs= data$Year
  start_yr=min(yrs)
  end_yr=max(yrs)
  win_length=length(unique(yrs))
  window_id=paste0(species_name, "_", start_yr, "_", end_yr, "_", win_length)

  arth_data= make_arth_data(sp_df, data)

  #create dir that contains both csv and model files
  fit_dir= file.path(out_dir, paste0("phenology_", window_id))

  species_csv_dir= file.path(fit_dir, "cmdstan_csv")

  #make sure output dir exists
  dir.create(species_csv_dir, recursive = TRUE)

  fit= model$sample(
    data= arth_data,
    seed= 123,
    chains= 4,
    parallel_chains= 4,
    iter_sampling= 2000,
    iter_warmup= 500,
    output_dir= species_csv_dir,
    save_warmup= FALSE) #to reduce size, not save warmup

  fit_draws=fit$draws(format= "draws_df") #save only draws to reduce size

  saveRDS(list(fit_draws= fit_draws,
               meta= list(
                 species_name= species_name,
                 start_yr= start_yr,
                 end_yr= end_yr,
                 tsl= win_length,
                 window_id = window_id)),
          file = file.path(out_dir, paste0("phenology_", window_id, ".rds")))
}

fit_arth_model <- function(
    species_name,
    data,
    model,
    out_dir = "models",
    csv_dir = "cmdstan_csv") {

  message("Fit model for: ", species_name)

  sp_df <- data %>% filter(HoyeTaxon == species_name)

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

chain_check <- function(fit) {
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

run_model_diagnostics <- function(species_name, data,
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

  sp_df <- data %>% dplyr::filter(HoyeTaxon == species_name)

  diag <- extract_diagnostics(fit)

  #Rhat / Neff
  pdf(file.path(out_dir, paste0(species_name, "_mcmc_diagnostics.pdf")),width = 12, height = 4)
  print(plot_mcmc_diagnostics(diag))
  dev.off()

  # chains
  pdf(file.path(out_dir, paste0(species_name, "_caterpillar_plots.pdf")),width = 8, height = 6)
  print(chain_check(fit))
  dev.off()

  # ppc
  pdf(file.path(out_dir, paste0(species_name, "_ppc_plots.pdf")),width = 10, height = 8)
  print(plot_ppc(fit, y_obs = sp_df$TotalCatch1))
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

  sp_df <- data %>% dplyr::filter(HoyeTaxon == species_name)

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
  sp_df <- dplyr::filter(data, HoyeTaxon == species_name)

  mu         <- fit$draws("mu", format="draws_matrix")
  width      <- fit$draws("width", format="draws_matrix")
  alpha_year <- fit$draws("alpha_year", format="draws_matrix")
  u_plot     <- fit$draws("u_plot", format="draws_matrix")
  u_plot_mu  <- fit$draws("u_plot_mu", format="draws_matrix")

  # Setup
  log_obs_days <- log(mean(sp_df$TrapDays, na.rm = TRUE))

  DOY <- seq(-2, 2, length.out=100)

  # PDF
  pdf(file.path(out_dir, paste0(species_name, "_preds.pdf")))
  par(mfrow=c(4, 4))

  plot_ids <- sort(unique(sp_df$plot_id))
  year_ids <- sort(unique(sp_df$Year))  # actual year numbers
  P <- length(plot_ids)
  cols <- rainbow(P)

  for (y_idx in seq_along(year_ids)) {
    y <- year_ids[y_idx]  # actual year

    for (pl in seq_along(plot_ids)) {
      p_id <- plot_ids[pl]
      pl_idx <- which(plot_ids == p_id)  # correct column index

      eta_hat <- numeric(length(DOY))

      for (nd in seq_along(DOY)) {
        eta_post <- exp(alpha_year[, y_idx] + u_plot[, pl_idx] -
                          ((DOY[nd] - (mu[, y_idx] + u_plot_mu[, pl_idx]))^2 / (width[, y_idx]^2)) +
                          log_obs_days)

        eta_hat[nd] <- mean(eta_post)
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
      points(sp_df$DOYs[ind], sp_df$TotalCatch1[ind],
             col=cols[pl], pch=16)
    }
  }

  dev.off()
}

