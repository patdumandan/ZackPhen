make_sliding_windows=function(dat, species_name) {

  dat=dat%>%filter(species== species_name)

  step_size=1
  year_var=dat$year
  nyrs=length(unique(year_var))
  years_all=sort(unique(year_var))

  min_win=5
  max_win=length(years_all)

  slide_data <- list()

  for (nyr in seq(min_win, max_win, by = step_size)) {
    for (start in seq(1, length(years_all) - nyr + 1)) {

      yrs <- years_all[start:(start + nyr - 1)]

      slide_data [[length(slide_data) + 1]]=dat%>%filter(year %in% yrs)
    }}
  return(slide_data)
}

dryas_slide=make_sliding_windows(dat=plant_datA, "Dryas")

fit_plant_mod=function(data, model,
                        out_dir = "spp_models",
                        csv_dir = "spp_cmdstan_csv") {

  species_name=unique(data$species)
  sp_df= data%>%filter(species == species_name)
  yrs= data$year
  start_yr=min(yrs)
  end_yr=max(yrs)
  win_length=length(unique(yrs))
  window_id=paste0(species_name, "_", start_yr, "_", end_yr, "_", win_length)

  plant_data <- make_plant_data(sp_df, data)

  # Ensure output directories exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  species_csv_dir <- file.path(csv_dir, species_name)
  if (!dir.exists(species_csv_dir)) dir.create(species_csv_dir, recursive = TRUE)

  fit <- model$sample(
    data = plant_data,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 200,
    iter_warmup = 50,
    output_dir = species_csv_dir,
    save_warmup = TRUE)

  saveRDS(list(fit= fit,
               meta= list(
                 species= species_name,
                 start_yr= start_yr,
                 end_yr= end_yr,
                 tsl= win_length,
                 window_id = window_id)),
    file = file.path(out_dir, paste0("phenology_", window_id, ".rds")))
}

dryas_slide_models=lapply(dryas_slide, fit_plant_mod,model=plant_mod)

extract_beta_mu_win <- function(model_file) {

  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
  }

  obj <- readRDS(model_file)

  mod_fit=obj$fit
  metadata=obj$meta


  beta_mu_draws <- mod_fit$draws(variables = "beta_mu",format = "df")

  beta_mu_draws=beta_mu_draws%>%
  mutate(species   = metadata$species,
         start_yr  = metadata$start_yr,
         end_yr    = metadata$end_yr,
         tsl       = metadata$tsl,
         window_id = metadata$window_id)

  return(beta_mu_draws)
}

mod_list=list.files(path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\spp_models",
                    full.names=T)

dryas_betamu <- lapply(mod_list,extract_beta_mu_win)
