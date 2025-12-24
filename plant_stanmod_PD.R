library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)

source()

#data
dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
dat_name=paste(dat_path, '\\plant_datA','.csv', sep = '')

#dat_path = paste0(getwd(), "/data")
#dat_name = paste(dat_path, '/plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

#model
plant_mod=cmdstan_model("plant_phen6b.stan")

species_list <- unique(plant_datA$species)

fits <- lapply(species_list, fit_species_model,data = plant_datA, model = plant_mod)

names(fits) <- species_list

#diagnostics

species_list <- unique(plant_datA$species)

pars_to_plot <- c("alpha", "mu", "beta_mu", "width_bar")

for (sp in species_list) {

  fit <- readRDS(
    file.path("models", paste0("phenology_", sp, ".rds"))
  )

  run_model_diagnostics(
    fit = fit,
    species_name = sp,
    data = plant_datA,
    pars = pars_to_plot
  )
}
