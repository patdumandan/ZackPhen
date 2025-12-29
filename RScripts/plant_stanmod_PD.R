library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)

source("plant_functions.R")

#data####
dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
dat_name=paste(dat_path, '\\plant_datA','.csv', sep = '')

#dat_path = paste0(getwd(), "/data")
#dat_name = paste(dat_path, '/plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)

species_list= unique(plant_datA$species)

#model####
plant_mod=cmdstan_model("RScripts/plant_phen6b.stan")

fits= lapply(species_list, fit_species_model,data= plant_datA, model= plant_mod)

names(fits)= species_list

#diagnostics####
diagnost_output=lapply(species_list, run_model_diagnostics, data=plant_datA, pars=pars_to_plot)

#param estimates####
params_output=lapply(species_list, plot_params, data=plant_datA)

#predictions####
preds_output=lapply(species_list, plot_preds, data=plant_datA)
