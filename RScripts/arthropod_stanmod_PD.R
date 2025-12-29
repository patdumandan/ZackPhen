library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)

#data####
dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
arth_name=paste(dat_path, '\\arth_datA','.csv', sep = '')

arth_datA=read.csv(arth_name, header=T, sep=',',  stringsAsFactors = F)
#or: arth_datA=read.csv("https://raw.githubusercontent.com/patdumandan/ZackPhen/refs/heads/main/data/arth_datA.csv")

species_list= unique(arth_datA$HoyeTaxon)

#model####
arth_mod=cmdstan_model("RScripts/arth_phen2.stan")

arth_fits= lapply(species_list, fit_arth_model,data= arth_datA, model= arth_mod)

names(arth_fits)= species_list

#diagnostics####
diagnost_arth=lapply(species_list, run_model_diagnostics, data=arth_datA)
