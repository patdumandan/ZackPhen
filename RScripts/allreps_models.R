library(cmdstanr)
set_cmdstan_path("/users/dumandan/cmdstan")
cmdstan_version()

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)
library(ggridges)

source("~/plant_functions.R")
source("~/arthropod_functions.R")

plant_mod=cmdstan_model("~/plant_phen6b.stan")
arth_mod=cmdstan_model("~/arth_phen2.stan")

setwd("/scratch/project_2017453")

#data####
dat_path="/projappl/project_2017453"
dat_name=paste(dat_path, '/plant_datA','.csv', sep = '')
arth_name=paste(dat_path, '/arth_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)
arth_datA=read.csv(arth_name, header=T, sep=',',  stringsAsFactors = F)

#make sliding windows
dryas_slide=make_sliding_windows(dat=plant_datA, "Dryas")
silene_slide=make_sliding_windows(dat=plant_datA, "Silene")
papaver_slide=make_sliding_windows(dat=plant_datA, "Papaver")
cassiope_slide=make_sliding_windows(dat=plant_datA, "Cassiope")
salix_slide=make_sliding_windows(dat=plant_datA, "Salix")
saxifraga_slide=make_sliding_windows(dat=plant_datA, "Saxifraga")

muscidae_slide=make_sliding_windows(dat=arth_datA, "Muscidae")
ichneumonidae_slide=make_sliding_windows(dat=arth_datA, "Ichneumonidae")
chironomidae_slide=make_sliding_windows(dat=arth_datA, "Chironomidae")
sciaridae_slide=make_sliding_windows(dat=arth_datA, "Sciaridae")
phoridae_slide=make_sliding_windows(dat=arth_datA, "Phoridae")
acari_slide=make_sliding_windows(dat=arth_datA, "Acari")
collembola_slide=make_sliding_windows(dat=arth_datA, "Collembola")
lycosidae_slide=make_sliding_windows(dat=arth_datA, "Lycosidae")
linyphiidae_slide=make_sliding_windows(dat=arth_datA, "Linyphiidae")
nymphalidae_slide=make_sliding_windows(dat=arth_datA, "Nymphalidae")
coccoidea_slide=make_sliding_windows(dat=arth_datA, "Coccoidea")

#fit models
dryas_slide_models=lapply(dryas_slide, fit_plant_mod,model=plant_mod)
silene_slide_models=lapply(silene_slide, fit_plant_mod,model=plant_mod)
papaver_slide_models=lapply(papaver_slide, fit_plant_mod,model=plant_mod)
cassiope_slide_models=lapply(cassiope_slide, fit_plant_mod,model=plant_mod)
salix_slide_models=lapply(salix_slide, fit_plant_mod,model=plant_mod)
saxifraga_slide_models=lapply(saxifraga_slide, fit_plant_mod,model=plant_mod)

muscidae_slide_models=lapply(muscidae_slide, fit_arth_mod,model=arth_mod)
ichneumonidae_slide_models=lapply(ichneumonidae_slide, fit_arth_mod,model=arth_mod)
chironomidae_slide_models=lapply(chironomidae_slide, fit_arth_mod,model=arth_mod)
sciaridae_slide_models=lapply(sciaridae_slide, fit_arth_mod,model=arth_mod)
phoridae_slide_models=lapply(phoridae_slide, fit_arth_mod,model=arth_mod)
acari_slide_models=lapply(acari_slide, fit_arth_mod,model=arth_mod)
collembola_slide_models=lapply(collembola_slide, fit_arth_mod,model=arth_mod)
lycosidae_slide_models=lapply(lycosidae_slide, fit_arth_mod,model=arth_mod)
linyphiidae_slide_models=lapply(linyphiidae_slide, fit_arth_mod,model=arth_mod)
nymphalidae_slide_models=lapply(nymphalidae_slide, fit_arth_mod,model=arth_mod)
coccoidea_slide_models=lapply(coccoidea_slide, fit_arth_mod,model=arth_mod)
