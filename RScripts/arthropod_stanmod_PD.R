library(cmdstanr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(bayesplot)
library(posterior)
library(ggridges)

source("RScripts/arthropod_functions.R")

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

#param estimates####
params_arth=lapply(species_list, plot_params, data=arth_datA)

#predictions####
preds_arth=lapply(species_list, plot_arth_preds, data=arth_datA)

#beta_mu####
beta_mu_arth <- purrr::map_dfr(species_list, extract_beta_mu, data = arth_datA)

consumer_order <- c(
  "Chironomidae", "Nymphalidae", "Coccoidea", "Linyphiidae", "Lycosidae",
  "Ichneumonidae", "Collembola", "Phoridae", "Acari","Sciaridae", "Muscidae")

arth_beta_mu_summary <- beta_mu_arth %>%
  group_by(species) %>%
  summarise(
    mean  = mean(beta_mu),
    lwr   = quantile(beta_mu, 0.025),
    upr   = quantile(beta_mu, 0.975))%>%
  mutate(consumer = factor(species, levels = consumer_order))

m3=ggplot(arth_beta_mu_summary,
          aes(x = mean, y = consumer)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = expression(beta[mu] ~ "(days per year)"),
    y = "Species",
    title = "Arthropods",
    subtitle = "mean and 95 % CI") +
  theme_classic(base_size = 13)
