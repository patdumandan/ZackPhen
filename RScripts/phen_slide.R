require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)
require(lubridate)
require(date)
library(broom)
require(posterior)

source("RScripts/plant_functions.R")

dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
dat_name=paste(dat_path, '\\plant_datA','.csv', sep = '')

plant_datA=read.csv(dat_name, header=T, sep=',',  stringsAsFactors = F)
DOY_mean <- mean(plant_datA$DOY)
DOY_sd   <- sd(plant_datA$DOY)

#Dryas####
dry_datA=plant_datA%>%filter(species=="Dryas")

#extract peak DOY
dry_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Dryas.rds")

drydraws_doy=dry_fit$draws(variable="mu", format="df")

dryyears <- sort(unique(dry_datA$year))

drypeak_df <- as_draws_df(drydraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = dryyears[year_index])

drysummary_peak <- drypeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(dry_datA$year))

years_all <- sort(unique(drypeak_df$year))

dryp_slide_draws <- list()

for (nyr in seq(5, 24, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    dryp_slide_draws[[length(dryp_slide_draws) + 1]] <-
      drypeak_df %>% filter(year %in% yrs)
  }}

dry_slopes_draws <- lapply(dryp_slide_draws, fit_slopes_per_window)

dryas_peak_slope_ci <- map2_dfr(dry_slopes_draws,dryp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

dryp=ggplot(dryas_peak_slope_ci,
       aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Dryas",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Cassiope####
cas_datA=plant_datA%>%filter(species=="Cassiope")

#extract peak DOY
cas_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Cassiope.rds")

casdraws_doy=cas_fit$draws(variable="mu", format="df")

casyears <- sort(unique(cas_datA$year))

caspeak_df <- as_draws_df(casdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = casyears[year_index])

cassummary_peak <- caspeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(cas_datA$year))

years_all <- sort(unique(caspeak_df$year))

casp_slide_draws <- list()

for (nyr in seq(5, 24, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    casp_slide_draws[[length(casp_slide_draws) + 1]] <-
      caspeak_df %>% filter(year %in% yrs)
  }}

cas_slopes_draws <- lapply(casp_slide_draws, fit_slopes_per_window)

cassiope_peak_slope_ci <- map2_dfr(cas_slopes_draws,casp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

casp=ggplot(cassiope_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Cassiope",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))
#Silene####
sil_datA=plant_datA%>%filter(species=="Silene")

#extract peak DOY
sil_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Silene.rds")

sildraws_doy=sil_fit$draws(variable="mu", format="df")

silyears <- sort(unique(sil_datA$year))

silpeak_df <- as_draws_df(sildraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = silyears[year_index])

silsummary_peak <- silpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(sil_datA$year))

years_all <- sort(unique(silpeak_df$year))

silp_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    silp_slide_draws[[length(silp_slide_draws) + 1]] <-
      silpeak_df %>% filter(year %in% yrs)
  }}

sil_slopes_draws <- lapply(silp_slide_draws, fit_slopes_per_window)

silene_peak_slope_ci <- map2_dfr(sil_slopes_draws,silp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

silp=ggplot(silene_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Silene",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Salix####
sal_datA=plant_datA%>%filter(species=="Salix")

#extract peak DOY
sal_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Salix.rds")

saldraws_doy=sal_fit$draws(variable="mu", format="df")

salyears <- sort(unique(sal_datA$year))

salpeak_df <- as_draws_df(saldraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = salyears[year_index])

salsummary_peak <- salpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(sal_datA$year))

years_all <- sort(unique(salpeak_df$year))

salp_slide_draws <- list()

for (nyr in seq(5, 14, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    salp_slide_draws[[length(salp_slide_draws) + 1]] <-
      salpeak_df %>% filter(year %in% yrs)
  }}

sal_slopes_draws <- lapply(salp_slide_draws, fit_slopes_per_window)

salix_peak_slope_ci <- map2_dfr(sal_slopes_draws,salp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

salp=ggplot(salix_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Salix",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Saxifraga####
sax_datA=plant_datA%>%filter(species=="Saxifraga")

#extract peak DOY
sax_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Saxifraga.rds")

saxdraws_doy=sax_fit$draws(variable="mu", format="df")

saxyears <- sort(unique(sax_datA$year))

saxpeak_df <- as_draws_df(saxdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = saxyears[year_index])

saxsummary_peak <- saxpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(sax_datA$year))

years_all <- sort(unique(saxpeak_df$year))

saxp_slide_draws <- list()

for (nyr in seq(5, 11, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    saxp_slide_draws[[length(saxp_slide_draws) + 1]] <-
      saxpeak_df %>% filter(year %in% yrs)
  }}

sax_slopes_draws <- lapply(saxp_slide_draws, fit_slopes_per_window)

saxifraga_peak_slope_ci <- map2_dfr(sax_slopes_draws,saxp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

saxp=ggplot(saxifraga_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Saxifraga",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Papaver####
pap_datA=plant_datA%>%filter(species=="Papaver")

#extract peak DOY
pap_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Papaver.rds")

papdraws_doy=pap_fit$draws(variable="mu", format="df")

papyears <- sort(unique(pap_datA$year))

pappeak_df <- as_draws_df(papdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = papyears[year_index])

papsummary_peak <- pappeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(pap_datA$year))

years_all <- sort(unique(pappeak_df$year))

papp_slide_draws <- list()

for (nyr in seq(5, 21, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    papp_slide_draws[[length(papp_slide_draws) + 1]] <-
      pappeak_df %>% filter(year %in% yrs)
  }}

pap_slopes_draws <- lapply(papp_slide_draws, fit_slopes_per_window)

papaver_peak_slope_ci <- map2_dfr(pap_slopes_draws,papp_slide_draws, summarize_slopes)%>%
  mutate(TSL = end_yr - start_yr + 1,
         CI_width = slope_upr - slope_lwr)

papp=ggplot(papaver_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Papaver",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#plot####
ggarrange(dryp,silp)
ggarrange(casp,papp)
ggarrange(salp,saxp)
