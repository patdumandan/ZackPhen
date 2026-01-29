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

dryas_slope_draws <- map2_dfr(dry_slopes_draws,dryp_slide_draws, summarize_slopes)

dryas_slope_df <- dryas_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

dryp=ggplot(dryas_slope_df,
       aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y =NULL,
       title = "Dryas")+
  #     subtitle="Trend uncertainty across different time windows") +
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

cassiope_slope_draws <- map2_dfr(cas_slopes_draws,casp_slide_draws, summarize_slopes)
 # mutate(TSL = end_yr - start_yr + 1,
  #       CI_width = slope_upr - slope_lwr)

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

silene_slope_draws <- map2_dfr(sil_slopes_draws,silp_slide_draws, summarize_slopes)

silene_slope_df <- silene_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

silp=ggplot(silene_slope_df,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y = NULL,
       title = "Silene")+
#       subtitle="Trend uncertainty across different time windows") +
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

salix_slope_draws <- map2_dfr(sal_slopes_draws,salp_slide_draws, summarize_slopes)

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

saxifraga_slope_draws <- map2_dfr(sax_slopes_draws,saxp_slide_draws, summarize_slopes)

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

papaver_slope_draws <- map2_dfr(pap_slopes_draws,papp_slide_draws, summarize_slopes)

papaver_slope_df <- papaver_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

papp=ggplot(papaver_slope_df,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y = NULL,
       title = "Papaver")+
   #    subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#PLANT PLOTS####
ggarrange(dryp,silp)
ggarrange(casp,papp)
ggarrange(salp,saxp)

dat_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
arth_name=paste(dat_path, '\\arth_datA','.csv', sep = '')

arth_datA=read.csv(arth_name, header=T, sep=',',  stringsAsFactors = F)
#or: arth_datA=read.csv("https://raw.githubusercontent.com/patdumandan/ZackPhen/refs/heads/main/data/arth_datA.csv")

#Muscidae####
mus_datA=arth_datA%>%filter(HoyeTaxon=="Muscidae")

#extract peak DOY
mus_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Muscidae.rds")

musdraws_doy=mus_fit$draws(variable="mu", format="df")

musyears <- sort(unique(mus_datA$Year))

muspeak_df <- as_draws_df(musdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = musyears[year_index])

mussummary_peak <- muspeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(mus_datA$year))

years_all <- sort(unique(muspeak_df$year))

musp_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    musp_slide_draws[[length(musp_slide_draws) + 1]] <-
      muspeak_df %>% filter(year %in% yrs)
  }}

mus_slopes_draws <- lapply(musp_slide_draws, fit_slopes_per_window)

muscidae_slope_draws <- map2_dfr(mus_slopes_draws,musp_slide_draws, summarize_slopes)

muscidae_slope_df <- muscidae_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

musp=ggplot(muscidae_slope_df,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y = NULL,
       title = "Muscidae")+
  #     subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Sciaridae####
sci_datA=arth_datA%>%filter(HoyeTaxon=="Sciaridae")

#extract peak DOY
sci_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Sciaridae.rds")

scidraws_doy=sci_fit$draws(variable="mu", format="df")

sciyears <- sort(unique(sci_datA$Year))

scipeak_df <- as_draws_df(scidraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = sciyears[year_index])

scisummary_peak <- scipeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(sci_datA$year))

years_all <- sort(unique(scipeak_df$year))

scip_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    scip_slide_draws[[length(scip_slide_draws) + 1]] <-
      scipeak_df %>% filter(year %in% yrs)
  }}

sci_slopes_draws <- lapply(scip_slide_draws, fit_slopes_per_window)

sciaridae_slope_draws <- map2_dfr(sci_slopes_draws,scip_slide_draws, summarize_slopes)
 # mutate(TSL = end_yr - start_yr + 1,
  #       CI_width = slope_upr - slope_lwr,
   #      year = sciyears[start_yr])

scip=ggplot(sciaridae_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Sciaridae",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Lycosidae####
lyc_datA=arth_datA%>%filter(HoyeTaxon=="Lycosidae")

#extract peak DOY
lyc_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Lycosidae.rds")

lycdraws_doy=lyc_fit$draws(variable="mu", format="df")

lycyears <- sort(unique(lyc_datA$Year))

lycpeak_df <- as_draws_df(lycdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = lycyears[year_index])

lycsummary_peak <- lycpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(lyc_datA$year))

years_all <- sort(unique(lycpeak_df$year))

lycp_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    lycp_slide_draws[[length(lycp_slide_draws) + 1]] <-
      lycpeak_df %>% filter(year %in% yrs)
  }}

lyc_slopes_draws <- lapply(lycp_slide_draws, fit_slopes_per_window)

lycosidae_slope_draws<- map2_dfr(lyc_slopes_draws,lycp_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = lycyears[start_yr])

lycp=ggplot(lycosidae_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Lycosidae",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Phoridae####
pho_datA=arth_datA%>%filter(HoyeTaxon=="Phoridae")

#extract peak DOY
pho_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Phoridae.rds")

phodraws_doy=pho_fit$draws(variable="mu", format="df")

phoyears <- sort(unique(pho_datA$Year))

phopeak_df <- as_draws_df(phodraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = phoyears[year_index])

phosummary_peak <- phopeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(pho_datA$year))

years_all <- sort(unique(phopeak_df$year))

phop_slide_draws <- list()

for (nyr in seq(5, 21, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    phop_slide_draws[[length(phop_slide_draws) + 1]] <-
      phopeak_df %>% filter(year %in% yrs)
  }}

pho_slopes_draws <- lapply(phop_slide_draws, fit_slopes_per_window)

phoridae_slope_draws <- map2_dfr(pho_slopes_draws,phop_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = phoyears[start_yr])

phop=ggplot(phoridae_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Phoridae",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Collembola####
col_datA=arth_datA%>%filter(HoyeTaxon=="Collembola")

#extract peak DOY
col_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Collembola.rds")

coldraws_doy=col_fit$draws(variable="mu", format="df")

colyears <- sort(unique(col_datA$Year))

colpeak_df <- as_draws_df(coldraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = colyears[year_index])

colsummary_peak <- colpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(col_datA$year))

years_all <- sort(unique(colpeak_df$year))

colp_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    colp_slide_draws[[length(colp_slide_draws) + 1]] <-
      colpeak_df %>% filter(year %in% yrs)
  }}

col_slopes_draws <- lapply(colp_slide_draws, fit_slopes_per_window)

collembola_slope_draws <- map2_dfr(col_slopes_draws,colp_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = colyears[start_yr])

colp=ggplot(collembola_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Collembola",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Ichneumonidae####
ich_datA=arth_datA%>%filter(HoyeTaxon=="Ichneumonidae")

#extract peak DOY
ich_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Ichneumonidae.rds")

ichdraws_doy=ich_fit$draws(variable="mu", format="df")

ichyears <- sort(unique(ich_datA$Year))

ichpeak_df <- as_draws_df(ichdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = ichyears[year_index])

ichsummary_peak <- ichpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(ich_datA$year))

years_all <- sort(unique(ichpeak_df$year))

ichp_slide_draws <- list()

for (nyr in seq(5, 24, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    ichp_slide_draws[[length(ichp_slide_draws) + 1]] <-
      ichpeak_df %>% filter(year %in% yrs)
  }}

ich_slopes_draws <- lapply(ichp_slide_draws, fit_slopes_per_window)

ichneumonidae_slope_draws <- map2_dfr(ich_slopes_draws,ichp_slide_draws, summarize_slopes)

ichneumonidae_slope_df <- ichneumonidae_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

ichp=ggplot(ichneumonidae_slope_df,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y = NULL,
       title = "Ichneumonidae")+
#       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Linyphiidae####
lin_datA=arth_datA%>%filter(HoyeTaxon=="Linyphiidae")

#extract peak DOY
lin_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Linyphiidae.rds")

lindraws_doy=lin_fit$draws(variable="mu", format="df")

linyears <- sort(unique(lin_datA$Year))

linpeak_df <- as_draws_df(lindraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = linyears[year_index])

linsummary_peak <- linpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(lin_datA$year))

years_all <- sort(unique(linpeak_df$year))

linp_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    linp_slide_draws[[length(linp_slide_draws) + 1]] <-
      linpeak_df %>% filter(year %in% yrs)
  }}

lin_slopes_draws <- lapply(linp_slide_draws, fit_slopes_per_window)

linyphiidae_slope_draws <- map2_dfr(lin_slopes_draws,linp_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = linyears[start_yr])

linp=ggplot(linyphiidae_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Linyphiidae",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Coccoidea####
coc_datA=arth_datA%>%filter(HoyeTaxon=="Coccoidea")

#extract peak DOY
coc_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Coccoidea.rds")

cocdraws_doy=coc_fit$draws(variable="mu", format="df")

cocyears <- sort(unique(coc_datA$Year))

cocpeak_df <- as_draws_df(cocdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = cocyears[year_index])

cocsummary_peak <- cocpeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(coc_datA$year))

years_all <- sort(unique(cocpeak_df$year))

cocp_slide_draws <- list()

for (nyr in seq(5, 16, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    cocp_slide_draws[[length(cocp_slide_draws) + 1]] <-
      cocpeak_df %>% filter(year %in% yrs)
  }}

coc_slopes_draws <- lapply(cocp_slide_draws, fit_slopes_per_window)

coccoidea_slope_draws <- map2_dfr(coc_slopes_draws,cocp_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = cocyears[start_yr])

cocp=ggplot(coccoidea_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Coccoidea",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Nymphalidae####
nym_datA=arth_datA%>%filter(HoyeTaxon=="Nymphalidae")

#extract peak DOY
nym_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Nymphalidae.rds")

nymdraws_doy=nym_fit$draws(variable="mu", format="df")

nymyears <- sort(unique(nym_datA$Year))

nympeak_df <- as_draws_df(nymdraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = nymyears[year_index])

nymsummary_peak <- nympeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(nym_datA$year))

years_all <- sort(unique(nympeak_df$year))

nymp_slide_draws <- list()

for (nyr in seq(5, 9, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    nymp_slide_draws[[length(nymp_slide_draws) + 1]] <-
      nympeak_df %>% filter(year %in% yrs)
  }}

nym_slopes_draws <- lapply(nymp_slide_draws, fit_slopes_per_window)

nymphalidae_slope_draws<- map2_dfr(nym_slopes_draws,nymp_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = nymyears[start_yr])

nymp=ggplot(nymphalidae_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Nymphalidae",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#Chironomidae####
chi_datA=arth_datA%>%filter(HoyeTaxon=="Chironomidae")

#extract peak DOY
chi_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Chironomidae.rds")

chidraws_doy=chi_fit$draws(variable="mu", format="df")

chiyears <- sort(unique(chi_datA$Year))

chipeak_df <- as_draws_df(chidraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = chiyears[year_index])

chisummary_peak <- chipeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(chi_datA$year))

years_all <- sort(unique(chipeak_df$year))

chip_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    chip_slide_draws[[length(chip_slide_draws) + 1]] <-
      chipeak_df %>% filter(year %in% yrs)
  }}

chi_slopes_draws <- lapply(chip_slide_draws, fit_slopes_per_window)

chironomidae_slope_draws<- map2_dfr(chi_slopes_draws,chip_slide_draws, summarize_slopes)

chironomidae_slope_df <- chironomidae_slope_draws %>%
  group_by(start_yr, end_yr, TSL) %>%
  summarise(
    slope_median = median(slope),
    slope_mean   = mean(slope),
    slope_sd= sd(slope),
    CI_width = quantile(slope, 0.975) - quantile(slope, 0.025),
    .groups = "drop")

chip=ggplot(chironomidae_slope_df,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = NULL,y = NULL,
       title = "Chironomidae")+
#       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#####for window categories
chironomidae_tsl <- chironomidae_slope_draws %>%
  group_by(start_yr, TSL) %>%
  summarise(
    Pr_pos = mean(slope > 0),
    Pr_neg = mean(slope < 0),
    .groups = "drop"
  )%>%
  mutate(
    direction = case_when(
      Pr_pos > 0.95 ~ "positive",
      Pr_neg > 0.95 ~ "negative",
      TRUE          ~ "weak"
    ),
    Pr_direction = pmax(Pr_pos, Pr_neg)
  )%>%
  filter(Pr_direction > 0.95)%>%
  arrange(desc(Pr_direction))

#Acari####
aca_datA=arth_datA%>%filter(HoyeTaxon=="Acari")

#extract peak DOY
aca_fit=readRDS("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\models\\phenology_Acari.rds")

acadraws_doy=aca_fit$draws(variable="mu", format="df")

acayears <- sort(unique(aca_datA$Year))

acapeak_df <- as_draws_df(acadraws_doy)%>%
  mutate(.draw = row_number())%>%
  tidyr::pivot_longer(
    cols = starts_with("mu"),
    names_pattern = "mu\\[(\\d+)\\]",
    names_to = "year_index",
    values_to = "DOY_peak_std")%>%select(year_index, DOY_peak_std, .draw)%>%
  mutate(DOY_peak=DOY_peak_std*DOY_sd+DOY_mean,
         year_index = as.integer(year_index),
         year = acayears[year_index])

acasummary_peak <- acapeak_df%>%
  group_by(year)%>%
  summarise(
    mean  = mean(DOY_peak),
    median= median(DOY_peak),
    lower= quantile(DOY_peak, 0.05),
    upper= quantile(DOY_peak, 0.95))

step_size=1
nyr=length(unique(aca_datA$year))

years_all <- sort(unique(acapeak_df$year))

acap_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    acap_slide_draws[[length(acap_slide_draws) + 1]] <-
      acapeak_df %>% filter(year %in% yrs)
  }}

aca_slopes_draws <- lapply(acap_slide_draws, fit_slopes_per_window)

acari_slope_draws <- map2_dfr(aca_slopes_draws,acap_slide_draws, summarize_slopes)
  # mutate(TSL = end_yr - start_yr + 1,
  #        CI_width = slope_upr - slope_lwr,
  #        year = acayears[start_yr])

acap=ggplot(acari_peak_slope_ci,
            aes(x = TSL, y = start_yr)) +
  geom_point(aes(size = CI_width,fill = slope_mean),
             shape = 21, color = "black",alpha = 0.7) +
  scale_size_continuous(name = "95% CI width",range = c(1, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,
                       name = "Trend (slope)") +
  theme_classic() +
  labs(x = "Time-series length (years)",y = "Start year",
       title = "Acari",
       subtitle="Trend uncertainty across different time windows") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(face = "bold"))

#ARTH PLOTS####

ggarrange(ichp, lycp)
ggarrange(linp, ncol=2)

ggarrange(cocp, nymp)
ggarrange(scip, musp)

ggarrange(acap, phop)
ggarrange(scip, colp)
ggarrange(colp, ncol=2)

ggarrange(chip, ncol=2)

#spring AT####
step_size=1
nyr=length(unique(airtemp_dat_apr_jul$year))

years_all <- sort(unique(airtemp_dat_apr_jul$year))

springtemp_slide_draws <- list()

for (nyr in seq(5, 27, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    springtemp_slide_draws[[length(springtemp_slide_draws) + 1]] <-
      airtemp_dat_apr_jul %>% filter(year %in% yrs)
  }}

#spring temp (Apr-May)
spring_draws <- lapply(springtemp_slide_draws,get_slope_and_pval_per_window,
                       covar = "spring_mean_temp")

spring_draws_df=bind_rows(lapply(springtemp_slide_draws,
                                 get_slope_and_pval_per_window,covar = "spring_mean_temp"))%>%
  mutate(significance=if_else(pval<0.05, "S", "NS"),
         slope_sign = if_else(slope > 0, "Positive", "Negative"),
         fill_slope = case_when(
           significance == "S" & slope > 0 ~ "Positive",
           significance == "S" & slope < 0 ~ "Negative",
           TRUE ~ NA_character_))

#summer AT####
#summer temp(June-Aug)
summertemp_slide_draws <- list()

for (nyr in seq(5, 27, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    summertemp_slide_draws[[length(summertemp_slide_draws) + 1]] <-
      airtemp_dat_apr_jul %>% filter(year %in% yrs)
  }}

summer_draws <- lapply(summertemp_slide_draws,get_slope_and_pval_per_window,
                       covar = "summer_mean_temp")

summer_draws_df=bind_rows(lapply(summertemp_slide_draws,
                                 get_slope_and_pval_per_window,covar = "summer_mean_temp"))%>%
  mutate(significance=if_else(pval<0.05, "S", "NS"),
         slope_sign = if_else(slope > 0, "Positive", "Negative"),
         fill_slope = case_when(
           significance == "S" & slope > 0 ~ "Positive",
           significance == "S" & slope < 0 ~ "Negative",
           TRUE ~ NA_character_))


summerp=ggplot(summer_draws_df, aes(x = n_years, y = start_year)) +
  geom_point(aes(size = abs(slope),color = slope_sign, fill=fill_slope),
             shape = 21,alpha = 0.8,stroke = 1) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue"),
                     name = "Trend (slope)") +
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue"),
                    guide = "none",na.value = NA) +
  scale_size_continuous(name = "Slope magnitude",range = c(1, 6)) +
  theme_classic() +
  labs(
    x = "Time-series length (years)",
    y = "Start year",
    title = "Summer air temperature",
    subtitle = "Trend uncertainty across different time windows") +
  theme(plot.title = element_text(face = "bold"))

springp=ggplot(spring_draws_df, aes(x = n_years, y = start_year)) +
  geom_point(aes(size = abs(slope),color = slope_sign, fill=fill_slope),
             shape = 21,alpha = 0.8,stroke = 1) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue"),
                     name = "Trend (slope)") +
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue"),
                    guide = "none",na.value = NA) +
  scale_size_continuous(name = "Slope magnitude",range = c(1, 6)) +
  theme_classic() +
  labs(
    x = "Time-series length (years)",
    y = "Start year",
    title = "Spring air temperature",
    subtitle = "Trend uncertainty across different time windows") +
  theme(plot.title = element_text(face = "bold"))
#snow cover####
step_size=1
nyr=length(unique(snow_ave$Year))

years_all <- sort(unique(snow_ave$Year))

snow_slide_draws <- list()

for (nyr in seq(5, 29, by = 1)) {
  for (start in seq(1, length(years_all) - nyr + 1)) {

    yrs <- years_all[start:(start + nyr - 1)]

    snow_slide_draws[[length(snow_slide_draws) + 1]] <-
      snow_ave %>% filter(Year %in% yrs)
  }}

#ave. snow cover in monitoring area (June 10)
snow_draws <- lapply(snow_slide_draws,get_slope_and_pval_per_window,
                     covar = "Jun10_cover")

snow_draws_df=bind_rows(lapply(snow_slide_draws,
                               get_slope_and_pval_per_window,covar = "Jun10_cover"))%>%
  mutate(significance=if_else(pval<0.05, "S", "NS"),
         slope_sign = if_else(slope > 0, "Positive", "Negative"),
         fill_slope = case_when(
           significance == "S" & slope > 0 ~ "Positive",
           significance == "S" & slope < 0 ~ "Negative",
           TRUE ~ NA_character_))

snowp=ggplot(snow_draws_df, aes(x = n_Years, y = start_Year)) +
  geom_point(aes(size = abs(slope),color = slope_sign, fill=fill_slope),
             shape = 21,alpha = 0.8,stroke = 1) +
  scale_color_manual(values = c("Positive" = "red", "Negative" = "blue"),
                     name = "Trend (slope)") +
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue"),
                    guide = "none",na.value = NA) +
  scale_size_continuous(name = "Slope magnitude",range = c(1, 6)) +
  theme_classic() +
  labs(
    x = "Time-series length (years)",
    y = "Start year",
    title = "Spring Snow Cover",
    subtitle = "Trend uncertainty across different time windows") +
  theme(plot.title = element_text(face = "bold"))

#COVARS PLOTS####
ggarrange(springp,summerp, snowp, ncol = 3)

allt=ggarrange(dryp,musp, silp, ichp, papp, chip,
                common.legend = T, ncol=2, nrow=3, legend = "right")

annotate_figure(allt,
                left = text_grob("Start year", rot = 90, size = 12),
                bottom = text_grob("Time-series length (years)", size = 12))
