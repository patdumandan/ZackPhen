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

ggplot(dryas_peak_slope_ci,
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

