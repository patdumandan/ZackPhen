require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)
require(lubridate)
require(date)

source("M:\\My Drive\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

phen_dat_all=read.csv("M:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_phenology_metrics_1996-2023.csv")

dryas_dat=phen_dat_all%>%filter(Species=="Dryas", !Plot%in%c("Dry7", "Dry8"))%>%select(-1)

dry10=dryas_dat%>%filter(metric==10)
dry10$Plot=as.factor(dry10$Plot)

#Dryas10%

step_size=1
nyr=length(unique(dry_datA$year))
dryp_5y_slide <- list()

# Create sliding windows
for (nyr in seq(5,20, by=1)) {

  window_size <- nyr  # size of the sliding window

  for (start in seq(1, nrow(summary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    dryp_5y_slide[[length(dryp_5y_slide) + 1]] <- summary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

slopes_list <- increasing_mod(dryp_5y_slide)%>%unlist()

#get data for plotting
get_start_years <- function(slide_list) {
  lapply(slide_list, function(df) {

    start_year = as.numeric(min(df$year, na.rm = TRUE))
    return(start_year)

  })
}

get_end_years <- function(slide_list) {
  lapply(slide_list, function(df) {

    end_year = as.numeric(max(df$year, na.rm = TRUE))
    return(end_year)

  })
}

dryp_syrs=get_start_years(dryp_5y_slide)%>%unlist()
dryp_eyrs=get_end_years(dryp_5y_slide)%>%unlist()


dryas_peak_slope=cbind(dryp_syrs,dryp_eyrs, slopes_list)%>%as.data.frame()

colnames(dryas_peak_slope)= c("start_yr", "end_yr", "slope")

dryas_peak_slope$start_yr=as.integer(dryas_peak_slope$start_yr)
dryas_peak_slope$end_yr=as.integer(dryas_peak_slope$end_yr)

dryas_peak_slope=dryas_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

ggpubr::ggviolin(dryas_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Dryas (peak)")+geom_hline(yintercept=0, lty=2)
