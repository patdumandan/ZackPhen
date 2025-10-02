require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)
require(lubridate)
require(date)

source("C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

silsummary_peak=silsummary_peak%>%mutate(species="Silene")
salsummary_peak=salsummary_peak%>%mutate(species="Salix")
saxsummary_peak=saxsummary_peak%>%mutate(species="Saxifraga")
papsummary_peak=papsummary_peak%>%mutate(species="Papaver")
cassummary_peak=cassummary_peak%>%mutate(species="Cassiope")
summary_peak=summary_peak%>%mutate(species="Dryas")

plant_phen_all=vctrs::vec_rbind(silsummary_peak, salsummary_peak,
                                saxsummary_peak, papsummary_peak,
                                cassummary_peak, summary_peak)%>%
  select(-highlight_group)
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

drytsl=ggpubr::ggviolin(dryas_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Dryas (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#cassiope
step_size=1
casnyr=length(unique(cassummary_peak$year))

casp_5y_slide <- list()

# Create sliding windows
for (casnyr in seq(5,20, by=1)) {

  window_size <- casnyr  # size of the sliding window

  for (start in seq(1, nrow(cassummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    casp_5y_slide[[length(casp_5y_slide) + 1]] <- cassummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

casslopes_list <- increasing_mod(casp_5y_slide)%>%unlist()

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

casp_syrs=get_start_years(casp_5y_slide)%>%unlist()
casp_eyrs=get_end_years(casp_5y_slide)%>%unlist()


cas_peak_slope=cbind(casp_syrs,casp_eyrs, casslopes_list)%>%as.data.frame()

colnames(cas_peak_slope)= c("start_yr", "end_yr", "slope")

cas_peak_slope$start_yr=as.integer(cas_peak_slope$start_yr)
cas_peak_slope$end_yr=as.integer(cas_peak_slope$end_yr)

cas_peak_slope=cas_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

castsl=ggpubr::ggviolin(cas_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Cassiope (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#silene
step_size=1
silnyr=length(unique(silsummary_peak$year))

silp_5y_slide <- list()

# Create sliding windows
for (silnyr in seq(5,20, by=1)) {

  window_size <- silnyr  # size of the sliding window

  for (start in seq(1, nrow(silsummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    silp_5y_slide[[length(silp_5y_slide) + 1]] <- silsummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

silslopes_list <- increasing_mod(silp_5y_slide)%>%unlist()

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

silp_syrs=get_start_years(silp_5y_slide)%>%unlist()
silp_eyrs=get_end_years(silp_5y_slide)%>%unlist()


sil_peak_slope=cbind(silp_syrs,silp_eyrs, silslopes_list)%>%as.data.frame()

colnames(sil_peak_slope)= c("start_yr", "end_yr", "slope")

sil_peak_slope$start_yr=as.integer(sil_peak_slope$start_yr)
sil_peak_slope$end_yr=as.integer(sil_peak_slope$end_yr)

sil_peak_slope=sil_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

siltsl=ggpubr::ggviolin(sil_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Silene (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#papaver
step_size=1
papnyr=length(unique(papsummary_peak$year))

papp_5y_slide <- list()

# Create sliding windows
for (papnyr in seq(5,20, by=1)) {

  window_size <- papnyr  # size of the sliding window

  for (start in seq(1, nrow(papsummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    papp_5y_slide[[length(papp_5y_slide) + 1]] <- papsummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

papslopes_list <- increasing_mod(papp_5y_slide)%>%unlist()

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

papp_syrs=get_start_years(papp_5y_slide)%>%unlist()
papp_eyrs=get_end_years(papp_5y_slide)%>%unlist()


pap_peak_slope=cbind(papp_syrs,papp_eyrs, papslopes_list)%>%as.data.frame()

colnames(pap_peak_slope)= c("start_yr", "end_yr", "slope")

pap_peak_slope$start_yr=as.integer(pap_peak_slope$start_yr)
pap_peak_slope$end_yr=as.integer(pap_peak_slope$end_yr)

pap_peak_slope=pap_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

paptsl=ggpubr::ggviolin(pap_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Papaver (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggarrange(drytsl, castsl, siltsl, paptsl, common.legend = T)

