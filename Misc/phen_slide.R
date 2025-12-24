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
#salsummary_peak=salsummary_peak%>%mutate(species="Salix")
#saxsummary_peak=saxsummary_peak%>%mutate(species="Saxifraga")
papsummary_peak=papsummary_peak%>%mutate(species="Papaver")
cassummary_peak=cassummary_peak%>%mutate(species="Cassiope")
drysummary_peak=summary_peak%>%mutate(species="Dryas")

plant_phen_all=vctrs::vec_rbind(silsummary_peak,papsummary_peak,
                                cassummary_peak,drysummary_peak)

#Dryas

step_size=1
nyr=length(unique(dry_datA$year))
dryp_5y_slide <- list()

# Create sliding windows
for (nyr in seq(5,20, by=1)) { #increase TS length by 1, starting from 5 points to 20 datapoints

  window_size <- nyr  #set the max window size

  #create increasing TS length
  for (start in seq(1, nrow(drysummary_peak) - window_size + 1, by = step_size)) {
   #for each row in the dataset, from 1 to the second to the last datapoint,
    #increase by step size (i.e.,1)
    end <- start + window_size - 1
    #slide through the dataset, one year at a time
    dryp_5y_slide[[length(dryp_5y_slide) + 1]] <- drysummary_peak[start:end, ]

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

ggpubr::ggarrange(drytsl, castsl, siltsl, paptsl, common.legend = T)

#muscid
step_size=1
musnyr=length(unique(mussummary_peak$year))

muspp_5y_slide <- list()

# Create sliding windows
for (musnyr in seq(5,20, by=1)) {

  window_size <- musnyr  # size of the sliding window

  for (start in seq(1, nrow(mussummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    muspp_5y_slide[[length(muspp_5y_slide) + 1]] <- mussummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

musslopes_list <- increasing_mod(muspp_5y_slide)%>%unlist()

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

musp_syrs=get_start_years(muspp_5y_slide)%>%unlist()
musp_eyrs=get_end_years(muspp_5y_slide)%>%unlist()


mus_peak_slope=cbind(musp_syrs,musp_eyrs, musslopes_list)%>%as.data.frame()

colnames(mus_peak_slope)= c("start_yr", "end_yr", "slope")

mus_peak_slope$start_yr=as.integer(mus_peak_slope$start_yr)
mus_peak_slope$end_yr=as.integer(mus_peak_slope$end_yr)

mus_peak_slope=mus_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

mustsl=ggpubr::ggviolin(mus_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Muscid (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Lycosidae
step_size=1
lycnyr=length(unique(lycsummary_peak$year))

lycpp_5y_slide <- list()

# Create sliding windows
for (lycnyr in seq(5,20, by=1)) {

  window_size <- lycnyr  # size of the sliding window

  for (start in seq(1, nrow(lycsummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    lycpp_5y_slide[[length(lycpp_5y_slide) + 1]] <- lycsummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

lycslopes_list <- increasing_mod(lycpp_5y_slide)%>%unlist()

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

lycp_syrs=get_start_years(lycpp_5y_slide)%>%unlist()
lycp_eyrs=get_end_years(lycpp_5y_slide)%>%unlist()


lyc_peak_slope=cbind(lycp_syrs,lycp_eyrs, lycslopes_list)%>%as.data.frame()

colnames(lyc_peak_slope)= c("start_yr", "end_yr", "slope")

lyc_peak_slope$start_yr=as.integer(lyc_peak_slope$start_yr)
lyc_peak_slope$end_yr=as.integer(lyc_peak_slope$end_yr)

lyc_peak_slope=lyc_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

lyctsl=ggpubr::ggviolin(lyc_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Lycosidae (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Phoridae
step_size=1
phocnyr=length(unique(phosummary_peak$year))

phocpp_5y_slide <- list()

# Create sliding windows
for (phocnyr in seq(5,20, by=1)) {

  window_size <- phocnyr  # size of the sliding window

  for (start in seq(1, nrow(phosummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    phocpp_5y_slide[[length(phocpp_5y_slide) + 1]] <- phosummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

phocslopes_list <- increasing_mod(phocpp_5y_slide)%>%unlist()

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

phocp_syrs=get_start_years(phocpp_5y_slide)%>%unlist()
phocp_eyrs=get_end_years(phocpp_5y_slide)%>%unlist()


phoc_peak_slope=cbind(phocp_syrs,phocp_eyrs, phocslopes_list)%>%as.data.frame()

colnames(phoc_peak_slope)= c("start_yr", "end_yr", "slope")

phoc_peak_slope$start_yr=as.integer(phoc_peak_slope$start_yr)
phoc_peak_slope$end_yr=as.integer(phoc_peak_slope$end_yr)

phoc_peak_slope=phoc_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

phoctsl=ggpubr::ggviolin(phoc_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Phoridae (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Sciaridae
step_size=1
scicnyr=length(unique(scisummary_peak$year))

scicpp_5y_slide <- list()

# Create sliding windows
for (scicnyr in seq(5,20, by=1)) {

  window_size <- scicnyr  # size of the sliding window

  for (start in seq(1, nrow(scisummary_peak) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    scicpp_5y_slide[[length(scicpp_5y_slide) + 1]] <- scisummary_peak[start:end, ]

  }
}

#fit model
increasing_mod=function(slide_list) {
  lapply(slide_list, function(dat) {
    slope <-as.numeric(coef(glm(mean ~ year, data = dat))[2])
    return(slope)
  })
}

scicslopes_list <- increasing_mod(scicpp_5y_slide)%>%unlist()

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

scicp_syrs=get_start_years(scicpp_5y_slide)%>%unlist()
scicp_eyrs=get_end_years(scicpp_5y_slide)%>%unlist()


scic_peak_slope=cbind(scicp_syrs,scicp_eyrs, scicslopes_list)%>%as.data.frame()

colnames(scic_peak_slope)= c("start_yr", "end_yr", "slope")

scic_peak_slope$start_yr=as.integer(scic_peak_slope$start_yr)
scic_peak_slope$end_yr=as.integer(scic_peak_slope$end_yr)

scic_peak_slope=scic_peak_slope%>%
  mutate(TSL=(end_yr-start_yr+1))

scictsl=ggpubr::ggviolin(scic_peak_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Sciaridae (peak)")+geom_hline(yintercept=0, lty=2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
