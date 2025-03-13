require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)
require(lubridate)
require(date)

source("L:\\My Drive\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

phen_dat_all=read.csv("L:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_phenology_metrics_1996-2023.csv")

dryas_dat=phen_dat_all%>%filter(Species=="Dryas", !Plot%in%c("Dry7", "Dry8"))%>%select(-1)

dry10=dryas_dat%>%filter(metric==10)
dry10$Plot=as.factor(dry10$Plot)

#Dryas10%

# Parameters
nplot=6
nyr=10

window_size <- nplot*nyr  # size of the sliding window
step_size <- nplot     # step size for the sliding window

# Initialize a list to store the results
dry10_5y_slide <- list()

# Create sliding windows
for (nyr in seq(5,15, by=1)) {

  window_size <- nplot*nyr  # size of the sliding window

  for (start in seq(1, nrow(dry10) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    dry10_5y_slide[[length(dry10_5y_slide) + 1]] <- dry10[start:end, ]

  }
}

#fit model
increasing_mod=function(dat) {

  fit_model= lme4::lmer(dat[, "DOY"]~dat[, "Year"]+
                          (1|Plot),data=dat)

}

#get data for plotting
get_start_yr=function(dat) {

  start_yr=min(dat[,"Year"])
}

get_end_yr=function(dat) {

  end_yr=max(dat[,"Year"])
}

dry10_5slide=map(dry10_5y_slide, increasing_mod)
yr_slide1=map(dry10_5y_slide, get_start_yr)%>%as.matrix()
yr_slide2=map(dry10_5y_slide, get_end_yr)%>%as.matrix()

dry10_5slide_slope=map(dry10_5slide, fixef)%>%map_dbl(2)%>%as.data.frame()%>%
  cbind(yr_slide1, yr_slide2)

colnames(dry10_5slide_slope)= c("slope", "start_yr", "end_yr")

dry10_5slide_slope$start_yr=as.integer(dry10_5slide_slope$start_yr)
dry10_5slide_slope$end_yr=as.integer(dry10_5slide_slope$end_yr)

dry10_5slide_slope=dry10_5slide_slope%>%
  mutate(TSL=(end_yr-start_yr+1), metric="10")
str(dry10_5slide_slope)


ggpubr::ggviolin(dry10_5slide_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Dryas (10%)")

#Dryas (50%)
dry50=dryas_dat%>%filter(metric==50)
dry50$Plot=as.factor(dry50$Plot)

#Dryas50%

# Parameters
nplot=6
nyr=10

window_size <- nplot*nyr  # size of the sliding window
step_size <- nplot     # step size for the sliding window

# Initialize a list to store the results
dry50_5y_slide <- list()

# Create sliding windows
for (nyr in seq(5,15, by=1)) {

  window_size <- nplot*nyr  # size of the sliding window

  for (start in seq(1, nrow(dry50) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    dry50_5y_slide[[length(dry50_5y_slide) + 1]] <- dry50[start:end, ]

  }
}

#fit model
increasing_mod=function(dat) {

  fit_model= lme4::lmer(dat[, "DOY"]~dat[, "Year"]+
                          (1|Plot),data=dat)

}

#get data for plotting
get_start_yr=function(dat) {

  start_yr=min(dat[,"Year"])
}

get_end_yr=function(dat) {

  end_yr=max(dat[,"Year"])
}

dry50_5slide=map(dry50_5y_slide, increasing_mod)
yr_slide1=map(dry50_5y_slide, get_start_yr)%>%as.matrix()
yr_slide2=map(dry50_5y_slide, get_end_yr)%>%as.matrix()

dry50_5slide_slope=map(dry50_5slide, fixef)%>%map_dbl(2)%>%as.data.frame()%>%
  cbind(yr_slide1, yr_slide2)

colnames(dry50_5slide_slope)= c("slope", "start_yr", "end_yr")

dry50_5slide_slope$start_yr=as.integer(dry50_5slide_slope$start_yr)
dry50_5slide_slope$end_yr=as.integer(dry50_5slide_slope$end_yr)

dry50_5slide_slope=dry50_5slide_slope%>%
  mutate(TSL=(end_yr-start_yr+1), metric="50")

str(dry50_5slide_slope)


ggpubr::ggviolin(dry50_5slide_slope, x="start_yr", y="slope", add="jitter", color = "TSL")+
  ggtitle("Dryas (50%)")

#Dryas (90%)
dry90=dryas_dat%>%filter(metric==90)
dry90$Plot=as.factor(dry90$Plot)

#Dryas90%

# Parameters
nplot=6
nyr=10

window_size <- nplot*nyr  # size of the sliding window
step_size <- nplot     # step size for the sliding window

# Initialize a list to store the results
dry90_5y_slide <- list()

# Create sliding windows
for (nyr in seq(5,15, by=1)) {

  window_size <- nplot*nyr  # size of the sliding window

  for (start in seq(1, nrow(dry90) - window_size + 1, by = step_size)) {

    end <- start + window_size - 1
    dry90_5y_slide[[length(dry90_5y_slide) + 1]] <- dry90[start:end, ]

  }
}

#fit model
increasing_mod=function(dat) {

  fit_model= lme4::lmer(dat[, "DOY"]~dat[, "Year"]+
                          (1|Plot),data=dat)

}

#get data for plotting
get_start_yr=function(dat) {

  start_yr=min(dat[,"Year"])
}

get_end_yr=function(dat) {

  end_yr=max(dat[,"Year"])
}

dry90_5slide=map(dry90_5y_slide, increasing_mod)
yr_slide1=map(dry90_5y_slide, get_start_yr)%>%as.matrix()
yr_slide2=map(dry90_5y_slide, get_end_yr)%>%as.matrix()

dry90_5slide_slope=map(dry90_5slide, fixef)%>%map_dbl(2)%>%as.data.frame()%>%
  cbind(yr_slide1, yr_slide2)

colnames(dry90_5slide_slope)= c("slope", "start_yr", "end_yr")

dry90_5slide_slope$start_yr=as.integer(dry90_5slide_slope$start_yr)
dry90_5slide_slope$end_yr=as.integer(dry90_5slide_slope$end_yr)

dry90_5slide_slope=dry90_5slide_slope%>%
  mutate(TSL=(end_yr-start_yr+1), metric="90")
str(dry90_5slide_slope)


#combine all results
dry_res=rbind(dry10_5slide_slope, dry50_5slide_slope, dry90_5slide_slope)

ggpubr::ggviolin(dry_res, x="start_yr", y="slope", add="jitter", color = "TSL", fill="metric")+
  ggtitle("Dryas")

