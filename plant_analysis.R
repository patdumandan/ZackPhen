require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)

source("J:\\My Drive\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

#focal species####
phen_dat_all=read.csv("J:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_phenology_metrics_1996-2023.csv")
foc_dat=phen_dat_all%>%filter(Species%in%c("Dryas", "Salix"))

#Dryas####

dryas_dat=foc_dat%>%filter(Species=="Dryas")%>%select(-1)

# ggplot(dryas_dat, aes(x=Year, y=DOY, col=as.factor(metric)))+geom_point()+theme_classic()+geom_line()+
#   stat_smooth(method="gam")+ggtitle("Dryas")
#
# ggplot(salix_dat, aes(x=Year, y=DOY, col=as.factor(metric)))+geom_point()+theme_classic()+geom_line()+
#   stat_smooth(method="gam")+ggtitle("Salix")
#
# ggplot(dryas_dat, aes(x=Year, y=DOY, col=as.factor(metric)))+geom_point()+theme_classic()+geom_line()+
#   stat_smooth(method="gam")+ggtitle("Dryas")+facet_wrap(~Plot)
#
# ggplot(salix_dat, aes(x=Year, y=DOY, col=as.factor(metric)))+geom_point()+theme_classic()+geom_line()+
#   stat_smooth(method="gam")+ggtitle("Salix")+facet_wrap(~Plot)

#Dryas
dry10=dryas_dat%>%filter(metric==10)
dry50=dryas_dat%>%filter(metric==50)
dry90=dryas_dat%>%filter(metric==90)

all_yr_mod=brm(DOY~Year+(1|Plot), data=dry10)
all_yr_mod50=brm(DOY~Year+(1|Plot), data=dry50)
all_yr_mod90=brm(DOY~Year+(1|Plot), data=dry90)

#Salix
sal10=salix_dat%>%filter(metric==10)
sal50=salix_dat%>%filter(metric==50)
sal90=salix_dat%>%filter(metric==90)

all_yrsal_mod=brm(DOY~Year+(1|Plot), data=sal10)
all_yrsal_mod50=brm(DOY~Year+(1|Plot), data=sal50)
all_yrsal_mod90=brm(DOY~Year+(1|Plot), data=sal90)

#apply sliding-index to create subsets of training data at different windows
#set 5 years of training data, rolling over every year (every 12 months)

#10% bloom date####
dry10=dryas_dat%>%filter(metric==10, !Plot%in%c("Dry7", "Dry8"))
dry10$Plot=as.factor(dry10$Plot)

multi_YR_dat5=sliding_index(
  data= dry10, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR_dat5$model=map(multi_YR_dat5$splits, rolling_mod)
multi_YR_dat5$slope=multi_YR_dat5$model%>%map(coef)%>%map_dbl(2)
multi_YR_dat5$intercept=multi_YR_dat5$model%>%map(coef)%>%map_dbl(1)
multi_YR_dat5$start_yr=map(multi_YR_dat5$splits, get_year)
multi_YR_dat5$plot=pmap(list(multi_YR_dat5$splits, multi_YR_dat5$slope, multi_YR_dat5$intercept), rolling_plot)

dry10_5=multi_YR_dat5%>%select(start_yr, slope)%>%mutate(metric=10)
dry10_5$start_yr=as.integer(dry10_5$start_yr)

#90% bloom date####
dry90=dryas_dat%>%filter(metric==90)
dry90$Plot=as.factor(dry90$Plot)

multi_YR90_dat5=sliding_index(
  data= dry90, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length (including starting year so -1)
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR90_dat5$model=map(multi_YR90_dat5$splits, rolling_mod)
multi_YR90_dat5$slope=multi_YR90_dat5$model%>%map(coef)%>%map_dbl(2)
multi_YR90_dat5$intercept=multi_YR90_dat5$model%>%map(coef)%>%map_dbl(1)
multi_YR90_dat5$start_yr=map(multi_YR90_dat5$splits, get_year)
multi_YR90_dat5$plot=pmap(list(multi_YR90_dat5$splits, multi_YR90_dat5$slope, multi_YR90_dat5$intercept), rolling_plot)

dry90_5=multi_YR90_dat5%>%select(start_yr, slope)%>%mutate(metric=90)
dry90_5$start_yr=as.integer(dry90_5$start_yr)

#50 bloom date####

dry50=dryas_dat%>%filter(metric==50)
dry50$Plot=as.factor(dry50$Plot)

multi_YR50_dat5=sliding_index(
  data= dry50, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length (including starting year so -1)
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR50_dat5$model=map(multi_YR50_dat5$splits, rolling_mod)
multi_YR50_dat5$slope=multi_YR50_dat5$model%>%map(coef)%>%map_dbl(2)
multi_YR50_dat5$intercept=multi_YR50_dat5$model%>%map(coef)%>%map_dbl(1)
multi_YR50_dat5$start_yr=map(multi_YR50_dat5$splits, get_year)
multi_YR50_dat5$plot=pmap(list(multi_YR50_dat5$splits, multi_YR50_dat5$slope, multi_YR50_dat5$intercept), rolling_plot)

dry50_5=multi_YR50_dat5%>%select(start_yr, slope)%>%mutate(metric=50)
dry50_5$start_yr=as.integer(dry50_5$start_yr)

#slopes####
dry_comb=rbind(dry10_5, dry90_5, dry50_5)
dry_comb$metric=as.factor(dry_comb$metric)

ggplot()+
  geom_point(data=dry_comb, aes(x=start_yr, y=slope,col=as.factor(metric)))+
  geom_line(data=dry_comb, aes(x=start_yr, y=slope, col=as.factor(metric)))+
  theme_classic()+ggtitle("Dryas")

#Salix####

salix_dat=foc_dat%>%filter(Species=="Salix")%>%select(-1)

#Salix
sal10=salix_dat%>%filter(metric==10)
sal50=salix_dat%>%filter(metric==50)
sal90=salix_dat%>%filter(metric==90)

# all_yrsal_mod=brm(DOY~Year+(1|Plot), data=sal10)
# all_yrsal_mod50=brm(DOY~Year+(1|Plot), data=sal50)
# all_yrsal_mod90=brm(DOY~Year+(1|Plot), data=sal90)

#apply sliding-index to create subsets of training data at different windows
#set 5 years of training data, rolling over every year (every 12 months)

#10% bloom date####
sal10=salix_dat%>%filter(metric==10)
sal10$Plot=as.factor(sal10$Plot)

multi_YR_saldat5=sliding_index(
  data= sal10, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR_saldat5$model=map(multi_YR_saldat5$splits, rolling_mod)
multi_YR_saldat5$slope=multi_YR_saldat5$model%>%map(coef)%>%map_dbl(2)
multi_YR_saldat5$intercept=multi_YR_saldat5$model%>%map(coef)%>%map_dbl(1)
multi_YR_saldat5$start_yr=map(multi_YR_saldat5$splits, get_year)
multi_YR_saldat5$plot=pmap(list(multi_YR_saldat5$splits, multi_YR_saldat5$slope, multi_YR_saldat5$intercept), rolling_plot)

sal10_5=multi_YR_saldat5%>%select(start_yr, slope)%>%mutate(metric=10)
sal10_5$start_yr=as.integer(sal10_5$start_yr)

#90% bloom date####
sal90=salix_dat%>%filter(metric==90)
sal90$Plot=as.factor(sal90$Plot)

multi_YR90_saldat5=sliding_index(
  data= sal90, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length (including starting year so -1)
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR90_saldat5$model=map(multi_YR90_saldat5$splits, rolling_mod)
multi_YR90_saldat5$slope=multi_YR90_saldat5$model%>%map(coef)%>%map_dbl(2)
multi_YR90_saldat5$intercept=multi_YR90_saldat5$model%>%map(coef)%>%map_dbl(1)
multi_YR90_saldat5$start_yr=map(multi_YR90_saldat5$splits, get_year)
multi_YR90_saldat5$plot=pmap(list(multi_YR90_saldat5$splits, multi_YR90_saldat5$slope, multi_YR90_saldat5$intercept), rolling_plot)

sal90_5=multi_YR90_saldat5%>%select(start_yr, slope)%>%mutate(metric=90)
sal90_5$start_yr=as.integer(sal90_5$start_yr)

#50 bloom date####

sal50=salix_dat%>%filter(metric==50)
sal50$Plot=as.factor(sal50$Plot)

multi_YR50_saldat5=sliding_index(
  data= sal50, #all plots except 7 and 8
  index=Year,
  lookback=10, #training data length (including starting year so -1)
  complete = TRUE,
  step=6 #every 6th row is a new year
)

multi_YR50_saldat5$model=map(multi_YR50_saldat5$splits, rolling_mod)
multi_YR50_saldat5$slope=multi_YR50_saldat5$model%>%map(coef)%>%map_dbl(2)
multi_YR50_saldat5$intercept=multi_YR50_saldat5$model%>%map(coef)%>%map_dbl(1)
multi_YR50_saldat5$start_yr=map(multi_YR50_saldat5$splits, get_year)
multi_YR50_saldat5$plot=pmap(list(multi_YR50_saldat5$splits, multi_YR50_saldat5$slope, multi_YR50_saldat5$intercept), rolling_plot)

sal50_5=multi_YR50_saldat5%>%select(start_yr, slope)%>%mutate(metric=50)
sal50_5$start_yr=as.integer(sal50_5$start_yr)

#slopes####
sal_comb=rbind(sal10_5, sal90_5, sal50_5)
sal_comb$metric=as.factor(sal_comb$metric)

ggplot()+
  geom_point(data=sal_comb, aes(x=start_yr, y=slope,col=as.factor(metric)))+
  geom_line(data=sal_comb, aes(x=start_yr, y=slope, col=as.factor(metric)))+
  theme_classic()+ggtitle("Salix")
