require(rsample)
require(tidyr)
require(dplyr)
require(purrr)
require(ggplot2)
require(brms)
require(lubridate)
require(date)

source("J:\\My Drive\\SLU\\phenology-project\\ZackPhen\\plant_functions.R")

phen_dat_all=read.csv("J:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_phenology_metrics_1996-2023.csv")
foc_dat=phen_dat_all%>%filter(Species%in%c("Dryas", "Salix"))

#moving window####

#Dryas####

dryas_dat=phen_dat_all%>%filter(Species=="Dryas", !Plot%in%c("Dry7", "Dry8"))%>%select(-1)

ggplot(dryas_dat, aes(x=Year, y=DOY, col=as.factor(metric)))+geom_point()+theme_classic()+geom_line()+
  stat_smooth(method="gam")+ggtitle("Dryas")
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

#Salix
sal10=salix_dat%>%filter(metric==10)
sal50=salix_dat%>%filter(metric==50)
sal90=salix_dat%>%filter(metric==90)

#apply sliding-index to create subsets of training data at different windows
#set 5 years of training data, rolling over every year (every 12 months)

#10% bloom date####
dry10=dryas_dat%>%filter(metric==10, !Plot%in%c("Dry7", "Dry8"))
dry10$Plot=as.factor(dry10$Plot)

multi_YR_dat5=sliding_index(
  data= dry10, #all plots except 7 and 8
  index=Year,
  lookback=9, #training data length
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
  lookback=9, #training data length (including starting year so -1)
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

#Summary: metric diffs####
foc_dat_diff=foc_dat%>%select(-1)%>%pivot_wider(names_from = "metric", values_from = "DOY")%>%
  rename("DOY_10"="10", "DOY_50"="50", "DOY_90"="90")%>%
  mutate(diffs_end=DOY_90-DOY_10,
         diffs_startmid=DOY_50-DOY_10,
         diffs_midend=DOY_90-DOY_50)%>%
  select(-DOY_10, -DOY_50, -DOY_90)%>%
  pivot_longer(cols=4:6, names_to="diffs", values_to="values")

ggplot(foc_dat_diff, aes(x=Year, y=values, col=as.factor(diffs)))+
  geom_point()+
  stat_smooth(method="gam")+facet_wrap(~Species)+
  theme_classic()

f1=foc_dat_diff%>%filter(diffs=="diffs_end", Species=="Salix")
f2=foc_dat_diff%>%filter(diffs=="diffs_startmid", Species=="Salix")
f3=foc_dat_diff%>%filter(diffs=="diffs_midend", Species=="Salix")

summary(lme4::lmer(f1$values~f1$Year+(1|Plot), data=f1))
summary(lm(f1$values~f1$Year))
summary(lm(f2$values~f2$Year))
summary(lm(f3$values~f3$Year))

#Wavelet####
#Dryas####
dryas_dat=phen_dat_all%>%filter(Species=="Dryas", !Plot%in%c("Dry7", "Dry8"))%>%select(-1)

dry1=dryas_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

dry2=dryas_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

dry3=dryas_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

dry_com10=analyze.wavelet(dry1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(dry_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1),
                                 labels = unique(dry2$Year)))
wt.image(dry_com10, color.key = "quantile",main="Dryas (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.avg(dry_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#50% bloom date
dry_com50=analyze.wavelet(dry2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(dry_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1),
                                 labels = unique(dry2$Year)))
wt.image(dry_com50, color.key = "quantile",main="Dryas (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.avg(dry_com50, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
dry_com90=analyze.wavelet(dry3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(dry_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(dry3$Year), by = 1),
                                 labels = unique(dry3$Year)))
wt.image(dry_com90, color.key = "quantile",main="Dryas (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry3$Year), by = 1), labels = unique(dry3$Year)))
wt.avg(dry_com90, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Salix####
salix_dat=phen_dat_all%>%filter(Species=="Salix")%>%select(-1)

sal1=salix_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sal2=salix_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sal3=salix_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

sal_com10=analyze.wavelet(sal1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sal_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1),
                                 labels = unique(sal2$Year)))
wt.image(sal_com10, color.key = "quantile",main="salix (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.avg(sal_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#50% bloom date
sal_com50=analyze.wavelet(sal2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sal_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1),
                                 labels = unique(sal2$Year)))
wt.image(sal_com50, color.key = "quantile",main="salix (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.avg(sal_com50, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
sal_com90=analyze.wavelet(sal3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sal_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sal3$Year), by = 1),
                                 labels = unique(sal3$Year)))
wt.image(sal_com90, color.key = "quantile",main="salix (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal3$Year), by = 1), labels = unique(sal3$Year)))
wt.avg(sal_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Saxifraga####
saxifraga_dat=phen_dat_all%>%filter(Species=="Saxifraga")%>%select(-1)

sax1=saxifraga_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sax2=saxifraga_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sax3=saxifraga_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

sax_com10=analyze.wavelet(sax1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sax_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1),
                                 labels = unique(sax2$Year)))
wt.image(sax_com10, color.key = "quantile",main="saxifraga (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.avg(sax_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#50% bloom date
sax_com50=analyze.wavelet(sax2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sax_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1),
                                 labels = unique(sax2$Year)))
wt.image(sax_com50, color.key = "quantile",main="saxifraga (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.avg(sax_com50, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
sax_com90=analyze.wavelet(sax3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sax_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sax3$Year), by = 1),
                                 labels = unique(sax3$Year)))
wt.image(sax_com90, color.key = "quantile",main="saxifraga (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax3$Year), by = 1), labels = unique(sax3$Year)))
wt.avg(sax_com90, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Silene####
silene_dat=phen_dat_all%>%filter(Species=="Silene")%>%select(-1)

sil1=silene_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sil2=silene_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

sil3=silene_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

sil_com10=analyze.wavelet(sil1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sil_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1),
                                 labels = unique(sil2$Year)))
wt.image(sil_com10, color.key = "quantile",main="silene (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.avg(sil_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#50% bloom date
sil_com50=analyze.wavelet(sil2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sil_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1),
                                 labels = unique(sil2$Year)))
wt.image(sil_com50, color.key = "quantile",main="silene (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.avg(sil_com50, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
sil_com90=analyze.wavelet(sil3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(sil_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(sil3$Year), by = 1),
                                 labels = unique(sil3$Year)))
wt.image(sil_com90, color.key = "quantile",main="silene (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil3$Year), by = 1), labels = unique(sil3$Year)))
wt.avg(sil_com90, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Cassiope####
cassiope_dat=phen_dat_all%>%filter(Species=="Cassiope", !Year<1996)%>%select(-1)

cas1=cassiope_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

cas2=cassiope_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

cas3=cassiope_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

cas_com10=analyze.wavelet(cas1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(cas_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1),
                                 labels = unique(cas2$Year)))
wt.image(cas_com10, color.key = "quantile",main="cassiope (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.avg(cas_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE)

#50% bloom date
cas_com50=analyze.wavelet(cas2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(cas_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1),
                                 labels = unique(cas2$Year)))
wt.image(cas_com50, color.key = "quantile",main="cassiope (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.avg(cas_com50, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
cas_com90=analyze.wavelet(cas3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(cas_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(cas3$Year), by = 1),
                                 labels = unique(cas3$Year)))
wt.image(cas_com90, color.key = "quantile",main="cassiope (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas3$Year), by = 1), labels = unique(cas3$Year)))
wt.avg(cas_com90, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Papaver####
papaver_dat=phen_dat_all%>%filter(Species=="Papaver")%>%select(-1)

pap1=papaver_dat%>%filter(metric==10)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

pap2=papaver_dat%>%filter(metric==50)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

pap3=papaver_dat%>%filter(metric==90)%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY))

#10% bloom date

pap_com10=analyze.wavelet(pap1, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(pap_com10, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1),
                                 labels = unique(pap2$Year)))
wt.image(pap_com10, color.key = "quantile",main="papaver (10% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.avg(pap_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#50% bloom date
pap_com50=analyze.wavelet(pap2, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(pap_com50, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1),
                                 labels = unique(pap2$Year)))
wt.image(pap_com50, color.key = "quantile",main="papaver (50% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.avg(pap_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#90% bloom date
pap_com90=analyze.wavelet(pap3, "mean_doy", make.pval = TRUE, n.sim = 10000, loess.span = 0.75)
reconstruct(pap_com90, "mean_doy", show.legend = F,only.coi = T,only.sig = T,
            spec.time.axis =list(at = seq(1, length(pap3$Year), by = 1),
                                 labels = unique(pap3$Year)))
wt.image(pap_com90, color.key = "quantile",main="papaver (90% bloom date)", col.contour = "black",plot.ridge = F,
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap3$Year), by = 1), labels = unique(pap3$Year)))
wt.avg(pap_com10, siglvl = c(0.05, 0.1), sigcol = c("red", "blue"),show.siglvl = TRUE,
       periodlab = "period (years)")

#Summary: Cyclicity####
cyc=read.csv("J:\\My Drive\\SLU\\phenology-project\\ZackPhen\\cyclicity_plants.csv")%>%filter(!Year<1996)

ggplot(cyc)+geom_line(aes(x=Year, y=Species,
                          size=cyclicity, col=Species),
                      show.legend = F)+facet_wrap(~Metric)+
  theme_classic()+geom_vline(xintercept=c(1998,2021), linetype="dashed", lwd=0.5)+
  scale_y_discrete()


cyc_grp=cyc%>%filter(!Year<1996)%>%
  group_by(Species, Metric)%>%
  summarise(tot_yr=length(Year) ,
            cyc_yrs=sum(cyclicity),
            prop_cyc=cyc_yrs/tot_yr)

cyc_grp_yr=cyc%>%
  group_by(Year, Metric)%>%
  summarise(tot_yr=28,
            cyc_yrs=sum(cyclicity),
            prop_cyc=cyc_yrs/tot_yr)

#increase TS length####
create_incremental_training_sets <- function(data, max_size, increment) {
  # Create an empty list to store the training sets
  training_sets <- list()

  # Loop to create training sets of increasing size
  for (size in seq(increment, max_size, by = increment)) {
    # Ensure the size does not exceed the number of available rows
    if (size <= nrow(data)) {

      training_sets[[length(training_sets) + 1]] <- data[1:size, ]
    }
  }

  return(training_sets)
}

# Create incrementally larger training sets

# Define parameters
max_size <- length(dry10$Year)  # Maximum size of the training set
increment <- 6*5  # Incremental increase in size (5 years at a time)

dry10_ts<- create_incremental_training_sets(dry10, max_size, increment)

increasing_mod=function(dat) {

  fit_model= lme4::lmer(dat[, "DOY"]~dat[, "Year"]+
                          (1|Plot),data=dat)

}

get_yr=function(dat) {

  end_yr=max(dat[,"Year"])

  }

year1=map(dry10_ts, get_yr)

model1=map(dry10_ts, increasing_mod)
train_dat=cbind(d1, training_sets)
slope1=map(d1, fixef)%>%map_dbl(2)

s1=d1%>%map(fixef)%>%map_dbl(2)%>%as.list()
yr_start=rep(1996, 5)
training_dat=cbind(d1, s1, y1, yr_start)

