#plants
require(tidyr)
require(dplyr)
require(lubridate)

dry_dat=read.csv("L:\\My Drive\\SLU\\phenology-project\\ZackPhen\\data\\raw\\Dryas phenology_10.17897_JSQ7-6355_data.txt",
                 sep="\t", header=T)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")

dry_dat_raw=dry_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Dryas")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date, -Larvae, -Eaten)%>%
  filter(!Plot%in%c("Dry7", "Dry8"))

dry_dat_raw[dry_dat_raw==-9999] <- NA
dry_dat_raw[dry_dat_raw==-999] <- NA

dry_dat_raw=dry_dat_raw%>%replace_na(list(Buds=0, Flowers=0, Senescent=0))

dry_dat_grp=dry_dat_raw%>%
  group_by(Plot, year, DOY)%>%summarise(tot_bud=sum(Buds),
                                     tot_flwr=sum(Flowers),
                                     tot_sen=sum(Senescent),
                                     tot_all=sum(tot_bud, tot_flwr, tot_sen))%>%
  group_by(Plot, year)%>%mutate(tot_yr=sum(tot_all),
                                tot_buds=sum(tot_bud),
                                tot_flwrs=sum(tot_flwr),
                                tot_sens=sum(tot_sen))%>%
  group_by(Plot, year, DOY)%>%
  mutate(prop_bud=tot_bud/tot_all,
         prop_flwr=tot_flwr/tot_all,
         prop_sen=tot_sen/tot_all,
         cumprop_10=(tot_bud)/tot_yr,
         cumprop_50=(tot_flwr)/tot_yr,
         cumprop_90=(tot_sen)/tot_yr,
         props_10=tot_bud/tot_buds,
         props_50=tot_flwr/tot_flwrs,
         props_90=tot_sen/tot_sens)

dry_datA=dry_dat_grp%>%
  select(,-c(15:20))%>%
  pivot_longer(cols=c(12:14),names_to="period")%>%
  replace_na(list(value=0))

require(ggplot2)

ggplot(dry_datA, aes(x=DOY, y=value, col=period))+geom_point()+facet_wrap(~Plot)+
  theme_classic()+stat_smooth(method="gam")+ggtitle("A")

#for scaling (look into 1 instead of 2 SDs)
dry_datA$value[dry_datA$value==-NaN] <- 0
dry_datA$value[dry_datA$value==-NA] <- 0

min(dry_datA$value)
max(dry_datA$value, na.rm = T)

#for scaling (look into 1 instead of 2 SDs)
dry_datA$years = (dry_datA$year - mean(dry_datA$year))/(1 *sd(dry_datA$year))
dry_datA$DOYs = (dry_datA$DOY - mean(dry_datA$DOY))/(1 *sd(dry_datA$DOY))
dry_datA$periods=as.factor(dry_datA$period)

require(lme4)

#Dryas 1####
dry_datA1=dry_datA%>%filter(Plot=="Dry1")

s1=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA1, weights=tot_all)

summary(s1)
AIC(s1)

s1_preds=as.vector(predict(s1, type="response"))%>%cbind(dry_datA1)
colnames(s1_preds)[1]="preds"

s1_preds_bud=s1_preds%>%filter(period=="prop_bud")
s1_preds_flwr=s1_preds%>%filter(period=="prop_flwr")
s1_preds_sen=s1_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s1b=ggplot(s1_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 1 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s1f=ggplot(s1_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 1 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s1s=ggplot(s1_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 1 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

#Dryas 2####
dry_datA2=dry_datA%>%filter(Plot=="Dry2")

s2=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA2, weights=tot_all)

summary(s2)
AIC(s2)

s2_preds=as.vector(predict(s2, type="response"))%>%cbind(dry_datA2)
colnames(s2_preds)[1]="preds"

s2_preds_bud=s2_preds%>%filter(period=="prop_bud")
s2_preds_flwr=s2_preds%>%filter(period=="prop_flwr")
s2_preds_sen=s2_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s2b=ggplot(s2_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 2 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s2f=ggplot(s2_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 2 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s2s=ggplot(s1_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 2 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

#Dryas 3####

dry_datA3=dry_datA%>%filter(Plot=="Dry3")

s3=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA3, weights=tot_all)

summary(s3)
AIC(s3)

s3_preds=as.vector(predict(s3, type="response"))%>%cbind(dry_datA3)
colnames(s3_preds)[1]="preds"

s3_preds_bud=s3_preds%>%filter(period=="prop_bud")
s3_preds_flwr=s3_preds%>%filter(period=="prop_flwr")
s3_preds_sen=s3_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s3b=ggplot(s3_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 3 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s3f=ggplot(s3_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 3 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s3s=ggplot(s3_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 3 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

#Dryas 4####

dry_datA4=dry_datA%>%filter(Plot=="Dry4")

s4=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA4, weights=tot_all)

summary(s4)
AIC(s4)

s4_preds=as.vector(predict(s4, type="response"))%>%cbind(dry_datA4)
colnames(s4_preds)[1]="preds"

s4_preds_bud=s4_preds%>%filter(period=="prop_bud")
s4_preds_flwr=s4_preds%>%filter(period=="prop_flwr")
s4_preds_sen=s4_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s4b=ggplot(s4_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 4 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s4f=ggplot(s4_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 4 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s4s=ggplot(s4_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 4 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

#Dryas 5####

dry_datA5=dry_datA%>%filter(Plot=="Dry5")

s5=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA5, weights=tot_all)

summary(s5)
AIC(s5)

s5_preds=as.vector(predict(s5, type="response"))%>%cbind(dry_datA5)
colnames(s5_preds)[1]="preds"

s5_preds_bud=s5_preds%>%filter(period=="prop_bud")
s5_preds_flwr=s5_preds%>%filter(period=="prop_flwr")
s5_preds_sen=s5_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s5b=ggplot(s5_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 5 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s5f=ggplot(s5_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 5 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s5s=ggplot(s5_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 5 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

#Dryas 6####

dry_datA6=dry_datA%>%filter(Plot=="Dry6")

s6=glm(value~ DOYs+ poly(DOYs,2) + years + periods + # base linear terms
         DOYs*years + # timing and long-term trend interaction
         DOYs*periods + # timing and life stage interaction
         years*periods + #long-term trend and life stage interaction
         poly(DOYs, 2)* years + #curve and long-term trend interaction
         poly(DOYs, 2) * periods + #curve and life stage interaction
         DOYs*years*periods + #timing, long-term trend and life stage interaction
         poly(DOYs, 2)*years*periods, #curve, long-term trend, and life-stage interaction
       family = "binomial", data=dry_datA6, weights=tot_all)

summary(s6)
AIC(s6)

s6_preds=as.vector(predict(s6, type="response"))%>%cbind(dry_datA6)
colnames(s6_preds)[1]="preds"

s6_preds_bud=s6_preds%>%filter(period=="prop_bud")
s6_preds_flwr=s6_preds%>%filter(period=="prop_flwr")
s6_preds_sen=s6_preds%>%filter(period=="prop_sen")

#marginalize over plots (to smooth out the liens in the ggplot)

s6b=ggplot(s6_preds_bud, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 6 (buds)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s6f=ggplot(s6_preds_flwr, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 6 (flowers)")
#  geom_point(aes(x=DOY, y=value, col=as.factor(year)))

s6s=ggplot(s6_preds_sen, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  theme_classic()+ggtitle("Dry 6 (senescent)")
# geom_point(aes(x=DOY, y=value, col=as.factor(year)))

require(ggpubr)

ggarrange(s1b, s1f, s1s,s2b, s2f, s2s,
          s3b, s3f, s3s,s4b, s4f, s4s,
          s5b, s5f, s5s,s6b, s6f, s6s,
          ncol=3,
          nrow=6, common.legend = T, legend="right")
