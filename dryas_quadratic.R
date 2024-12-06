#plants
require(tidyr)
require(dplyr)
require(lubridate)

dry_dat=read.csv("k:\\My Drive\\SLU\\phenology-project\\ZackPhen\\data\\raw\\Dryas phenology_10.17897_JSQ7-6355_data.txt",
                 sep="\t", header=T)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")
#remove before 1997 because totals were across all sections
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
dry_datA$years = (dry_datA$year - mean(dry_datA$year))/(2 *sd(dry_datA$year))
dry_datA$DOYs = (dry_datA$DOY - mean(dry_datA$DOY))/(2 *sd(dry_datA$DOY))
dry_datA$periods=as.factor(dry_datA$period)

s1=lme4::glmer(value~ DOYs+ poly(DOYs, 2) + years + periods +
              DOYs*years+
              years*periods+
              DOYs*periods +
              poly(DOYs, 2)* periods+
              poly(DOYs, 2)* years+
              DOYs*years*periods+
              poly(DOYs, 2)*years*periods+
               (1|Plot),
               family = "binomial",data=dry_datA)

summary(s1)
AIC(s1)

s1_preds=as.vector(predict(s1, type="response"))%>%cbind(dry_datA)
colnames(s1_preds)[1]="preds"

#marginalize over plots (to smooth out the liens in the ggplot)
ggplot(s1_preds, aes(x=DOY, y=preds, col=period))+geom_line()+facet_wrap(~Plot)+
  theme_classic()+ggtitle("quadratic model")+
  geom_point(aes(x=DOY, y=value, col=period))
