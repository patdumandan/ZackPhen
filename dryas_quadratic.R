#Plant analysis

#load packages
require(tidyr)
require(dplyr)
require(lubridate)
require(ggplot2)
require(lme4)

#Data Wrangling####
dry_dat=read.csv("https://raw.githubusercontent.com/patdumandan/ZackPhen/refs/heads/main/Dryas%20phenology_10.17897_JSQ7-6355/Dryas%20phenology_10.17897_JSQ7-6355_data.txt",
                 sep="\t", header=T)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")

dry_dat_raw=dry_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Dryas")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date, -Larvae, -Eaten)%>%
  filter(!Plot%in%c("Dry7", "Dry8"), !year<1996)

dry_dat_raw[dry_dat_raw==-9999] <- NA
dry_dat_raw[dry_dat_raw==-999] <- NA

dry_dat_raw=dry_dat_raw%>%replace_na(list(Buds=0, Flowers=0, Senescent=0))

#calculate proportions (unused in binomial model)
dry_dat_grp=dry_dat_raw%>%
  group_by(Plot, year, DOY)%>%summarise(tot_bud=sum(Buds),
                                        tot_flwr=sum(Flowers),
                                        tot_sen=sum(Senescent),
                                        tot_NF=sum(tot_bud+tot_sen),
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
  #select(,-c(15:20))%>%
  replace_na(list(value=0))

#Data Viz####
ggplot(dry_datA, aes(x=DOY, y=tot_flwr, col=year))+geom_point()+facet_wrap(~Plot)+
  theme_classic()+stat_smooth(method="gam", col="black")+ggtitle("Dryas")

#Data Analysis###

#dry_datA$doys_sq11=poly(dry_datA$DOY,2, raw = T)[,2]
#dry_datA$doys_sq=poly(dry_datA$DOY,2, raw = T)

#standardize variables
dry_datA$yearc = as.factor(dry_datA$year) # - mean(dry_datA$year))/sd(dry_datA$year)

#dry_datA$years = (dry_datA$year - mean(dry_datA$year))/sd(dry_datA$year)
dry_datA$DOYs = scale(dry_datA$DOY, center = TRUE, scale = TRUE)
#polynomial term for DOY, no need to standardize to avoid distortion of polynomial relationship,
#and interfere interpretability of interaction terms
dry_datA$DOYsqs=dry_datA$DOYs^2

#dry_datA$doys_sqs1 = (dry_datA$doys_sq[,1] - mean(dry_datA$doys_sq[,1]))/sd(dry_datA$doys_sq[,1])
#dry_datA$doys_sqs2 = (dry_datA$doys_sq[,2] - mean(dry_datA$doys_sq[,2]))/sd(dry_datA$doys_sq[,2])
#dry_datA$doys_sqs11 = (dry_datA$doys_sq11 - mean(dry_datA$doys_sq11))/sd(dry_datA$doys_sq11)

#Data Analysis####

s1=glmer(cbind(tot_flwr, tot_NF)~  DOYs+ DOYsqs + yearc +
           DOYs * yearc+
           DOYsqs * yearc+
           (1 |Plot),
         family = "binomial", data=dry_datA)

#model diagnostics
AIC(s1)
ModelMetrics::rmse(s1_preds$preds, s1_preds$prop_flwr)

#plot predictions
s1_preds=as.vector(predict(s1, type="response"))%>%cbind(dry_datA)%>%
  replace_na(list(prop_flwr=0))

colnames(s1_preds)[1]="preds"

ggplot(s1_preds, aes(x=DOY, y=preds, col=yearc))+geom_line()+
  ylab("predicted proportions")+xlab("Day of Year (DOY)")+
  theme_classic()+ggtitle("Dryas (flowers)")+scale_color_viridis_d()+facet_wrap(~Plot)+
  geom_point(aes(x=DOY, y=prop_flwr, col=yearc), size=0.85)

highlight_years <- c("1997", "2021")

ggplot(s1_preds, aes(x = DOY, y = preds, group = year, col = as.factor(year))) +
  geom_line(linewidth = 0.6, alpha = 0.5) +  # default lines for all years
  geom_line(data = subset(s1_preds, year %in% highlight_years),
            aes(x = DOY, y = preds, group = year, color = as.factor(year)),
            linewidth = 1.2) +  # bold lines for highlighted years
  geom_point(aes(y = prop_flwr), size = 0.85) +  # points for observed data
  facet_wrap(~Plot) +
  scale_color_manual(
    values = c("1997" = "orange", "2021" = "darkgreen"),
    breaks = highlight_years,
    guide = guide_legend(title = "Highlighted Years")
  ) +
  theme_classic() +
  labs(
    title = "Predicted Flowering Curve per Year",
    y = "Predicted Flower Proportion",
    color = "Year"
  )

