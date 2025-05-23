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

#polynomial term for DOY
dry_datA$DOYsq=dry_datA$DOY^2
#dry_datA$doys_sq11=poly(dry_datA$DOY,2, raw = T)[,2]
#dry_datA$doys_sq=poly(dry_datA$DOY,2, raw = T)

#standardize variables
dry_datA$years = (dry_datA$year - mean(dry_datA$year))/sd(dry_datA$year)
dry_datA$DOYs = (dry_datA$DOY - mean(dry_datA$DOY))/sd(dry_datA$DOY)
dry_datA$DOYsqs = (dry_datA$DOYsq - mean(dry_datA$DOYsq))/sd(dry_datA$DOYsq)
#dry_datA$doys_sqs1 = (dry_datA$doys_sq[,1] - mean(dry_datA$doys_sq[,1]))/sd(dry_datA$doys_sq[,1])
#dry_datA$doys_sqs2 = (dry_datA$doys_sq[,2] - mean(dry_datA$doys_sq[,2]))/sd(dry_datA$doys_sq[,2])
#dry_datA$doys_sqs11 = (dry_datA$doys_sq11 - mean(dry_datA$doys_sq11))/sd(dry_datA$doys_sq11)

#Data Analysis####

s1=glmer(cbind(tot_flwr, tot_NF)~  DOYs+ DOYsqs + years +
           DOYs * years+
           DOYsqs * years+
           (1 |Plot),
         family = "binomial", data=dry_datA)

#plot predictions
s1_preds=as.vector(predict(s1, type="response"))%>%cbind(dry_datA)%>%
  replace_na(list(prop_flwr=0))

colnames(s1_preds)[1]="preds"

ggplot(s1_preds, aes(x=DOY, y=preds, col=as.factor(year)))+geom_line()+
  ylab("predicted proportions")+xlab("Day of Year (DOY)")+
  theme_classic()+ggtitle("Dryas (flowers)")+scale_color_viridis_d()+facet_wrap(~Plot)+
  geom_point(aes(x=DOY, y=prop_flwr, col=as.factor(year)), size=0.85)

#model diagnostics
AIC(s1)
ModelMetrics::rmse(s1_preds$preds, s1_preds$prop_flwr)

#determine change in peak dates using derivatives

#given y=ax^2+bx+c, we can assume that change in (a) over years can
#tell us if the curve is widening or narrowing,
#and the point where the 2nd derivative (2ax+b)=0 is the peak

#to get the 2nd derivative wrt to DOY, we assume that
#ax^2=DOY^2(B1)+DOY^2:year(B4)*year, bx=DOY(B2)+DOY:year(B5)*year

yrs=1996:2024
# Standardize years the same way as in the model
yrs_std = (yrs - mean(dry_datA$year)) / sd(dry_datA$year)

coefs=fixef(s1)

b0 <- coefs['(Intercept)']
b1 <- coefs['DOYs']
b2 <- coefs['DOYsqs']
b3 <- coefs['years']
b4 <- coefs['DOYs:years']
b5 <- coefs['DOYsqs:years']

bt = b1 + b4 * yrs_std
ct = b2 + b5 * yrs_std
peak_std = -bt / (2 * ct)

# Convert back from standardized DOY to actual DOY

peak_doy = peak_std * sd(dry_datA$DOY) + mean(dry_datA$DOY)

par(mfrow=c(1,2))
plot(peak_doy ~ yrs, type="l", ylab="Peak DOY", xlab="Year", main="Dryas Estimated Phenological Peak")
plot(ct~ yrs, type="l", ylab="c", xlab="Year", main="Dryas Phenological Shape")

# Check if predicted peaks match observed peaks roughly
ggplot(s1_preds, aes(x = DOY, y = preds, group = year, col = as.factor(year))) +
  geom_line() +
  facet_wrap(~Plot) +
  theme_minimal() +
  labs(title = "Predicted Flowering Curve per Year",
       y = "Predicted Flower Proportion")+
  theme_classic()+scale_color_viridis_d()+
  geom_point(aes(x=DOY, y=prop_flwr, col=as.factor(year)), size=0.85)
