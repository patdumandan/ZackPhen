#arthropod analysis

#load packages
require(dplyr)
require(tidyr)
require(lme4)
require(ggplot2)
require(ggpubr)

#load cleana arthropod data
arth_dat=read.csv("https://raw.githubusercontent.com/patdumandan/ZackPhen/refs/heads/main/arth_raw_dat.csv", header=T)

#inspect data crudely
#plots per taxon

arth_plots=arth_dat%>%group_by(HoyeTaxon, Year)%>%summarise(plot_ids=list(unique(Plot.ID)), n_plots=length(unique(Plot.ID)))

#remove Art1, Art4 and Art6

#Art1: window trap
#Art4 and 6:discontinuous opening
#Arth 7 missing from 1996-1998?
arth_plots=arth_dat%>%filter(!Plot.ID%in%c("Art1","Art4", "Art6"))%>%group_by(HoyeTaxon, Year)%>%
  summarise(plot_ids=list(unique(Plot.ID)), n_plots=length(unique(Plot.ID)))%>%
  arrange((Year))

#final dataset

arth_data=arth_dat%>%filter(!Plot.ID%in%c("Art1","Art4", "Art6"), !Year<1996)%>%
  select(-CatchID, -UnitID, -X)%>%
  mutate(DOYsq=DOY^2)

#rescale data

arth_data$years = (arth_data$Year - mean(arth_data$Year))/sd(arth_data$Year)
arth_data$DOYs = (arth_data$DOY - mean(arth_data$DOY))/sd(arth_data$DOY)
arth_data$DOYsqs = (arth_data$DOYsq - mean(arth_data$DOYsq))/sd(arth_data$DOYsq)

unique(arth_data$Plot.ID)
unique(arth_data$HoyeTaxon)

bom_dat=arth_data%>%filter(HoyeTaxon=="Bombus")
col_dat=arth_data%>%filter(HoyeTaxon=="Collembola")
ich_dat=arth_data%>%filter(HoyeTaxon=="Ichneumonidae")
aca_dat=arth_data%>%filter(HoyeTaxon=="Acari")
chi_dat=arth_data%>%filter(HoyeTaxon=="Chironomidae")
coc_dat=arth_data%>%filter(HoyeTaxon=="Coccoidea")
cul_dat=arth_data%>%filter(HoyeTaxon=="Culicidae")
lin_dat=arth_data%>%filter(HoyeTaxon=="Linyphiidae")
mus_dat=arth_data%>%filter(HoyeTaxon=="Muscidae")
nym_dat=arth_data%>%filter(HoyeTaxon=="Nymphalidae")
pho_dat=arth_data%>%filter(HoyeTaxon=="Phoridae")
sci_dat=arth_data%>%filter(HoyeTaxon=="Sciaridae")
lyc_dat=arth_data%>%filter(HoyeTaxon=="Lycosidae")

#models####

#Collembola

col_mod_nb=glmer.nb(TotalCatch1~ DOYs+ DOYsqs + years + # base linear terms
                DOYs*years + # timing and long-term trend interaction
                DOYsqs* years+ #curve and long-term trend interaction
                (1|Plot.ID), #random effects to account for plot-level differences in trends
               data=col_dat)

col_mod_pois=glmer(TotalCatch1~ DOYs+ DOYsqs + years + # base linear terms
                      DOYs*years + # timing and long-term trend interaction
                      DOYsqs* years+ #curve and long-term trend interaction
                      (1|Plot.ID), #random effects to account for plot-level differences in trends
                    family="poisson", data=col_dat)

AIC(col_mod_nb)
AIC(col_mod_pois)

col_preds=as.vector(predict(col_mod_pois, type="response"))%>%cbind(col_dat)
colnames(col_preds)[1]="preds"

col_preds1=col_preds%>%
  rename(predicted=preds, actual=TotalCatch1)%>%
  pivot_longer(cols=c("predicted", "actual"), names_to="value_type")

ggplot(col_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Collembola")+scale_color_viridis_d()+
# geom_point(aes(x=DOY, y=TotalCatch1), col="black", size=1)+
  facet_wrap(~Plot.ID)

ggplot(col_preds1, aes(y=value, x=DOY, col=value_type))+geom_point()+facet_wrap(~Plot.ID)+
  theme_classic()

#determine change in peak dates using derivatives

#given y=ax^2+bx+c, we can assume that change in (a) over years can
#tell us if the curve is widening or narrowing,
#and the point where the 2nd derivative (2ax+b)=0 is the peak

#to get the 2nd derivative wrt to DOY, we assume that
#ax^2=DOY^2(B1)+DOY^2:year(B4)*year, bx=DOY(B2)+DOY:year(B5)*year

yrs=1996:2024

coefs=fixef(col_mod_pois)

b0 <- coefs['(Intercept)']
b1 <- coefs['DOYs']
b2 <- coefs['DOYsqs']
b3 <- coefs['years']
b4 <- coefs['DOYs:years']
b5 <- coefs['DOYsqs:years']

bt=b1+b4*yrs
ct= b2+b5*yrs
peak=bt/(2*ct)
plot(peak~yrs)

