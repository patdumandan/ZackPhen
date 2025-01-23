require(dplyr)
require(tidyr)
require(lme4)

#load cleana arthropod data
arth_dat=read.csv("L:\\My Drive\\SLU\\phenology-project\\ZackPhen\\arth_raw_dat.csv", header=T)

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

arth_data=arth_dat%>%filter(!Plot.ID%in%c("Art1","Art4", "Art6"))%>%
  select(-CatchID, -UnitID, -X)

#rescale data

arth_data$years = (arth_data$Year - mean(arth_data$Year))/(1 *sd(arth_data$Year))
arth_data$DOYs = (arth_data$DOY - mean(arth_data$DOY))/(1 *sd(arth_data$DOY))

unique(arth_data$Plot.ID)
unique(arth_data$HoyeTaxon)

bom_dat=arth_data%>%filter(HoyeTaxon=="Bombus")
col_dat=arth_data%>%filter(HoyeTaxon=="Collembola")
ich_dat=arth_data%>%filter(HoyeTaxon=="Ichneumonidae")
aca_dat=arth_data%>%filter(HoyeTaxon=="Acari")

b1=ggplot(bomb_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Bombus")

c1=ggplot(col_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Collembola")

i1=ggplot(ich_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Ichneumonidae")

a1=ggplot(aca_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Acari")

ggarrange(b1,c1,i1,a1)


#models

#Collembola

col_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
                 DOYs*years + # timing and long-term trend interaction
                 poly(DOYs, 2)* years, #curve and long-term trend interaction
               family = "poisson", data=col_dat)

col_preds=as.vector(predict(col_mod, type="response"))%>%cbind(col_dat)
colnames(col_preds)[1]="preds"

col_plot1=ggplot(col_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Collembola")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Acari
aca_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=aca_dat)

aca_preds=as.vector(predict(aca_mod, type="response"))%>%cbind(aca_dat)
colnames(aca_preds)[1]="preds"

aca_plot1=ggplot(aca_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Acari")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Ichneumonidae
ich_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=ich_dat)

ich_preds=as.vector(predict(ich_mod, type="response"))%>%cbind(ich_dat)
colnames(ich_preds)[1]="preds"

ich_plot1=ggplot(ich_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Ichneumonidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Bombus
#Acari
bom_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=bom_dat)

bom_preds=as.vector(predict(bom_mod, type="response"))%>%cbind(bom_dat)
colnames(bom_preds)[1]="preds"

bom_plot1=ggplot(bom_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Bombus")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")
ggarrange(bom_plot1, col_plot1, ich_plot1, aca_plot1)
