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
chi_dat=arth_data%>%filter(HoyeTaxon=="Chironomidae")
coc_dat=arth_data%>%filter(HoyeTaxon=="Coccoidea")
cul_dat=arth_data%>%filter(HoyeTaxon=="Culicidae")
lin_dat=arth_data%>%filter(HoyeTaxon=="Linyphiidae")
mus_dat=arth_data%>%filter(HoyeTaxon=="Muscidae")
nym_dat=arth_data%>%filter(HoyeTaxon=="Nymphalidae")
pho_dat=arth_data%>%filter(HoyeTaxon=="Phoridae")
sci_dat=arth_data%>%filter(HoyeTaxon=="Sciaridae")
lyc_dat=arth_data%>%filter(HoyeTaxon=="Lycosidae")

b1=ggplot(bomb_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Bombus")

c1=ggplot(col_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Collembola")

i1=ggplot(ich_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Ichneumonidae")

a1=ggplot(aca_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Acari")

ch1=ggplot(chi_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Chironomidae")

co1=ggplot(coc_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Coccoidea")

cu1=ggplot(cul_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Culicidae")

l1=ggplot(lin_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Linyphiidae")

m1=ggplot(mus_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Muscidae")

n1=ggplot(nym_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Nymphalidae")

p1=ggplot(pho_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Phoridae")

s1=ggplot(sci_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
  geom_smooth(method="gam")+ggtitle("Sciaridae")

ly1=ggplot(lyc_dat, aes(x=DOY, y=log(TotalCatch1)))+geom_point()+facet_wrap(~Plot.ID)+
    geom_smooth(method="gam")+ggtitle("Lycosidae")

ggarrange(b1,c1,i1,a1,ly1, s1, p1, n1, m1, l1, cu1, co1, ch1)

#models####

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
bom_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=bom_dat)

bom_preds=as.vector(predict(bom_mod, type="response"))%>%cbind(bom_dat)
colnames(bom_preds)[1]="preds"

bom_plot1=ggplot(bom_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Bombus")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Lycosidae
lyc_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=lyc_dat)

lyc_preds=as.vector(predict(lyc_mod, type="response"))%>%cbind(lyc_dat)
colnames(lyc_preds)[1]="preds"

lyc_plot1=ggplot(lyc_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Lycosidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Sciaridae
sci_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=sci_dat)

sci_preds=as.vector(predict(sci_mod, type="response"))%>%cbind(sci_dat)
colnames(sci_preds)[1]="preds"

sci_plot1=ggplot(sci_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Sciaridae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Phoridae
pho_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=pho_dat)

pho_preds=as.vector(predict(pho_mod, type="response"))%>%cbind(pho_dat)
colnames(pho_preds)[1]="preds"

pho_plot1=ggplot(pho_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Phoridae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Nymphalidae
nym_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=nym_dat)

nym_preds=as.vector(predict(nym_mod, type="response"))%>%cbind(nym_dat)
colnames(nym_preds)[1]="preds"

nym_plot1=ggplot(nym_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Nymphalidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Muscidae
mus_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=mus_dat)

mus_preds=as.vector(predict(mus_mod, type="response"))%>%cbind(mus_dat)
colnames(mus_preds)[1]="preds"

mus_plot1=ggplot(mus_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Muscidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Linyphiidae
lin_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=lin_dat)

lin_preds=as.vector(predict(lin_mod, type="response"))%>%cbind(lin_dat)
colnames(lin_preds)[1]="preds"

lin_plot1=ggplot(lin_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Linyphiidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Culicidae
cul_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=cul_dat)

cul_preds=as.vector(predict(cul_mod, type="response"))%>%cbind(cul_dat)
colnames(cul_preds)[1]="preds"

cul_plot1=ggplot(cul_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Culicidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Chironomidae
chi_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=chi_dat)

chi_preds=as.vector(predict(chi_mod, type="response"))%>%cbind(chi_dat)
colnames(chi_preds)[1]="preds"

chi_plot1=ggplot(chi_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Chironomidae")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

#Coccoidea
coc_mod=glm(TotalCatch1~ DOYs+ poly(DOYs,2) + years + # base linear terms
              DOYs*years + # timing and long-term trend interaction
              poly(DOYs, 2)* years, #curve and long-term trend interaction
            family = "poisson", data=coc_dat)

coc_preds=as.vector(predict(coc_mod, type="response"))%>%cbind(coc_dat)
colnames(coc_preds)[1]="preds"

coc_plot1=ggplot(coc_preds, aes(x=DOY, y=preds, col=as.factor(Year)))+geom_line()+
  theme_classic()+ggtitle("Coccoidea")+scale_color_viridis_d()
#  geom_point(aes(x=DOY, y=TotalCatch1), col="black")

ggarrange(bom_plot1, col_plot1, ich_plot1, aca_plot1,
          lyc_plot1, sci_plot1, pho_plot1, nym_plot1,
          mus_plot1, lin_plot1, cul_plot1, coc_plot1, chi_plot1,
          common.legend = T, legend = "right")
