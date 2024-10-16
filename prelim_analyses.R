#arthropod data
arth_raw_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\arth_dat_clean.csv")
require(dplyr)
require(ggplot2)
require(tidyr)
require(WaveletComp)

arth_dat_clean=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\arth_dat_clean.csv")
arth_dat_clean=arth_dat_clean%>%filter(!Order%in%c("sp.", "unidentified"))

#all pooled
arth_dat_grp=arth_dat_clean%>%
  group_by(year, Order, doy)%>%
  mutate(wk_tot=sum(Total_Art, na.rm=T))%>%
  group_by(Order, year)%>%
  mutate(yr_tot=sum(Total_Art, na.rm=T))%>%
  filter(!year<1996)

arth_dat_prop10=arth_dat_grp%>%
  group_by(year, Order)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot, Order)

arth_dat_prop50=arth_dat_grp%>%
  group_by(year, Order)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_50=if_else(cum_pt>=0.5,"1", "0"))%>%
  filter(PD_50==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot, Order)

arth_dat_prop90=arth_dat_grp%>%
  group_by(year, Order)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_90=if_else(cum_pt>=0.9,"1", "0"))%>%
  filter(PD_90==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot, Order)

par(mfrow=c(3,1))
ggplot(arth_dat_prop10, aes(x=year,y=doy))+geom_line(aes(x=year,y=doy, col=Order))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("onset")

ggplot(arth_dat_prop50, aes(x=year,y=doy))+geom_line(aes(x=year,y=doy, col=Order))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("mid-season")

ggplot(arth_dat_prop90, aes(x=year,y=doy))+geom_line(aes(x=year,y=doy, col=Order))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("end")

#by environment
mesic_dat=arth_dat_clean%>%filter(Plot.ID%in%c("Art3", "Art4"))
arid_dat=arth_dat_clean%>%filter(Plot.ID%in%c("Art5", "Art7"))

#start
mesic_dat_grp=mesic_dat%>%filter(!Genus%in%c("sp.", "unidentified"))%>%
  group_by(year, Genus, doy)%>%
  mutate(wk_tot=sum(Total_Art, na.rm=T))%>%
  group_by(Genus, year)%>%
  mutate(yr_tot=sum(Total_Art, na.rm=T))%>%
  filter(!year<1996)

arid_dat_grp=mesic_dat%>%filter(!Genus%in%c("sp.", "unidentified"))%>%
  group_by(year, Genus, doy)%>%
  mutate(wk_tot=sum(Total_Art, na.rm=T))%>%
  group_by(Genus, year)%>%
  mutate(yr_tot=sum(Total_Art, na.rm=T))%>%
  filter(!year<1996)

mesic_dat_prop10=mesic_dat_grp%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot, Genus)

arid_dat_prop10=arid_dat_grp%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot, Genus)

#mean timing
mesic_dat_prop50=mesic_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_50=if_else(cum_pt>=0.5,"1", "0"))%>%
  filter(PD_50==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

arid_dat_prop50=arid_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_50=if_else(cum_pt>=0.5,"1", "0"))%>%
  filter(PD_50==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

#end
mesic_dat_prop90=mesic_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_90=if_else(cum_pt>=0.9,"1", "0"))%>%
  filter(PD_90==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

arid_dat_prop90=arid_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_90=if_else(cum_pt>=0.9,"1", "0"))%>%
  filter(PD_90==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

ggplot(mesic_dat_prop10)+geom_line(aes(x=year,y=doy, col=Genus))+
  theme_classic()+
  ylab("onset")+xlab("day of year")+ggtitle("mesic heath")

ggplot(arid_dat_prop10)+geom_line(aes(x=year,y=doy))+
  theme_classic()+
  ylab("mean timing")+xlab("day of year")+ggtitle("arid heath")

wv_com90=analyze.wavelet(mesic_dat_prop90, "doy", make.pval = TRUE, n.sim = 10)
wt.image(wv_com90, color.key = "quantile", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(wv_com)
reconstruct(wv_com10, "doy")

wv2_com=analyze.wavelet(arid_dat_prop90, "doy",make.pval = TRUE, n.sim = 10)
wt.image(wv2_com, color.key = "quantile", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(wv2_com)
reconstruct(wv2_com, "doy")

#plant data

plant_dat=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\data\\ZAC_plant_raw.csv")

plant_grp=plant_dat%>%
  group_by(year, group, DOY)%>%
  mutate(wk_tot=sum(Flowers, na.rm=T))%>%
  group_by(group, year)%>%
  mutate(yr_tot=sum(Flowers, na.rm=T))%>%
  filter(!year<1996)

#plant 10
plant_prop10=plant_grp%>%
  group_by(year, group)%>%
  mutate(cum_tot=cumsum(Flowers),
         cum_pt=cum_tot/sum(Flowers),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(group, year, DOY, cum_pt, cum_tot)

#plant 10
plant_prop10=plant_grp%>%
  group_by(year, group)%>%
  mutate(cum_tot=cumsum(Flowers),
         cum_pt=cum_tot/sum(Flowers),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(group, year, DOY, cum_pt, cum_tot)

ggplot(plant_prop10, aes(x=year,y=DOY))+geom_line(aes(x=year,y=DOY, col=group))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("flowering onset")

ggplot(plant_prop50, aes(x=year,y=DOY))+geom_line(aes(x=year,y=DOY, col=group))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("flowering peak")

ggplot(plant_prop90, aes(x=year,y=DOY))+geom_line(aes(x=year,y=DOY, col=group))+
  theme_classic()+geom_smooth(method="gam")+
  ylab("day of year")+xlab("year")+ggtitle("flowering end")


#plant data####
phne_dat=Phenology_Plants_Zackenberg_1996_2020%>%group_by(Species, Year)%>%
  summarise(mean_doy=mean(DOY_50_flowering))

phne_dat2=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_plant_phenology_1996-2023.csv")%>%
  group_by(Species, Year)%>%summarise(mean_doy=mean(DOY_50_flowering))

##50% bloom date####
#compare results when using 1996-2020 and 1996-2023

#Dryas####
dry1=phne_dat%>%filter(Species=="Dryas")
dry2=phne_dat2%>%filter(Species=="Dryas")

dry_com1=analyze.wavelet(dry1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(dry_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(dry1$Year), by = 1), labels = unique(dry1$Year)))
wt.image(dry_com1, color.key = "quantile",main="Dryas (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry1$Year), by = 1), labels = unique(dry1$Year)))

dry_com2=analyze.wavelet(dry2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(dry_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.image(dry_com2, color.key = "quantile",main="Dryas (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))

#Salix####
sal1=phne_dat%>%filter(Species=="Salix")
sal2=phne_dat2%>%filter(Species=="Salix")

sal_com1=analyze.wavelet(sal1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sal_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sal1$Year), by = 1), labels = unique(sal1$Year)))
wt.image(sal_com1, color.key = "quantile",main="Salix (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal1$Year), by = 1), labels = unique(sal1$Year)))

sal_com2=analyze.wavelet(sal2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sal_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.image(sal_com2, color.key = "quantile",main="Salix (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))


#Cassiope####
cas1=phne_dat%>%filter(Species=="Cassiope")
cas2=phne_dat2%>%filter(Species=="Cassiope")

cas_com1=analyze.wavelet(cas1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(cas_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(cas1$Year), by = 1), labels = unique(cas1$Year)))
wt.image(cas_com1, color.key = "quantile",main="Cassiope (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas1$Year), by = 1), labels = unique(cas1$Year)))

cas_com2=analyze.wavelet(cas2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(cas_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.image(cas_com2, color.key = "quantile",main="Cassiope (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))

#Silene####
sil1=phne_dat%>%filter(Species=="Silene")
sil2=phne_dat2%>%filter(Species=="Silene")

sil_com1=analyze.wavelet(sil1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sil_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sil1$Year), by = 1), labels = unique(sil1$Year)))
wt.image(sil_com1, color.key = "quantile",main="Silene (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil1$Year), by = 1), labels = unique(sil1$Year)))

sil_com2=analyze.wavelet(sil2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sil_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.image(sil_com2, color.key = "quantile",main="Silene (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))

#Saxifraga####

sax1=phne_dat%>%filter(Species=="Saxifraga")
sax2=phne_dat2%>%filter(Species=="Saxifraga")

sax_com1=analyze.wavelet(sax1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sax_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sax1$Year), by = 1), labels = unique(sax1$Year)))
wt.image(sax_com1, color.key = "quantile",main="Saxifraga (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax1$Year), by = 1), labels = unique(sax1$Year)))

sax_com2=analyze.wavelet(sax2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(sax_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.image(sil_com2, color.key = "quantile",main="Saxifraga (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))

#Papaver####

pap1=phne_dat%>%filter(Species=="Papaver")
pap2=phne_dat2%>%filter(Species=="Papaver")

pap_com1=analyze.wavelet(pap1, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(pap_com1, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(pap1$Year), by = 1), labels = unique(pap1$Year)))
wt.image(pap_com1, color.key = "quantile",main="Papaver (1996-2020)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap1$Year), by = 1), labels = unique(pap1$Year)))

pap_com2=analyze.wavelet(pap2, "mean_doy", make.pval = TRUE, n.sim = 10)
reconstruct(pap_com2, "mean_doy", show.legend = F,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.image(pap_com2, color.key = "quantile",main="Papaver (1996-2023)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))

phne_dat10=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_plant_phenology_1996-2023.csv")%>%
  group_by(Species, Year)%>%summarise(mean_doy10=mean(DOY_10_flowering))

##10% bloom date####
#Dryas####
dry2=phne_dat10%>%filter(Species=="Dryas")

dry_com2=analyze.wavelet(dry2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(dry_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.image(dry_com2, color.key = "quantile",main="Dryas (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))

#Salix####
sal2=phne_dat10%>%filter(Species=="Salix")

sal_com2=analyze.wavelet(sal2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(sal_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.image(sal_com2, color.key = "quantile",main="Salix (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))


#Cassiope####
cas2=phne_dat10%>%filter(Species=="Cassiope")

cas_com2=analyze.wavelet(cas2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(cas_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.image(cas_com2, color.key = "quantile",main="Cassiope (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))

#Silene####
sil2=phne_dat10%>%filter(Species=="Silene")

sil_com2=analyze.wavelet(sil2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(sil_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.image(sil_com2, color.key = "quantile",main="Silene (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))

#Saxifraga####
sax2=phne_dat10%>%filter(Species=="Saxifraga")

sax_com2=analyze.wavelet(sax2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(sax_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.image(sax_com2, color.key = "quantile",main="Saxifraga (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))

#Papaver####
pap2=phne_dat10%>%filter(Species=="Papaver")

pap_com2=analyze.wavelet(pap2, "mean_doy10", make.pval = TRUE, n.sim = 10)
reconstruct(pap_com2, "mean_doy10", show.legend = F,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.image(pap_com2, color.key = "quantile",main="Papaver (10% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))

##90% bloom date####
phne_dat90=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_plant_phenology_1996-2023.csv")%>%
  group_by(Species, Year)%>%summarise(mean_doy90=mean(DOY_90_flowering))

#Dryas####
dry2=phne_dat90%>%filter(Species=="Dryas")

dry_com2=analyze.wavelet(dry2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(dry_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.image(dry_com2, color.key = "quantile",main="Dryas (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))

#Salix####
sal2=phne_dat90%>%filter(Species=="Salix")

sal_com2=analyze.wavelet(sal2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(sal_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.image(sal_com2, color.key = "quantile",main="Salix (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))

#Cassiope####
cas2=phne_dat90%>%filter(Species=="Cassiope")

cas_com2=analyze.wavelet(cas2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(cas_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.image(cas_com2, color.key = "quantile",main="Cassiope (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))

#Silene####
sil2=phne_dat90%>%filter(Species=="Silene")

sil_com2=analyze.wavelet(sil2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(sil_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.image(sil_com2, color.key = "quantile",main="Silene (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))

#Saxifraga####
sax2=phne_dat90%>%filter(Species=="Saxifraga")

sax_com2=analyze.wavelet(sax2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(sax_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.image(sax_com2, color.key = "quantile",main="Saxifraga (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))

#Papaver####
pap2=phne_dat90%>%filter(Species=="Papaver")

pap_com2=analyze.wavelet(pap2, "mean_doy90", make.pval = TRUE, n.sim = 10)
reconstruct(pap_com2, "mean_doy90", show.legend = F,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.image(pap_com2, color.key = "quantile",main="Papaver (90% bloom)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))

##duration####
plant_timing=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_plant_phenology_1996-2023.csv")

plant_timings=plant_timing%>%
  mutate(duration=DOY_90_flowering-DOY_10_flowering)%>%
  group_by(Species, Year)%>%summarise(mean_duration=mean(duration))

par(mfrow=c(2,1))

#Dryas####
dry2=plant_timings%>%filter(Species=="Dryas")

dry_com2=analyze.wavelet(dry2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(dry_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))
wt.image(dry_com2, color.key = "quantile",main="Dryas (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(dry2$Year), by = 1), labels = unique(dry2$Year)))

#Salix####
sal2=plant_timings%>%filter(Species=="Salix")

sal_com2=analyze.wavelet(sal2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(sal_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))
wt.image(sal_com2, color.key = "quantile",main="Salix (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sal2$Year), by = 1), labels = unique(sal2$Year)))

#Silene####
sil2=plant_timings%>%filter(Species=="Silene")

sil_com2=analyze.wavelet(sil2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(sil_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
wt.image(sil_com2, color.key = "quantile",main="Silene (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sil2$Year), by = 1), labels = unique(sil2$Year)))
#Cassiope####
cas2=plant_timings%>%filter(Species=="Cassiope")

cas_com2=analyze.wavelet(cas2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(cas_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))
wt.image(cas_com2, color.key = "quantile",main="Cassiope (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(cas2$Year), by = 1), labels = unique(cas2$Year)))

#Saxifraga####
sax2=plant_timings%>%filter(Species=="Saxifraga")

sax_com2=analyze.wavelet(sax2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(sax_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))
wt.image(sax_com2, color.key = "quantile",main="Saxifraga (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(sax2$Year), by = 1), labels = unique(sax2$Year)))

#Papaver####
pap2=plant_timings%>%filter(Species=="Papaver")

pap_com2=analyze.wavelet(pap2, "mean_duration", make.pval = TRUE, n.sim = 10)
reconstruct(pap_com2, "mean_duration", show.legend = F,
            spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))
wt.image(pap_com2, color.key = "quantile",main="Papaver (duration)",
         n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7),
         spec.time.axis =list(at = seq(1, length(pap2$Year), by = 1), labels = unique(pap2$Year)))

#focal species####
phen_dat_all=read.csv("I:\\My Drive\\SLU\\phenology-project\\ZackPhen\\ZAC_phenology_metrics_1996-2023.csv")
foc_dat=phen_dat_all%>%filter(Species%in%c("Dryas", "Salix"))

require(ggplot2)

ggplot(foc_dat, aes(x=DOY, col=Year, fill=Year))+
  geom_histogram(position="dodge", binwidth = 50)+
  facet_wrap(vars(Species, metric))+theme_classic()

#moving window analyses####
#Dryas####
dry2_dat=foc_dat%>%filter(Species=="Dryas", Plot=="Dry2", metric=="50")

dry2_dat_ro <-
  rolling_origin(
    data       = dry2_dat, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = TRUE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat_ro$model=map(dry2_dat_ro$splits, rolling_mod)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat_ro$slope=dry2_dat_ro$model%>%map(coef)%>%map_dbl(2)
dry2_dat_ro$pval=map(dry2_dat_ro$model, get_pvalue)


require(ggpubr)

ds=dry2_dat_ro%>%select(slope, id, pval)%>%mutate(metric="50")
ggdensity(ds, x="slope", add="median")

#10% bloom####
dry2_dat10=foc_dat%>%filter(Species=="Dryas", Plot=="Dry2", metric=="10")

dry2_dat10_ro <-
  rolling_origin(
    data       = dry2_dat10, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = TRUE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat10_ro$model=map(dry2_dat10_ro$splits, rolling_mod)

dry2_dat10_ro$slope=dry2_dat10_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat10_ro$pval=map(dry2_dat10_ro$model, get_pvalue)


require(ggpubr)

ds10=dry2_dat10_ro%>%select(slope, id, pval)%>%mutate(metric="10")
ggdensity(ds10, x="slope", add="median")

dry2_dat90=foc_dat%>%filter(Species=="Dryas", Plot=="Dry2", metric=="90")

dry2_dat90_ro <-
  rolling_origin(
    data       = dry2_dat90, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = TRUE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat90_ro$model=map(dry2_dat90_ro$splits, rolling_mod)

dry2_dat90_ro$slope=dry2_dat90_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat90_ro$pval=map(dry2_dat90_ro$model, get_pvalue)


require(ggpubr)

ds90=dry2_dat90_ro%>%select(slope, id, pval)%>%mutate(metric="90")
ggdensity(ds90, x="slope", add="median")


dry2_dens=bind_rows(ds, ds10, ds90)
ggdensity(dry2_dens, x="slope", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0)+ggtitle("Dryas (Dry 2)")
#moving window analyses####

#MEAN ACROSS ALL PLOTS

dry2_dat=foc_dat%>%filter(Species=="Dryas", metric=="50")

dry2_dat=foc_dat%>%filter(Species=="Dryas", metric=="50")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

dry2_dat=as.data.frame(dry2_dat)

dry2_dat_ro <-
  rolling_origin(
    data       = dry2_dat, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat_ro$model=map(dry2_dat_ro$splits, rolling_mod)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat_ro$slope=dry2_dat_ro$model%>%map(coef)%>%map_dbl(2)
dry2_dat_ro$pval=map(dry2_dat_ro$model, get_pvalue)


require(ggpubr)

ds=dry2_dat_ro%>%select(slope, id, pval)%>%mutate(metric="50")
ggdensity(ds, x="slope", add="median")

#10% bloom####
dry2_dat10=foc_dat%>%filter(Species=="Dryas", metric=="10")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

dry2_dat10=as.data.frame(dry2_dat10)

dry2_dat10_ro <-
  rolling_origin(
    data       = dry2_dat10, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat10_ro$model=map(dry2_dat10_ro$splits, rolling_mod)

dry2_dat10_ro$slope=dry2_dat10_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat10_ro$pval=map(dry2_dat10_ro$model, get_pvalue)

require(ggpubr)

ds10=dry2_dat10_ro%>%select(slope, id, pval)%>%mutate(metric="10")
ggdensity(ds10, x="slope", add="median")

dry2_dat90=foc_dat%>%filter(Species=="Dryas", Plot=="Dry2", metric=="90")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

dry2_dat90=as.data.frame(dry2_dat90)

dry2_dat90_ro <-
  rolling_origin(
    data       = dry2_dat90, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

dry2_dat90_ro$model=map(dry2_dat90_ro$splits, rolling_mod)

dry2_dat90_ro$slope=dry2_dat90_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

dry2_dat90_ro$pval=map(dry2_dat90_ro$model, get_pvalue)

require(ggpubr)

ds90=dry2_dat90_ro%>%select(slope, id, pval)%>%mutate(metric="90")
ggdensity(ds90, x="slope", add="median")

dry2_dens=bind_rows(ds, ds10, ds90)

dry2_dens$pval=as.vector(as.numeric(dry2_dens$pval))

ggdensity(dry2_dens, x="slope", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0)+ggtitle("Dryas")

ggdensity(dry2_dens, x="pval", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0.05)+ggtitle("Dryas")

##Salix####
###50% bloom####
sa=foc_dat%>%filter(Species=="Salix")

sal1_dat=foc_dat%>%filter(Species=="Salix", Plot=="Sal1", metric=="50")

sal1_dat_ro <-
  rolling_origin(
    data       = sal1_dat, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

sal1_dat_ro$model=map(sal1_dat_ro$splits, rolling_mod)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

sal1_dat_ro$slope=sal1_dat_ro$model%>%map(coef)%>%map_dbl(2)
sal1_dat_ro$pval=map(sal1_dat_ro$model, get_pvalue)


require(ggpubr)

ds=sal1_dat_ro%>%select(slope, id, pval)%>%mutate(metric="50")
ggdensity(ds, x="slope", add="median")

##10% bloom####
sal1_dat10=foc_dat%>%filter(Species=="Salix", Plot=="Sal1", metric=="10")

sal1_dat10_ro <-
  rolling_origin(
    data       = sal1_dat10, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

sal1_dat10_ro$model=map(sal1_dat10_ro$splits, rolling_mod)

sal1_dat10_ro$slope=sal1_dat10_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

sal1_dat10_ro$pval=map(sal1_dat10_ro$model, get_pvalue)


require(ggpubr)

ds10=sal1_dat10_ro%>%select(slope, id, pval)%>%mutate(metric="10")
ggdensity(ds10, x="slope", add="median")

###90% bloom####
sal1_dat90=foc_dat%>%filter(Species=="Salix", Plot=="Sal1", metric=="90")

sal1_dat90_ro <-
  rolling_origin(
    data       = sal1_dat90, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

sal1_dat90_ro$model=map(sal1_dat90_ro$splits, rolling_mod)

sal1_dat90_ro$slope=sal1_dat90_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

sal1_dat90_ro$pval=map(sal1_dat90_ro$model, get_pvalue)


require(ggpubr)

ds90=sal1_dat90_ro%>%select(slope, id, pval)%>%mutate(metric="90")
ggdensity(ds90, x="slope", add="median")


sal1_dens=bind_rows(ds, ds10, ds90)
sal1_dens$pval=as.vector(as.numeric(sal1_dens$pval))

ggdensity(sal1_dens, x="slope", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0)+ggtitle("Salix (Sal1)")

ggdensity(sal1_dens, x="pval", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0.05)+ggtitle("Salix (Sal1)")

#MEAN ACROSS ALL PLOTS####
###50% bloom####
sal1_dat1=foc_dat%>%filter(Species=="Salix", metric=="50")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

sal1_dat1=as.data.frame(sal1_dat1)

sal1_dat1_ro <-
  rolling_origin(
    data       = sal1_dat1, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

sal1_dat1_ro$model=map(sal1_dat1_ro$splits, rolling_mod)

sal1_dat1_ro$slope=sal1_dat1_ro$model%>%map(coef)%>%map_dbl(2)
sal1_dat1_ro$pval=map(sal1_dat1_ro$model, get_pvalue)

###10% bloom####
sal1_dat10=foc_dat%>%filter(Species=="Salix", metric=="10")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

sal1_dat10=as.data.frame(sal1_dat10)

sal1_dat10_ro <-
  rolling_origin(
    data       = sal1_dat10, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

rolling_mod=function(split) {

  analysis_set= analysis(split) #get dataframe

  fit_model= lm(analysis_set[,"DOY"]~analysis_set[,"Year"])
}

sal1_dat10_ro$model=map(sal1_dat10_ro$splits, rolling_mod)

sal1_dat10_ro$slope=sal1_dat10_ro$model%>%map(coef)%>%map_dbl(2)

get_pvalue=function(model) {

  coefficients <- summary(model)$coefficients
  pval=coefficients[8]
}

sal1_dat10_ro$pval=map(sal1_dat10_ro$model, get_pvalue)
s10=sal1_dat10_ro%>%select(slope, id, pval)%>%mutate(metric="10")

###90% bloom####
sal1_dat90=foc_dat%>%filter(Species=="Salix", Plot=="Sal1", metric=="90")%>%
  group_by(Year, metric)%>%
  summarise(DOY=as.numeric(mean(DOY)))

sal1_dat90=as.data.frame(sal1_dat90)

sal1_dat90_ro <-
  rolling_origin(
    data       = sal1_dat90, #all PB control data (1999-2009)
    initial    = 10, #samples used for modelling (training)
    assess     = 3, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

sal1_dat90_ro$model=map(sal1_dat90_ro$splits, rolling_mod)

sal1_dat90_ro$slope=sal1_dat90_ro$model%>%map(coef)%>%map_dbl(2)

sal1_dat90_ro$pval=map(sal1_dat90_ro$model, get_pvalue)

s90=sal1_dat90_ro%>%select(slope, id, pval)%>%mutate(metric="90")

sal1_dens=bind_rows(s, s10, s90)

sal1_dens$pval=as.vector(as.numeric(sal1_dens$pval))

ggdensity(sal1_dens, x="slope", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0)+ggtitle("Salix")

ggdensity(sal1_dens, x="pval", add="median", size=1, col="metric", fill="metric")+
  geom_vline(xintercept = 0.05)+ggtitle("Salix")


