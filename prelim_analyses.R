#arthropod data
arth_raw_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\arth_dat_clean.csv")
require(dplyr)
require(ggplot2)
require(WaveletComp)

arth_dat_clean=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\arth_dat_clean.csv")
mesic_dat=arth_dat_clean%>%filter(Plot.ID%in%c("Art3", "Art4"))
arid_dat=arth_dat_clean%>%filter(Plot.ID%in%c("Art5", "Art7"))

#start
mesic_dat_prop10=mesic_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

arid_dat_prop10=arid_dat%>%select(-1)%>%
  group_by(year)%>%
  mutate(cum_tot=cumsum(Total_Art),
         cum_pt=cum_tot/sum(Total_Art),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(Plot.ID, year, doy, cum_pt, cum_tot)

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

ggplot(mesic_dat_prop10)+geom_line(aes(x=year,y=doy))+
  theme_classic()+
  ylab("mean timing")+xlab("day of year")+ggtitle("mesic heath")

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

ggplot(plant_prop10)+geom_line(aes(x=year,y=DOY))+facet_wrap(~group)+
  theme_classic()+
  ylab("day of year")+xlab("year")+ggtitle("flowering onset")

#plant 10
plant_prop10=plant_grp%>%
  group_by(year, group)%>%
  mutate(cum_tot=cumsum(Flowers),
         cum_pt=cum_tot/sum(Flowers),
         PD_10=if_else(cum_pt>=0.1,"1", "0"))%>%
  filter(PD_10==1)%>%slice_head(n=1)%>%
  select(group, year, DOY, cum_pt, cum_tot)

ggplot(plant_prop10)+geom_line(aes(x=year,y=DOY))+facet_wrap(~group)+
  theme_classic()+
  ylab("day of year")+xlab("year")+ggtitle("flowering onset")

#saxifraga
sax10=plant_prop10%>%filter(group=="Saxifraga")
sax_com10=analyze.wavelet(sax10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sax_com10, color.key = "quantile",main="Saxifraga", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sax_com10)
reconstruct(sax_com10, "DOY")

#salix
sal10=plant_prop10%>%filter(group=="Salix")
sal_com10=analyze.wavelet(sal10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sal_com10, color.key = "quantile",main="Salix", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sal_com10)
reconstruct(sal_com10, "DOY")

#dryas
dry10=plant_prop10%>%filter(group=="Dryas")
dry_com10=analyze.wavelet(dry10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(dry_com10, color.key = "quantile", n.levels = 250,main="Dryas",  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(dry_com10)
reconstruct(dry_com10, "DOY")

#cassiope

cas10=plant_prop10%>%filter(group=="Cassiope")
cas_com10=analyze.wavelet(cas10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(cas_com10, color.key = "quantile",main="Cassiope", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com10)
reconstruct(cas_com10, "DOY")

#Silene
sil10=plant_prop10%>%filter(group=="Silene")
sil_com10=analyze.wavelet(sil10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sil_com10, color.key = "quantile",main="Silene", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com10)
reconstruct(cas_com10, "DOY")

#Papaver
pap10=plant_prop10%>%filter(group=="Papaver")
pap_com10=analyze.wavelet(pap10, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(pap_com10, color.key = "quantile",main="Papaver", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com10)
reconstruct(cas_com10, "DOY")

#MEAN TIMING
#plant 50
plant_prop50=plant_grp%>%
  group_by(year, group)%>%
  mutate(cum_tot=cumsum(Flowers),
         cum_pt=cum_tot/sum(Flowers),
         PD_50=if_else(cum_pt>=0.5,"1", "0"))%>%
  filter(PD_50==1)%>%slice_head(n=1)%>%
  select(group, year, DOY, cum_pt, cum_tot)

par(mfrow=c(2,3))

#saxifraga
sax50=plant_prop50%>%filter(group=="Saxifraga")
sax_com50=analyze.wavelet(sax50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sax_com50, color.key = "quantile",main="Saxifraga", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sax_com50)
reconstruct(sax_com50, "DOY")

#salix???
sal50=plant_prop50%>%filter(group=="Salix")
sal_com50=analyze.wavelet(sal50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sal_com50, color.key = "quantile",main="Salix", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sal_com50)
reconstruct(sal_com50, "DOY")

#dryas
dry50=plant_prop50%>%filter(group=="Dryas")
dry_com50=analyze.wavelet(dry50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(dry_com50, color.key = "quantile", n.levels = 250,main="Dryas",  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(dry_com50)
reconstruct(dry_com50, "DOY")

#cassiope

cas50=plant_prop50%>%filter(group=="Cassiope")
cas_com50=analyze.wavelet(cas50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(cas_com50, color.key = "quantile",main="Cassiope", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com50)
reconstruct(cas_com50, "DOY")

#Silene
sil50=plant_prop50%>%filter(group=="Silene")
sil_com50=analyze.wavelet(sil50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sil_com50, color.key = "quantile",main="Silene", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com50)
reconstruct(cas_com50, "DOY")

#Papaver
pap50=plant_prop50%>%filter(group=="Papaver")
pap_com50=analyze.wavelet(pap50, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(pap_com50, color.key = "quantile",main="Papaver", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com50)
reconstruct(cas_com50, "DOY")

#plant 90
plant_prop90=plant_grp%>%
  group_by(year, group)%>%
  mutate(cum_tot=cumsum(Flowers),
         cum_pt=cum_tot/sum(Flowers),
         PD_90=if_else(cum_pt>=0.9,"1", "0"))%>%
  filter(PD_90==1)%>%slice_head(n=1)%>%
  select(group, year, DOY, cum_pt, cum_tot)

par(mfrow=c(2,3))

#saxifraga
sax90=plant_prop90%>%filter(group=="Saxifraga")
sax_com90=analyze.wavelet(sax90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sax_com90, color.key = "quantile",main="Saxifraga", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sax_com90)
reconstruct(sax_com90, "DOY")

#saliX
sal90=plant_prop90%>%filter(group=="Salix")
sal_com90=analyze.wavelet(sal90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sal_com90, color.key = "quantile",main="Salix", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(sal_com90)
reconstruct(sal_com90, "DOY")

#dryas
dry90=plant_prop90%>%filter(group=="Dryas")
dry_com90=analyze.wavelet(dry90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(dry_com90, color.key = "quantile", n.levels = 250,main="Dryas",  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(dry_com90)
reconstruct(dry_com90, "DOY")

#cassiope

cas90=plant_prop90%>%filter(group=="Cassiope")
cas_com90=analyze.wavelet(cas90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(cas_com90, color.key = "quantile",main="Cassiope", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com90)
reconstruct(cas_com90, "DOY")

#Silene
sil90=plant_prop90%>%filter(group=="Silene")
sil_com90=analyze.wavelet(sil90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(sil_com90, color.key = "quantile",main="Silene", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com90)
reconstruct(cas_com90, "DOY")

#Papaver
pap90=plant_prop90%>%filter(group=="Papaver")
pap_com90=analyze.wavelet(pap90, "DOY", make.pval = TRUE, n.sim = 10)
wt.image(pap_com90, color.key = "quantile",main="Papaver", n.levels = 250,  legend.params = list(lab = "wavelet power levels", mar = 4.7))
wt.avg(cas_com90)
reconstruct(cas_com90, "DOY")
