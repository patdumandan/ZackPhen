#plants
require(tidyr)
require(dplyr)
require(lubridate)

#data####
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data\\raw"

###cassiope####
cas_name=paste(file_path, '\\Cassiope phenology_10.17897_X0MY-K003_data','.txt', sep = '')
cas_dat=read.csv(cas_name, header=T, sep='\t',  stringsAsFactors = F)

cas_dat$Date=as.POSIXct(cas_dat$Date, tz="GMT", format = "%Y-%m-%d")

cas_dat_raw=cas_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Cassiope")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)%>%
  filter(Plot%in%c("Cas1", "Cas2", "Cas3", "Cas4"))

###dryas####
dry_name=paste(file_path, '\\Dryas phenology_10.17897_JSQ7-6355_data','.txt', sep = '')
dry_dat=read.csv(dry_name, header=T, sep='\t',  stringsAsFactors = F)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")

dry_dat_raw=dry_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Dryas")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)%>%
  filter(Plot%in%c("Dry1", "Dry2", "Dry3", "Dry4", "Dry5", "Dry6"))

###papaver####
pap_name=paste(file_path, '\\Papaver phenology_10.17897_NK32-H804_data','.txt', sep = '')
pap_dat=read.csv(pap_name, header=T, sep='\t',  stringsAsFactors = F)

pap_dat$Date=as.POSIXct(pap_dat$Date, tz="GMT", format = "%Y-%m-%d")

pap_dat_raw=pap_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Papaver")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)%>%
  filter(Plot%in%c("Pap1", "Pap2", "Pap3", "Pap4"))

###salix####
sal_name=paste(file_path, '\\Salix phenology_10.17897_NS7W-JT18_data','.txt', sep = '')
sal_dat=read.csv(sal_name, header=T, sep='\t',  stringsAsFactors = F)

sal_dat$Date=as.POSIXct(sal_dat$Date, tz="GMT", format = "%Y-%m-%d")

sal_dat_raw=sal_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Salix")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)%>%
  filter(Plot%in%c("Sal1", "Sal2", "Sal3", "Sal4", "Sal5", "Sal6", "Sal7"))

###saxifraga####
sax_name=paste(file_path, '\\Saxifraga phenology_10.17897_YXH1-ZB25_data','.txt', sep = '')
sax_dat=read.csv(sax_name, header=T, sep='\t',  stringsAsFactors = F)

sax_dat$Date=as.POSIXct(sax_dat$Date, tz="GMT", format = "%Y-%m-%d")

sax_dat_raw=sax_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Saxifraga")%>%
  mutate(DOY=yday(Date))%>%
  select(-General_remarks, -Date)%>%
  filter(Plot%in%c("Sax1", "Sax2", "Sax3"))

###silene####
sil_name=paste(file_path, '\\Silene phenology_10.17897_6GVG-QH42_data','.txt', sep = '')
sil_dat=read.csv(sil_name, header=T, sep='\t',  stringsAsFactors = F)

sil_dat$Date=as.POSIXct(sil_dat$Date, tz="GMT", format = "%Y-%m-%d")

sil_dat_raw=sil_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Silene")%>%
  mutate(DOY=yday(Date))%>%
  select(-General_Remarks, -Field_Remarks, -Date)%>%
  filter(Plot%in%c("Sil1", "Sil2", "Sil3", "Sil4"))

#all plants####
plant_dat_raw=vctrs::vec_rbind(sil_dat_raw, sax_dat_raw, sal_dat_raw, pap_dat_raw, dry_dat_raw, cas_dat_raw)

#(C): Replace -9999 or -999 with NA:
plant_dat_raw[plant_dat_raw==-9999] <- 0
plant_dat_raw[plant_dat_raw==-999] <- 0

plant_dat_clean=plant_dat_raw%>%
  mutate(tot_flower=rowSums(across(c(Male_flowers, Female_flowers, Flowers, Open)), na.rm=T),
         tot_bud=Buds,
         tot_senescent=rowSums(across(c(Senescent, Seed_hairs)), na.rm=T))%>%
  select(Plot, Section, year, month, dia, species, DOY, tot_flower, tot_bud, tot_senescent)%>%
  filter(!year<1997, month%in%c(6:9))

#data reshaping####
plant_dat=plant_dat_clean%>%
  group_by(Plot,species, year, DOY)%>%summarise(tot_buds=sum(tot_bud),
                                                tot_flowers=sum(tot_flower),
                                                tot_sen=sum(tot_senescent))

plant_datA <- plant_dat %>%
  mutate(
    tot_F = tot_flowers,
    tot_NF = tot_buds + tot_sen
  ) %>%
  group_by(species) %>%
  mutate(plot_id = dense_rank(Plot)) %>%
  ungroup()


#standardize variables
plant_datA$yearc = as.factor(plant_datA$year) # - mean(plant_datA$year))/sd(plant_datA$year)

#plant_datA$years = (plant_datA$year - mean(plant_datA$year))/sd(plant_datA$year)
plant_datA$DOYs = scale(plant_datA$DOY, center = TRUE, scale = TRUE)[,1]
#polynomial term for DOY, no need to standardize to avoid distortion of polynomial relationship,
#and interfere interpretability of interaction terms
plant_datA$DOYsqs=plant_datA$DOYs^2

