#plants
require(tidyr)
require(dplyr)
require(lubridate)

###cassiope

cas_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Cassiope phenology_10.17897_X0MY-K003_data.txt",
                   sep="\t", header=T)
#str(cas_dat)

cas_dat$Date=as.POSIXct(cas_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
cas_dat_raw=cas_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Cassiope")%>%
  mutate(DOY=yday(Date))%>%
 filter(!year<1997)%>%
  select(-Buds, -Field_remarks, -General_remarks, -Date)

###dryas

dry_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Dryas phenology_10.17897_JSQ7-6355_data.txt",
                 sep="\t", header=T)
#str(cas_dat)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
dry_dat_raw=dry_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Dryas")%>%
  mutate(DOY=yday(Date))%>%
  filter(!year<1997)%>%
  select(-Buds, -Field_remarks, -General_remarks, -Date)

###Papaver
pap_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Papaver phenology_10.17897_NK32-H804_data.txt",
                 sep="\t", header=T)

pap_dat$Date=as.POSIXct(pap_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
pap_dat_raw=pap_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Papaver")%>%
  mutate(DOY=yday(Date))%>%
  filter(!year<1997)%>%
  select(-Buds, -Field_remarks, -General_remarks, -Date)

###Salix
sal_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Salix phenology_10.17897_NS7W-JT18_data.txt",
                 sep="\t", header=T)

sal_dat$Date=as.POSIXct(sal_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sal_dat_raw=sal_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Salix")%>%
  mutate(DOY=yday(Date))%>%
  filter(!year<1997)%>%
  select(-Buds, -Field_remarks, -General_remarks, -Date)

###Saxifraga
sax_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Saxifraga phenology_10.17897_YXH1-ZB25_data.txt",
                 sep="\t", header=T)

sax_dat$Date=as.POSIXct(sax_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sax_dat_raw=sax_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Saxifraga")%>%
  mutate(DOY=yday(Date))%>%
  filter(!year<1997)%>%
  select(-Buds, -General_remarks, -Date)

###Silene
sil_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Silene phenology_10.17897_6GVG-QH42_data.txt",
                 sep="\t", header=T)

sil_dat$Date=as.POSIXct(sil_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sil_dat_raw=sil_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Silene")%>%
  mutate(DOY=yday(Date))%>%
  filter(!year<1997)%>%
  select(-Buds,  -Date)


plant_dat=vec_rbind(sil_dat_raw, sax_dat_raw, sal_dat_raw, pap_dat_raw, dry_dat_raw, cas_dat_raw)%>%
  select(Plot, Section, Flowers, year, month, dia, DOY, species)%>%
  rename(group=species)

#(C): Replace -9999 or -999 with NA:
plant_dat[plant_dat==-9999] <- NA
plant_dat[plant_dat==-999] <- NA

write.csv(plant_dat, "ZAC_plant_raw.csv")
