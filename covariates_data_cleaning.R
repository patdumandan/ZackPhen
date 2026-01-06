#Zack snow and ice raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)
require(tidyr)
require(reshape2)
require(readr)

file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data\\raw"

snow_name=paste(file_path, '\\Snow and ice cover_10.17897_BKEX-JR94_data','.txt', sep = '')
snow_raw_dat=read.csv(snow_name, header=T, sep='\t',  stringsAsFactors = F)

airtemp_name=paste(file_path, '\\Air temperature, 200cm @ 60_30min sample (°C)_10.17897_G5WS-0W04_data','.txt', sep = '')
airtemp_raw_dat=read.csv(airtemp_name, header=T, sep='\t',  stringsAsFactors = F)

#Replace -9999 or -999 with NA:
snow_raw_dat[snow_raw_dat==-9999] <- NA
snow_raw_dat[snow_raw_dat==-999] <- NA

#separate combined plots for snow data
plotd1=snow_raw_dat%>%filter(Plot%in%c("Dry2Sal7", "Dry6Pap4", "Pap2Sal5", "Sax2Sil2", "Sax1Sil1", "Sax3Sil3"))%>%
  mutate(across(c('Plot'), substr, 6, nchar(Plot)))
plotd2=snow_raw_dat%>%filter(Plot%in%c("Dry2Sal7", "Dry6Pap4", "Pap2Sal5", "Sax2Sil2", "Sax1Sil1", "Sax3Sil3"))%>%
  mutate(across(c('Plot'), substr, 1, nchar(Plot)-4))

sep_plots=rbind(plotd1, plotd2)
og_plots=snow_raw_dat%>%filter(Plot%in%c("Cas1", "Cas3", "Cas2", "Cas4",
                                         "Sal3", "Sal1", "Sal2", "Sal4",
                                         "Sil4", "Art1", "Art2", "Art3",
                                         "Dry1","Dry4" ,"Dry5", "Pap3","Art4",
                                         "Art5","Art6","Dry3", "Dry7","Dry8","Cas5",
                                         "Cas6","Art7","Sal7","Sal6"))

snow_plots=rbind(sep_plots,og_plots)%>%select(-Field_remarks, -General_remarks)%>%
  mutate(SnowCoverFraction=as.integer(SnowCoverFraction))%>%
  group_by(Date, Plot)%>%
  summarise(mean_snow=mean(SnowCoverFraction, na.rm=T))%>%
  mutate(DOY=yday(Date),
         year=year(Date))%>%
  filter(!year<1996)


snow_plant_plots=snow_plots%>%filter(!Plot%in%c("Art1", "Art2", "Art3", "Art4",
                                          "Art5", "Art6", "Art7"))

snow_arth_plots=snow_plots%>%filter(Plot%in%c("Art1", "Art2", "Art3", "Art4",
                                          "Art5", "Art6", "Art7"))

library(dplyr)
library(lubridate)

airtemp_dat_daily=airtemp_raw_dat%>%
  filter(Quality.Flag != "missing")%>%
  mutate(Date = as.Date(Date),
         DOY  = yday(Date),
         month= month(Date),
         year = year(Date))%>%        # ensure daily grouping if Date is datetime
  group_by(Date) %>%
  summarise(min_temp  = min(Temp, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            max_temp  = max(Temp, na.rm = TRUE),
            .groups = "drop")%>%
  filter(year >= 1996)%>%
  unique()

colnames(airtemp_dat)[3]="Temp"

airtemp_dat_spring=airtemp_raw_dat%>%mutate(year = year(Date),
                                            month= month(Date),
                                            DOY  = yday(Date))%>%
  filter(month%in%c(4:5))%>%group_by(year)%>%
  summarise(min_temp  = min(Temp, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            max_temp  = max(Temp, na.rm = TRUE),
            .groups = "drop")%>%
  filter(year >= 1996)%>%
  unique()

airtemp_dat_summer=airtemp_raw_dat%>%mutate(year = year(Date),
                                            month= month(Date),
                                            DOY  = yday(Date))%>%
  filter(month%in%c(6:8))%>%group_by(year)%>%
  summarise(min_temp  = min(Temp, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            max_temp  = max(Temp, na.rm = TRUE),
            .groups = "drop")%>%
  filter(year >= 1996)%>%
  unique()
