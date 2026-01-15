#Zack snow and ice raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)
require(tidyr)
require(reshape2)
require(readr)

file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data\\raw"


snow_name=paste(file_path, '\\Spring_snow_cover_10.17897_M7Y2-TK96_data','.txt', sep = '')
snow_raw_dat=read.csv(snow_name, header=T, sep='\t',  stringsAsFactors = F)

airtemp_name=paste(file_path, '\\Air temperature, 200cm @ 60_30min sample (°C)_10.17897_G5WS-0W04_data','.txt', sep = '')
airtemp_raw_dat=read.csv(airtemp_name, header=T, sep='\t',  stringsAsFactors = F)

colnames(airtemp_raw_dat)[3]="Temp"
airtemp_raw_dat$Date=as.Date(airtemp_raw_dat$Date)

airtemp_raw_dat[airtemp_raw_dat==-9999]=NA

#snow####
snow_ave=snow_raw_dat%>%
  group_by(Year)%>%
  summarise(Jun10_cover=mean(PercentSnowCover))%>%
  filter(!Year<1996)

#air temp####
airtemp_dat_daily= airtemp_raw_dat %>%
  filter(Quality.Flag != "missing") %>%
  mutate(year = year(Date),
         month= month(Date),
         DOY  = yday(Date)) %>%
  group_by(Date, year) %>%
  summarise(
    min_temp  = min(Temp, na.rm = TRUE),
    mean_temp = mean(Temp, na.rm = TRUE),
    max_temp  = max(Temp, na.rm = TRUE),
    .groups = "drop")%>%
  filter(year >= 1996)

airtemp_dat_spring=airtemp_raw_dat%>%mutate(year = year(Date),
                                            month= month(Date),
                                            DOY  = yday(Date))%>%
  filter(month%in%c(4:5))%>%group_by(year)%>%
  summarise(spring_min_temp  = min(Temp, na.rm = TRUE),
            spring_mean_temp = mean(Temp, na.rm = TRUE),
            spring_max_temp  = max(Temp, na.rm = TRUE),
            .groups = "drop")%>%
  filter(year >= 1996)%>%
  unique()

airtemp_dat_summer=airtemp_raw_dat%>%mutate(year = year(Date),
                                            month= month(Date),
                                            DOY  = yday(Date))%>%
  filter(month%in%c(6:7))%>%group_by(year)%>%
  summarise(summer_min_temp  = min(Temp, na.rm = TRUE),
            summer_mean_temp = mean(Temp, na.rm = TRUE),
            summer_max_temp  = max(Temp, na.rm = TRUE),
            .groups = "drop")%>%
  filter(year >= 1996)%>%
  unique()

airtemp_dat_apr_jul=full_join(airtemp_dat_spring,airtemp_dat_summer)
