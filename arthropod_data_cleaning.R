#Zack arthropod raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)

arth_raw_dat=read.table("G:\\My Drive\\SLU\\project\\GEM-datasets\\Arthropod emergence_10.17897_V285-Z265\\Arthropod emergence_10.17897_V285-Z265_data.txt",
                      sep="\t", header=T)

arth_raw_dat1=arth_raw_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date))%>%
  mutate(across(18:25, as.integer))%>%
  mutate(across(18:25, ~na_if(.,-9999)))%>%
  mutate(across(18:25))%>%
  mutate(TotalSpecimens=rowSums(.[18:25]))%>%
  select(-(5:12))%>%
  mutate(date=as.Date(paste(year, month, dia), "%Y %m %d"), DOY=yday(date))

#Exclude taxon year samples with less than lim1 specimens or less than lim2 for ichneumonidae and nymphalidea before 2007
#lim3 and 4 are the same but for after 2006 when number of traps was halved

  lim1 <- 100
  lim2 <- 50
  lim3 <- 50
  lim4 <- 25

#remove nymphalids & ichneumonids with less than 50 individuals
#before 1997 below 50
arth_dat2=arth_raw_dat1[-which(arth_raw_dat1$Family==c("Ichneumonidae", "Nymphalidae") &
                                arth_raw_dat1$TotalSpecimens< lim2 &
                                arth_raw_dat1$year<2007),]

#After 1997 below 25
arth_dat3=arth_dat2[-which(arth_raw_dat1$Family==c("Ichneumonidae", "Nymphalidae") &
                                arth_raw_dat1$TotalSpecimens< lim4 &
                                arth_raw_dat1$year>2006),]

arth_raw_dat3=anti_join(arth_raw_dat2,arth_dat3)

#remove the rest
#before 1997 below 100
arth_dat4=arth_raw_dat3[-which(arth_raw_dat3$Family==c("Ichneumonidae", "Nymphalidae") &
                                arth_raw_dat3$TotalSpecimens< lim1 &
                                arth_raw_dat3$year<2007),]

#after 1997 below 50
arth_dat5=arth_dat4[which(arth_dat4$Family==c("Ichneumonidae", "Nymphalidae") &
                                arth_dat4$TotalSpecimens< lim3 &
                                arth_dat4$year>2006),]

#Exclude bombus if still in data and 2020?
arth_raw_dat5=anti_join(arth_dat4,arth_dat5)%>%filter(!Genus=="Bombus", !year==2020)

#exclude mites for 1996 because cups were not calculated separately
arth_dat6=arth_raw_dat5[which(arth_raw_dat5$Order=="Acari" &
                              arth_raw_dat5$year==1996 &
                              arth_raw_dat5$Plot.ID=="Art1"),]

arth_raw_dat6=anti_join(arth_raw_dat5,arth_dat6)

library(ggplot2)

ggplot(data=arth_raw_dat6, aes(x=DOY, y=TotalSpecimens, col=Family))+
  geom_line()+facet_wrap(~Order)+theme_classic()

str(arth_raw_dat6)
