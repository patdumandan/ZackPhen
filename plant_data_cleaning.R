#plants
require(tidyr)
require(dplyr)
require(lubridate)

#load data####
###cassiope

cas_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Cassiope phenology_10.17897_X0MY-K003_data.txt",
                   sep="\t", header=T)
#str(cas_dat)

cas_dat$Date=as.POSIXct(cas_dat$Date, tz="GMT", format = "%Y-%m-%d")
#remove before 1997 because totals were across all sections
cas_dat_raw=cas_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Cassiope")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)

###dryas

dry_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Dryas phenology_10.17897_JSQ7-6355_data.txt",
                 sep="\t", header=T)
#str(cas_dat)

dry_dat$Date=as.POSIXct(dry_dat$Date, tz="GMT", format = "%Y-%m-%d")
#remove before 1997 because totals were across all sections
dry_dat_raw=dry_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Dryas")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)

###Papaver
pap_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Papaver phenology_10.17897_NK32-H804_data.txt",
                 sep="\t", header=T)

pap_dat$Date=as.POSIXct(pap_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
pap_dat_raw=pap_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Papaver")%>%
  mutate(DOY=yday(Date))%>%
  select(-Field_remarks, -General_remarks, -Date)

###Salix
sal_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Salix phenology_10.17897_NS7W-JT18_data.txt",
                 sep="\t", header=T)

sal_dat$Date=as.POSIXct(sal_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sal_dat_raw=sal_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Salix")%>%
  mutate(DOY=yday(Date))%>%rename(Flowers=Female_flowers)%>%
  select(-Field_remarks, -General_remarks, -Date, -Male_flowers)

###Saxifraga
sax_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Saxifraga phenology_10.17897_YXH1-ZB25_data.txt",
                 sep="\t", header=T)

sax_dat$Date=as.POSIXct(sax_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sax_dat_raw=sax_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Saxifraga")%>%
  mutate(DOY=yday(Date))%>%
  select(-General_remarks, -Date)

###Silene
sil_dat=read.csv("G:\\My Drive\\SLU\\project\\ZackPhen\\data\\raw\\Silene phenology_10.17897_6GVG-QH42_data.txt",
                 sep="\t", header=T)

sil_dat$Date=as.POSIXct(sil_dat$Date, tz="GMT", format = "%Y-%m-%d")
#move before 1997 because totals were across all sections
sil_dat_raw=sil_dat%>%
  mutate(year=year(Date), month=month(Date), dia=day(Date), species="Silene")%>%
  mutate(DOY=yday(Date))%>%
  select( -Date)

#combine raw datasets####
plant_dat=vctrs::vec_rbind(sil_dat_raw, sax_dat_raw, sal_dat_raw, pap_dat_raw, dry_dat_raw, cas_dat_raw)%>%
  select(Plot, Section, Buds, Flowers, Senescent, year, month, dia, DOY, species)%>%
  rename(group=species)

#(C): Replace -9999 or -999 with NA:
plant_dat[plant_dat==-9999] <- NA
plant_dat[plant_dat==-999] <- NA

#calculate phenological metrics####
FlowersEx <- plyr::ddply(plant_dat, c("Plot","group","year","DOY"), summarise,
                         Flower_heads = sum(Buds,0.00001), Flowers= sum(Flowers,0.000001))
#ratio of open flowers to buds and open flowers
FlowersEx<-mutate(FlowersEx, PercentFlower = (Flowers/Flower_heads*100))
FlowersEx$UniqueID <- as.factor(interaction(FlowersEx$group,
                                            FlowersEx$Plot,
                                            FlowersEx$year,
                                            sep="_"))

##50% bloom####

Samlet50Flow = data.frame()
for (y in unique(FlowersEx$UniqueID)){
  No3<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower<=50)
  No4<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower>50)
  if (nrow(No3)==0 || nrow(No4)==0) {
    NA_dataframeF<-data.frame(a = NA)
    NA_plotF<-print(paste(y))
    NA_fileF<-cbind(NA_plotF,NA_dataframeF)
    colnames(NA_fileF)<-c("Plot","DOY_50%_flowering")
    #    filenameF <- paste(y, ".csv", sep="_flowering")
    #    write.csv(NA_fileF, file=filenameF, row.names=T)
  } else {
    LastDateF<-head(No4, n=1) #overvej noget med order i first and last date
    FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1) #FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1)
    SamletF<-rbind(FirstDateF,LastDateF)
    x1F = SamletF$DOY
    y1F = SamletF$PercentFlower
    #plot(x1F,y1F)
    model_F<- lm(y1F ~ x1F)
    #abline(model_F)
    xpred_F<-t(sapply(50,function(y1F) chemCal::inverse.predict(model_F,y1F)[1]))
    #points(xpred_F[,1],pch=16,cex=2,50,col="blue")
    #xpred_F
    ResultsF <- rbind(xpred_F)
    NamesF<-(print(paste(y)))
    CombinedF <- cbind(NamesF,ResultsF)
    colnames(CombinedF) <- c("Plot","DOY_50_flowering")
    df50<-as.data.frame(CombinedF)
    Samlet50Flow <- rbind(Samlet50Flow,df50)
  }
}

Flowering_data_50<-Samlet50Flow %>% separate(Plot, c("Species","Plot", "Year"), "_")
Flowering_data_50$DOY_50_flowering <- unlist(Flowering_data_50$DOY_50_flowering)

##10% bloom#####

Samlet10Flow = data.frame()
for (y in unique(FlowersEx$UniqueID)){
  No3<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower<=10)
  No4<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower>10)
  if (nrow(No3)==0 || nrow(No4)==0) {
    NA_dataframeF<-data.frame(a = NA)
    NA_plotF<-print(paste(y))
    NA_fileF<-cbind(NA_plotF,NA_dataframeF)
    colnames(NA_fileF)<-c("Plot","DOY_10%_flowering")
    #    filenameF <- paste(y, ".csv", sep="_flowering")
    #    write.csv(NA_fileF, file=filenameF, row.names=T)
  } else {
    LastDateF<-head(No4, n=1) #overvej noget med order i first and last date
    FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1) #FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1)
    SamletF<-rbind(FirstDateF,LastDateF)
    x1F = SamletF$DOY
    y1F = SamletF$PercentFlower
    #plot(x1F,y1F)
    model_F<- lm(y1F ~ x1F)
    #abline(model_F)
    xpred_F<-t(sapply(10,function(y1F) chemCal::inverse.predict(model_F,y1F)[1]))
    #points(xpred_F[,1],pch=16,cex=2,10,col="blue")
    #xpred_F
    ResultsF <- rbind(xpred_F)
    NamesF<-(print(paste(y)))
    CombinedF <- cbind(NamesF,ResultsF)
    colnames(CombinedF) <- c("Plot","DOY_10_flowering")
    df10<-as.data.frame(CombinedF)
    Samlet10Flow <- rbind(Samlet10Flow,df10)
  }
}

Flowering_data_10<-Samlet10Flow %>% separate(Plot, c("Species","Plot", "Year"), "_")
Flowering_data_10$DOY_10_flowering <- unlist(Flowering_data_10$DOY_10_flowering)

##90% bloom####
Samlet90Flow = data.frame()
for (y in unique(FlowersEx$UniqueID)){
  No3<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower<=90)
  No4<-subset(FlowersEx, FlowersEx$UniqueID==(print(paste(y))) & FlowersEx$PercentFlower>90)
  if (nrow(No3)==0 || nrow(No4)==0) {
    NA_dataframeF<-data.frame(a = NA)
    NA_plotF<-print(paste(y))
    NA_fileF<-cbind(NA_plotF,NA_dataframeF)
    colnames(NA_fileF)<-c("Plot","DOY_90%_flowering")
    #    filenameF <- paste(y, ".csv", sep="_flowering")
    #    write.csv(NA_fileF, file=filenameF, row.names=T)
  } else {
    LastDateF<-head(No4, n=1) #overvej noget med order i first and last date
    FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1) #FirstDateF<-tail(subset(No3, No3$DOY<max(No4$DOY)), n=1)
    SamletF<-rbind(FirstDateF,LastDateF)
    x1F = SamletF$DOY
    y1F = SamletF$PercentFlower
    #plot(x1F,y1F)
    model_F<- lm(y1F ~ x1F)
    #abline(model_F)
    xpred_F<-t(sapply(90,function(y1F) chemCal::inverse.predict(model_F,y1F)[1]))
    #points(xpred_F[,1],pch=16,cex=2,10,col="blue")
    #xpred_F
    ResultsF <- rbind(xpred_F)
    NamesF<-(print(paste(y)))
    CombinedF <- cbind(NamesF,ResultsF)
    colnames(CombinedF) <- c("Plot","DOY_90_flowering")
    df90<-as.data.frame(CombinedF)
    Samlet90Flow <- rbind(Samlet90Flow,df90)
  }
}
Flowering_data_90<-Samlet90Flow %>% separate(Plot, c("Species","Plot", "Year"), "_")
Flowering_data_90$DOY_90_flowering <- unlist(Flowering_data_90$DOY_90_flowering)

plant_time_dat=merge(Flowering_data_10, Flowering_data_50, all.x = T, all.y = T)
plant_time_dat2=merge(plant_time_dat, Flowering_data_90, all.x = T, all.y = T)

#identify where values need to be interpolated
check_missing_data = function(df,fields) {
  missing = c()
  for (n in seq(length(df[,1]))) {
    if (any(is.na(df[n,fields])))
      missing = rbind(missing, df[n,])
  }
  return(as.numeric(row.names(missing)))
}

df_col=colnames(plant_time_dat2)

na_d10=check_missing_data(plant_time_dat2, df_col)

mis_dat=plant_time_dat2[row(plant_time_dat2) %in% na_d10,]
mis_dat=mis_dat[1:length(na_d10),]

plant_timing=rbind(plant_time_dat2, mis_dat)%>%arrange(Species, Plot, Year)

plant_timing=plant_timing %>%
  distinct(DOY_10_flowering, DOY_50_flowering, DOY_90_flowering, .keep_all = TRUE)

#interpolate values####
plant_timing$DOY_10_flowering= zoo::na.approx(plant_timing$DOY_10_flowering)
plant_timing$DOY_50_flowering= zoo::na.approx(plant_timing$DOY_50_flowering)
plant_timing$DOY_90_flowering= zoo::na.approx(plant_timing$DOY_90_flowering)

#subset data to focal plots####
target_plots=c("Cas1", "Cas2", "Cas3", "Cas4", "Cas5", "Cas6",
               "Dry1", "Dry2", "Dry3", "Dry4", "Dry5", "Dry6", "Dry7", "Dry8",
               "Pap1", "Pap2", "Pap3", "Pap4",
               "Sal1", "Sal2", "Sal3", "Sal4", "Sal5", "Sal6",
               "Sax1", "Sax2", "Sax3",
               "Sil1", "Sil2", "Sil3", "Sil4")

plant_timing=plant_timing%>%filter(Plot%in%target_plots)

#write.csv(plant_timing, "ZAC_plant_phenology_1996-2023.csv")

plant_timings=plant_timing%>%
  pivot_longer(cols=c(DOY_10_flowering:DOY_90_flowering), names_to="period")


require(WaveletComp)

pap90=plant_timings%>%filter(Species=="Papaver", period=="DOY_90_flowering")%>%
  group_by(Year)%>%summarise(ave_doy90=mean(value))

pap_com90=analyze.wavelet(pap90, "ave_doy90", make.pval = TRUE, n.sim = 10)
