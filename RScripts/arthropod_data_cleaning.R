#Zack arthropod raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)
require(tidyr)
require(reshape2)
require(readr)

file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data\\raw"

arth_name=paste(file_path, '\\Arthropod emergence_10.17897_V285-Z265_data','.txt', sep = '')
arth_raw_dat=read.csv(arth_name, header=T, sep='\t',  stringsAsFactors = F)

#date-time manipulation

colnames(arth_raw_dat)[1] <- "Date"
colnames(arth_raw_dat)[2] <- "Time"
arth_raw_dat$Date <- as.Date(arth_raw_dat$Date, "%Y-%m-%d")
arth_raw_dat=arth_raw_dat%>%
  mutate(Year=year(Date), Month=month(Date), Day=day(Date))%>%
  mutate(across(18:25, as.integer))%>%
  mutate(DOY=yday(Date))

colnames(arth_raw_dat)[colnames(arth_raw_dat) %in% c("A", "B", "C", "D", "E", "F", "G", "H")] <- c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")

#Replace -9999 or -999 with NA:
arth_raw_dat[arth_raw_dat==-9999] <- NA
arth_raw_dat[arth_raw_dat==-999] <- NA

#renaming taxa
# Gracilariidae to Gracillariidae
# Arigades franklinii to Agriades glandon

arth_raw_dat=arth_raw_dat%>%
  mutate(Family=if_else(Family=="Gracilariidae", "Gracillariidae", Family),
        Species= if_else(Species=="Arigades franklinii", "Arigades glandon", Species))

#creating a sampling unit ID
arth_raw_dat$UnitID <- paste0(arth_raw_dat$Plot.ID,"_",arth_raw_dat$Year, "_", arth_raw_dat$DOY)

#Making a new variable to sort taxa of various level as they were in
#Schmidt et al. 2016 ecography paper. This will combine some families, which
#were not separated during the first few years of monitoring into one throughout the time series

arth_raw_dat$HoyeTaxon <- "Other"
arth_raw_dat[which(arth_raw_dat$Family=="Chironomidae"),"HoyeTaxon"] <- "Chironomidae"
arth_raw_dat[which(arth_raw_dat$Family=="Ceratopogonidae"),"HoyeTaxon"] <- "Chironomidae"
arth_raw_dat[which(arth_raw_dat$Family=="Coccoidea"),"HoyeTaxon"] <- "Coccoidea"
arth_raw_dat[which(arth_raw_dat$Order=="Collembola"),"HoyeTaxon"] <- "Collembola"
arth_raw_dat[which(arth_raw_dat$Order=="Acari"),"HoyeTaxon"] <- "Acari"
arth_raw_dat[which(arth_raw_dat$Family=="Culicidae"),"HoyeTaxon"] <- "Culicidae"
arth_raw_dat[which(arth_raw_dat$Family=="Ichneumonidae"),"HoyeTaxon"] <- "Ichneumonidae"
arth_raw_dat[which(arth_raw_dat$Family=="Linyphiidae"),"HoyeTaxon"] <- "Linyphiidae"
arth_raw_dat[which(arth_raw_dat$Family=="Lycosidae"),"HoyeTaxon"] <- "Lycosidae"
arth_raw_dat[which(arth_raw_dat$Family=="Muscidae"),"HoyeTaxon"] <- "Muscidae"
arth_raw_dat[which(arth_raw_dat$Family=="Anthomyiidae"),"HoyeTaxon"] <- "Muscidae"
arth_raw_dat[which(arth_raw_dat$Family=="Anthomyzidae"),"HoyeTaxon"] <- "Muscidae"
arth_raw_dat[which(arth_raw_dat$Family=="Nymphalidae"),"HoyeTaxon"] <- "Nymphalidae"
arth_raw_dat[which(arth_raw_dat$Family=="Phoridae"),"HoyeTaxon"] <- "Phoridae"
arth_raw_dat[which(arth_raw_dat$Family=="Sciaridae"),"HoyeTaxon"] <- "Sciaridae"
arth_raw_dat[which(arth_raw_dat$Family=="Mycetophilidae"),"HoyeTaxon"] <- "Sciaridae"

#+bumblebees for controlling mites
arth_raw_dat[which(arth_raw_dat$Genus=="Bombus"),"HoyeTaxon"] <- "Bombus"

# mark lycosid egg sacks for removal
arth_raw_dat[which(arth_raw_dat$Species=="Lycosidae egg sac"),"HoyeTaxon"] <- "Other"

#subset data
arth_raw_data=arth_raw_dat%>%filter(Month%in%c(6:8), !HoyeTaxon=="Other")

#Making an even finer ID to saparate all taxa
arth_raw_data$CatchID <- paste0(arth_raw_data$Plot.ID,"_",arth_raw_data$Year,"_",arth_raw_data$DOY,"_",arth_raw_data$HoyeTaxon)
#NAs to zeros as some rows have NA for catches and trap days and zero trapdays is the same, but doesn't phase the code
arth_raw_data[is.na(arth_raw_data)] <-0
arth_raw_data=arth_raw_data%>%mutate_at(c(23:30), ~replace_na(.,0),
                                        c(10:17), ~replace_na(.,0))

#We sum the catches in each trap for each taxa in each trapping-period
Counts <- aggregate(cbind(Plot_A,Plot_B,Plot_C,Plot_D,Plot_E,Plot_F,Plot_G,Plot_H)~CatchID,drop=FALSE, data=arth_raw_data,  FUN = sum, na.action = na.omit)
#And take the maximum of trap days for the same trapping periods
TrapDays <- aggregate(cbind(Days.A,Days.B,Days.C,Days.D,Days.E,Days.F,Days.G,Days.H)~CatchID, drop=FALSE, data=arth_raw_data, FUN = max,na.action = na.omit)

newdat=left_join(Counts, TrapDays, by="CatchID")

arth_raw_data<- cbind(data.frame(do.call("rbind", strsplit(as.character(newdat$CatchID), "_", fixed = TRUE))), newdat)
colnames(arth_raw_data)[1:4] <- c("Plot.ID", "Year", "DOY", "HoyeTaxon")
arth_raw_data$UnitID <- paste0(arth_raw_data$Plot.ID,"_",arth_raw_data$Year, "_", arth_raw_data$DOY)

###    Lycosids and Mites ----
# Correcting for clusters of mites carried by bubmlebees and
# baby lycosids carried by their mothers

arth_raw_data<-mutate(arth_raw_data,TrapDays = rowSums(across(Days.A:Days.H), na.rm = T))
arth_raw_data<-mutate(arth_raw_data,TotalCatch1 = rowSums(across(Plot_A:Plot_H), na.rm = T))
arth_raw_data<-mutate(arth_raw_data,CatchPerTrapDay = (TotalCatch1/TrapDays))
arth_raw_data$TotalCatch1[is.na(arth_raw_data$TotalCatch1)] <- 0
arth_raw_data[is.na(arth_raw_data$CatchPerTrapDay),"CatchPerTrapDay"] <- 0

#### lycosids ####
#Correcting for clusters of juveline lycosids. We disregard any counts of lycosids over 3x standard deviation in trap catches

Lycosids <- arth_raw_data[which(arth_raw_data$HoyeTaxon=="Lycosidae"),]
Lh <- Lycosids[Lycosids$Year>1996 & Lycosids$Year<2006 & Lycosids$Plot.ID!="Art1",]
Lhm <- as.matrix(Lh[,c("Plot_A","Plot_B","Plot_C","Plot_D","Plot_E","Plot_F","Plot_G","Plot_H")])

THPooled <- 100
TH <- 3*sd(Lhm)+mean(Lhm)

Lycosids2 <- Lycosids[0,]

Samples <- unique(Lycosids$UnitID)
PitCatches <-c("Plot_A","Plot_B","Plot_C","Plot_D", "Plot_E", "Plot_F","Plot_G","Plot_H")
Pitfalls <- c("Days.A", "Days.B", "Days.C", "Days.D" ,"Days.E", "Days.F", "Days.G", "Days.H")

for(i in Samples){
  al  <- Lycosids[Lycosids$UnitID==i,]
  if(nrow(al)>0){
    #  if(al$Year[1]==1996){  #different rule for the pooled 1996
    #    for(j in 1:8){
    #      if(al[,PitCatches[j]]>THPooled){
    #        al[,PitCatches[j]]<-0
    #        al[,Pitfalls[j]]<-0
    #      }
    #    }
    #  }
    if(al$Year[1]>1996){ #in the first year catches from pitfalls were combined but number of juvenile spiders was low so we skip the year here
      for(j in 1:8){
        if(al[,PitCatches[j]]>TH){
          al[,PitCatches[j]]<-0
          al[,Pitfalls[j]]<-0
        }
      }
    }
  }
  Lycosids2 <- rbind(Lycosids2,al)
}

Lycosids2<-mutate(Lycosids2,TrapDays = rowSums(across(Days.A:Days.H), na.rm = T))
Lycosids2<-mutate(Lycosids2,TotalCatch1 = rowSums(across(Plot_A:Plot_H), na.rm = T))
Lycosids2<-mutate(Lycosids2,CatchPerTrapDay = (TotalCatch1/TrapDays))
Lycosids2$TotalCatch1[is.na(Lycosids2$TotalCatch1)] <- 0
Lycosids2[is.na(Lycosids2$CatchPerTrapDay),"CatchPerTrapDay"] <- 0

#### Mites ####
#For mites, we set those pitfall catches of mites which also contain a bumblebee to zero count and zero trap-days
# essentially eliminating them from affecting the indiciduals per trap day accumulation from the mean of other traps

Mites <- arth_raw_data[which(arth_raw_data$HoyeTaxon=="Acari" | arth_raw_data$HoyeTaxon=="Bombus"),]
Bombi <- arth_raw_data[which(arth_raw_data$HoyeTaxon=="Bombus" & arth_raw_data$TotalCatch1>0) ,]
Mites2 <- Mites[0,]

Mites[is.na(Mites)] <- 0

Samples <- unique(Mites$UnitID)
PitCatches <-c("Plot_A","Plot_B","Plot_C","Plot_D", "Plot_E", "Plot_F","Plot_G","Plot_H")
Pitfalls <- c("Days.A", "Days.B", "Days.C", "Days.D" ,"Days.E", "Days.F", "Days.G", "Days.H")

for(i in Samples){
  ab <- Mites[Mites$UnitID==i & Mites$HoyeTaxon=="Bombus" & Mites$TotalCatch1>0,]
  am  <- Mites[Mites$UnitID==i & Mites$HoyeTaxon=="Acari",]
  if(nrow(am)>0 & nrow(ab)>0){
    for(j in 1:8){
      if(ab[,PitCatches[j]]>0){
        am[,PitCatches[j]]<-0
        am[,Pitfalls[j]]<-0
      }
    }
  }
  Mites2 <- rbind(Mites2,am)
}

Mites2<-mutate(Mites2,TrapDays = rowSums(across(Days.A:Days.H), na.rm = T))
Mites2<-mutate(Mites2,TotalCatch1 = rowSums(across(Plot_A:Plot_H), na.rm = T))
Mites2<-mutate(Mites2,CatchPerTrapDay = (TotalCatch1/TrapDays))
Mites2$TotalCatch1[is.na(Mites2$TotalCatch1)] <- 0
Mites2[is.na(Mites2$CatchPerTrapDay),"CatchPerTrapDay"] <- 0

## Compiling corrected arthropod data ####

#Removing the uncorrected mites and lycosids from the main data
#and putting the corrected rows back in
arth_raw_data2<-arth_raw_data[which(arth_raw_data$HoyeTaxon!="Acari"& arth_raw_data$HoyeTaxon!="Lycosidae"),]

arth_raw_data2 <- rbind(arth_raw_data2,Mites2)
arth_raw_data2 <- rbind(arth_raw_data2,Lycosids2)

arth_raw_data <- arth_raw_data2

#write.csv(arth_raw_data, "arth_dat_clean.csv")

#arthropod analysis

#load packages
require(dplyr)
require(tidyr)
require(lme4)
require(ggplot2)
require(ggpubr)

#load cleana arthropod data
file_path="C:\\pdumandanSLU\\PatD-SLU\\SLU\\phenology-project\\ZackPhen\\data"
art_name=paste(file_path, '\\arth_dat_clean','.csv', sep = '')

arth_dat=read.csv(art_name, header=T, sep=',',  stringsAsFactors = F)
arth_dat=arth_dat%>%select(-c(6:23))

#inspect data crudely
#plots per taxon

arth_plots=arth_dat%>%group_by(HoyeTaxon, Year)%>%
  summarise(plot_ids=list(unique(Plot.ID)), n_plots=length(unique(Plot.ID)))

#remove Art1, Art4 and Art6

#Art1: window trap
#Art4 and 6:discontinuous opening
#Arth 7 missing from 1996-1998?

arth_data <- arth_dat %>%
  filter(!Plot.ID %in% c("Art1", "Art4", "Art6"),Year >= 1996)

#check plot and species totals, remove those with <50, and those with <5 years of data

arth_tot=arth_data%>%
  group_by(Plot.ID, Year, HoyeTaxon)%>%
  summarise(tot=sum(TotalCatch1), .groups = "drop")%>%
  filter(tot>=50)%>%
  group_by(Plot.ID, HoyeTaxon)%>%
  mutate(nyr = n_distinct(Year), .groups = "drop")%>%
  filter(nyr>=5)

arth_data50 <- arth_data %>%
  semi_join(arth_tot, by = c("Plot.ID", "Year", "HoyeTaxon"))

#check for unimodality
#include only plots and years that pass Hartigan's dip test for unimodality
library(diptest)

unimod_results <- arth_data50 %>%
  group_by(HoyeTaxon, Plot.ID, Year) %>%
  summarise(
    D_statistic = tryCatch({
      if (length(TotalCatch1) > 5) dip.test(TotalCatch1)$statistic else NA
    }, error = function(e) NA),
    p_value = tryCatch({
      if (length(TotalCatch1) > 5) dip.test(TotalCatch1)$p.value else NA
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    unimodal = ifelse(!is.na(p_value), p_value > 0.05, NA),
    Include = ifelse(unimodal == TRUE, 1, 0))

arth_df=arth_data50%>%
  semi_join(unimod_results, by=c("HoyeTaxon", "Year", "Plot.ID"))%>%
  group_by(HoyeTaxon)%>%
  mutate(plot_id = dense_rank(Plot.ID)) %>%
  ungroup()

#rescale data
arth_df$yearc=as.integer(factor(arth_df$Year))
#arth_data$years = (arth_data$Year - mean(arth_data$Year))/sd(arth_data$Year)
arth_df$DOYs = scale(arth_df$DOY, center = TRUE, scale = TRUE)[,1]
arth_df$DOYsqs = arth_df$DOYs^2
arth_df$TrapDays[arth_df$TrapDays <= 0] <- 0.001

#write.csv(arth_df, "arth_datA.csv")
# ### Calculate CatchPerTrapDay per trap per census per year ####
# arth_raw_dat_full<- plyr::ddply(arth_raw_data, c("Year","Plot.ID","DOY","HoyeTaxon"), summarise,
#                      TotCatchTrapDay = sum(CatchPerTrapDay),
#                      TotCatchInd=sum(TotalCatch1),
#                      TotalTrapDays=mean(TrapDays))#summarize per census round per year
#
# arth_raw_dat_full$UniqueID <- as.factor(interaction(arth_raw_dat_full$Year, arth_raw_dat_full$Plot.ID,sep="_"))
# arth_raw_dat_full$FULL_name <- as.factor(interaction(arth_raw_dat_full$UniqueID, arth_raw_dat_full$HoyeTaxon,sep="_"))
#
# arth_raw_dat_full$CumSum2<-ave(arth_raw_dat_full$TotCatchTrapDay,arth_raw_dat_full$FULL_name,FUN=cumsum)
# arth_raw_dat_full$Total<-ave(arth_raw_dat_full$TotCatchTrapDay,arth_raw_dat_full$FULL_name,FUN=sum)
# arth_raw_dat_full<-mutate(arth_raw_dat_full,CumFreq = 100*(CumSum2/Total))
#
# arth_raw_dat_full <- arth_raw_dat_full[!is.na(arth_raw_dat_full$CumFreq),]
