#Zack arthropod raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)
require(tidyr)
require(reshape2)
require(readr)

arth_raw_dat=read.table("G:\\My Drive\\SLU\\project\\GEM-datasets\\Arthropod emergence_10.17897_V285-Z265\\Arthropod emergence_10.17897_V285-Z265_data.txt",
                      sep="\t", header=T)

colnames(arth_raw_dat)[1] <- "Date"
colnames(arth_raw_dat)[2] <- "Time"
colnames(arth_raw_dat)[colnames(arth_raw_dat) %in% c("A", "B", "C", "D", "E", "F", "G", "H")] <- c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")
arth_raw_dat$Date <- as.Date(arth_raw_dat$Date, tz="GMT")

#(C): Replace -9999 or -999 with NA:
arth_raw_dat[arth_raw_dat==-9999] <- NA
arth_raw_dat[arth_raw_dat==-999] <- NA

#renaming taxa
# Gracilariidae to Gracillariidae
# Arigades franklinii to Agriades glandon

arth_raw_dat=arth_raw_dat%>%
  mutate(Family=if_else(Family=="Gracilariidae", "Gracillariidae", Family),
        Species= if_else(Species=="Arigades franklinii", "Arigades glandon", Species))

arth_raw_dat=arth_raw_dat%>%
  mutate(Year=year(Date), Month=month(Date), Day=day(Date))%>%
  mutate(across(18:25, as.integer))%>%
  mutate(DOY=yday(Date))

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

#We sum the catches in each trap for each taxa in each trapping-period
Counts <- aggregate(cbind(Plot_A,Plot_B,Plot_C,Plot_D,Plot_E,Plot_F,Plot_G,Plot_H)~CatchID,drop=FALSE, data=arth_raw_data,  FUN = sum, na.action = na.omit)
#And take the maximum of trap days for the same trapping periods
TrapDays <- aggregate(cbind(Days.A,Days.B,Days.C,Days.D,Days.E,Days.F,Days.G,Days.H)~CatchID, drop=FALSE, data=arth_raw_data, FUN = max,na.action = na.omit)

newdat=left_join(Counts, TrapDays, by="CatchID")

arth_raw_data<- cbind(data.frame(do.call("rbind", strsplit(as.character(newdat$CatchID), "_", fixed = TRUE))), newdat)
colnames(arth_raw_data)[1:4] <- c("Plot.ID", "Year", "DOY", "HoyeTaxon")
arth_raw_data$UnitID <- paste0(arth_raw_data$Plot.ID,"_",arth_raw_data$Year, "_", arth_raw_data$DOY)

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

tail(Mites2)

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

### Calculate CatchPerTrapDay per trap per census per year ####
arth_raw_dat_full<- plyr::ddply(arth_raw_data2, c("Year","Plot.ID","DOY","HoyeTaxon"), summarise,
                     TotCatchTrapDay = sum(CatchPerTrapDay),
                     TotCatchInd=sum(TotalCatch1),
                     TotalTrapDays=mean(TrapDays))#summarize per census round per year

arth_raw_dat_full$UniqueID <- as.factor(interaction(arth_raw_dat_full$Year, arth_raw_dat_full$Plot.ID,sep="_"))
arth_raw_dat_full$FULL_name <- as.factor(interaction(arth_raw_dat_full$UniqueID, arth_raw_dat_full$HoyeTaxon,sep="_"))


arth_raw_dat_full$CumSum2<-ave(arth_raw_dat_full$TotCatchTrapDay,arth_raw_dat_full$FULL_name,FUN=cumsum)
arth_raw_dat_full$Total<-ave(arth_raw_dat_full$TotCatchTrapDay,arth_raw_dat_full$FULL_name,FUN=sum)
arth_raw_dat_full<-mutate(arth_raw_dat_full,CumFreq = 100*(CumSum2/Total))

arth_raw_dat_full <- arth_raw_dat_full[!is.na(arth_raw_dat_full$CumFreq),]

ARTPhenologyPlot<-ggplot(arth_raw_dat_full[arth_raw_dat_full$Year==2004 & arth_raw_dat_full$Plot.ID=="Art2",], aes(DOY, CumFreq), color="red") +
  #geom_abline(slope=0, intercept=50, linetype=2, color="gray60")+
  geom_line(size=1.3)+
  #geom_line(aes(DOY, PercentBuds), color="blue",size=1.3)+
  #geom_line(aes(DOY, PercentSenescent), color="brown",size=1.3)+
  #geom_smooth(method="lm", se=FALSE)+
  theme_bw()+
  #geom_abline(slope=1, intercept=0, linetype=1, color="gray60")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

ARTPhenologyPlot

#10% EMERGENCE
arth_dat_prop10=arth_raw_dat_full%>%
  group_by(Year, FULL_name)%>%
  mutate(cum_tot=cumsum(TotCatchTrapDay),
         cum_pt=cum_tot/Total)

#choose last date when approximately 10% emergence was reached
arth_grp1=arth_dat_prop10%>%
  group_by(Year, Plot.ID, HoyeTaxon, FULL_name)%>%
  mutate(PD_10=if_else(cum_pt<=0.1,"1", "2"))%>%
  filter(PD_10==1)%>%slice_tail(n=1)%>%rename(DOY_1=DOY)%>%
  select(Year, Plot.ID, DOY_1, HoyeTaxon, FULL_name)

arth_grp1$DOY_1=as.integer(arth_grp1$DOY_1)

#choose first date when approximately 10% emergence was reached
arth_grp2=arth_dat_prop10%>%
    group_by(Year, Plot.ID, HoyeTaxon, FULL_name)%>%
  mutate(PD_10=if_else(cum_pt<=0.1,"1", "2"))%>%
    filter(PD_10==2)%>%slice_head(n=1)%>%rename(DOY_2=DOY)%>%
  select(Year, Plot.ID, DOY_2, HoyeTaxon, FULL_name)

arth_grp2$DOY_2=as.integer(arth_grp2$DOY_2)

#approximate timing, date between the first and last date when >/< 10% emergence was reached
approx_dat=left_join(arth_grp1, arth_grp2)%>%mutate(interp_DOY=median(c(DOY_1,DOY_2)))

