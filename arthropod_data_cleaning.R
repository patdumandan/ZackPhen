#Zack arthropod raw data processing

require(lubridate)
require(dplyr)
require(ggplot2)

arth_raw_dat=read.table("G:\\My Drive\\SLU\\project\\GEM-datasets\\Arthropod emergence_10.17897_V285-Z265\\Arthropod emergence_10.17897_V285-Z265_data.txt",
                      sep="\t", header=T)

colnames(arth_raw_dat)[1] <- "Date"
colnames(arth_raw_dat)[2] <- "Time"
colnames(arth_raw_dat)[colnames(arth_raw_dat) %in% c("A", "B", "C", "D", "E", "F", "G", "H")] <- c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")
arth_raw_dat$Date <- as.Date(arth_raw_dat$Date)

#(B): Change date and time into a datetime column

#If time is not specified, set to noon
arth_raw_dat$Time[is.na(arth_raw_dat$Time)] <- 12
arth_raw_dat$Time[arth_raw_dat$Time < 0] <- 12

#Change time annotation:
arth_raw_dat$Time[arth_raw_dat$Time<10 & !is.na(arth_raw_dat$Time)] <- paste0(0, arth_raw_dat$Time[arth_raw_dat$Time<10 & !is.na(arth_raw_dat$Time)] )
arth_raw_dat$Time[!is.na(arth_raw_dat$Time)] <-  paste0(arth_raw_dat$Time[!is.na(arth_raw_dat$Time)], ":00:00")

#Create datetime:
arth_raw_dat$datetime <- paste0(arth_raw_dat$Date, " ", arth_raw_dat$Time)
arth_raw_dat$datetime <- as.POSIXct(arth_raw_dat$datetime, tz="GMT", format = "%Y-%m-%d %H:%M:%S")

#(C): Replace -9999 or -999 with NA:
arth_raw_dat[arth_raw_dat==-9999] <- NA
arth_raw_dat[arth_raw_dat==-999] <- NA

#1996 data adjustments
index_1996 <- grepl( "1996", arth_raw_dat$Date, fixed = TRUE)
tmp_days <- round(arth_raw_dat[index_1996, "Days.A"] / 4, 0) #make sure these are integers
tmp_counts <- as.numeric(arth_raw_dat[index_1996, "Plot_A"] / 4) #double precision
arth_raw_dat[index_1996, "Days.A"] <- tmp_days
arth_raw_dat[index_1996, "Days.B"] <- tmp_days
arth_raw_dat[index_1996, "Days.C"] <- tmp_days
arth_raw_dat[index_1996, "Days.D"] <- tmp_days
arth_raw_dat[index_1996, "Plot_A"] <- tmp_counts
arth_raw_dat[index_1996, "Plot_B"] <- tmp_counts
arth_raw_dat[index_1996, "Plot_C"] <- tmp_counts
arth_raw_dat[index_1996, "Plot_D"] <- tmp_counts
rm(index_1996, tmp_days, tmp_counts)

#2010 data adjustments
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Days.C"] <-
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Days.B"] / 2
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Days.B"] <-
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Days.B"] / 2

#Adjust caught arthropods (equally divide the number in Plot_B over the columns Plot_B and Plot_C, but only allow for integers)
art_tmp <- arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Plot_B"] / 2
set.seed(23)
random <- runif(length(art_tmp),0,1)
index <- (random > 0.5) & (art_tmp != 0)
index2 <- (random <= 0.5) & (art_tmp != 0)
art_tmp[index] <- ceiling(art_tmp[index])
art_tmp[index2] <- floor(art_tmp[index2])
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Plot_C"] <- art_tmp
arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Plot_B"] <-
  (arth_raw_dat[arth_raw_dat$Sorting_remarks %in% "Samples B & C lumped together by mistake. Trap days corrected accordingly (NMS)", "Plot_B"] - art_tmp)
rm(art_tmp, random, index, index2)

arth_raw_dat$use.A <- !is.na(arth_raw_dat$Days.A) & !is.na(arth_raw_dat$Plot_A) & arth_raw_dat$Days.A!=0
arth_raw_dat$use.B <- !is.na(arth_raw_dat$Days.B) & !is.na(arth_raw_dat$Plot_B) & arth_raw_dat$Days.B!=0
arth_raw_dat$use.C <- !is.na(arth_raw_dat$Days.C) & !is.na(arth_raw_dat$Plot_C) & arth_raw_dat$Days.C!=0
arth_raw_dat$use.D <- !is.na(arth_raw_dat$Days.D) & !is.na(arth_raw_dat$Plot_D) & arth_raw_dat$Days.D!=0
arth_raw_dat$use.E <- !is.na(arth_raw_dat$Days.E) & !is.na(arth_raw_dat$Plot_E) & arth_raw_dat$Days.E!=0
arth_raw_dat$use.F <- !is.na(arth_raw_dat$Days.F) & !is.na(arth_raw_dat$Plot_F) & arth_raw_dat$Days.F!=0
arth_raw_dat$use.G <- !is.na(arth_raw_dat$Days.G) & !is.na(arth_raw_dat$Plot_G) & arth_raw_dat$Days.G!=0
arth_raw_dat$use.H <- !is.na(arth_raw_dat$Days.H) & !is.na(arth_raw_dat$Plot_H) & arth_raw_dat$Days.H!=0

#keep the row in the dataframe if any of the traps contains ccounts of arthropos ANd info on trap-days
arth_raw_dat <- arth_raw_dat %>%
  rowwise() %>%
  mutate(keep=any(use.A, use.B, use.C, use.D, use.E, use.F, use.G, use.H))

#Exclude rows where none of the traps were open, or insects were not counted:
#tmp <- arth_raw_dat[arth_raw_dat$keep==FALSE,]
arth_raw_dat <- arth_raw_dat[arth_raw_dat$keep==TRUE,]

#Exclude unused columns:
use <- arth_raw_dat[,c("use.A", "use.B", "use.C", "use.D", "use.E", "use.F", "use.G", "use.H")]
arth_raw_dat <- arth_raw_dat[,c("datetime", "Plot.ID", "Days.A", "Days.B", "Days.C", "Days.D", "Days.E", "Days.F", "Days.G", "Days.H",
            "Phylum", "Order", "Family", "Genus", "Species", "Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E",
            "Plot_F", "Plot_G", "Plot_H")]

#Change columns with counts to numeric:
arth_raw_dat[c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")] <-
  sapply(arth_raw_dat[c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")], as.numeric)

#major data adjustments: ensure even distribution of captures in a week

#Create an empty list:
arth_raw_dat_list <- list()

#Iterate through all rows of the dataframe and create separate entries for each day that traps were open:
for(i in 1:nrow(arth_raw_dat)){

  #for debugging:
  #i=2964

  #Create a temporary empty dataframe:
  arth_raw_dat_tmp <- arth_raw_dat[0,]

  #Store the row for which all steps have to be conducted:
  tmp_row <- arth_raw_dat[i,]

  #Create an index parameter that shows which columns both have arthropod counts and trap-days registered:
  use_column <- use[i,]
  index <- which(use_column==TRUE)+2

  #Store sum of counts before processing (for debugging)
  sum_before <- sum(tmp_row[1, index+13] %>% replace(is.na(.), 0))

  #what is the maximum number of trap days among all traps of this row of the dataframe:
  max_trap_days <- max(tmp_row[,index])

  #Create extra rows containing the datetime for all trap days:
  for(c in 1:max_trap_days){
    row <- max_trap_days - c
    arth_raw_dat_tmp[c, 'datetime'] <- tmp_row$datetime - (row*60*60*24)
  }

  #update all other rows that should contain the same information as our starting row:
  arth_raw_dat_tmp[1:max_trap_days, c("Plot.ID", "Phylum", "Order", "Family", "Genus", "Species")] <-
    tmp_row[rep(1, each = max_trap_days), c("Plot.ID", "Phylum", "Order", "Family", "Genus", "Species")]

  #Add the average number of arthropods caught for each separate trapping day (use average values):

  #Select columns that have arthropod counts:
  for(r in index){

    #Determine number of trap-days
    days <- tmp_row[,r]

    #If number of trap-days is 1, then do not take an average
    if(days==1){
      arth_raw_dat_tmp[max_trap_days, r+13] <- tmp_row[,r+13]
      arth_raw_dat_tmp[max_trap_days, r] <- tmp_row[,r]
    }

    #If number of trap days is larger than 1, then take averages for each trap-day:
    if(days>1){
      arth_raw_dat_tmp[as.numeric(max_trap_days-days+1):max_trap_days, r+13] <- as.numeric(rep(tmp_row[, r+13], days))/days
      arth_raw_dat_tmp[max_trap_days, r] <- tmp_row[,r]
    }

  }

  #Store sum of counts after processing (for debugging)
  sum_after <- sum(arth_raw_dat_tmp[1:max_trap_days, index+13] %>% replace(is.na(.), 0))

  #Check if sum_before==sum_after
  if(round(sum_before, 0)!=round(sum_after, 0)){
    print("Stop! Sums not equal")
    print(i)
    break
  }

  #Store dataframe in the pre-specified list:
  arth_raw_dat_list[[i]] <- arth_raw_dat_tmp

  #Print progress:
  print(paste0(round((100*i)/nrow(arth_raw_dat),3), " %"))

}

#Combine list into dataframe:
arth_raw_dat_2 = do.call(rbind, arth_raw_dat_list)

#remove obsolete variables:
rm(days, arth_raw_dat_tmp, tmp_row, use_column, c, i, index, max_trap_days, r, row, arth_raw_dat_list, arth_raw_dat, use)

#Rename dataframe:
arth_raw_dat <- arth_raw_dat_2 ; rm(arth_raw_dat_2, sum_before, sum_after)

#Add total number of individuals caught per day in all open traps:
arth_raw_dat <- arth_raw_dat %>%
  rowwise() %>%
  mutate(Total_Art = sum(Plot_A, Plot_B, Plot_C, Plot_D, Plot_E, Plot_F, Plot_G, Plot_H, na.rm=TRUE))

#Specify per day the number of open traps where arthropods were counted:
arth_raw_dat <- arth_raw_dat %>%
  rowwise() %>%
  mutate(Total_Traps = sum(!is.na(Plot_A), !is.na(Plot_B), !is.na(Plot_C), !is.na(Plot_D),
                           !is.na(Plot_E), !is.na(Plot_F), !is.na(Plot_G), !is.na(Plot_H)))

#Calculate average number of individuals caught of each taxon per trap per day
arth_raw_dat$Mean_Art <- arth_raw_dat$Total_Art / arth_raw_dat$Total_Traps

# #Add day of year (Julian Day, taking into account leap years):
# arth_raw_dat$doy <- strftime(arth_raw_dat$datetime, format = "%j", tz="GMT")
#
# #add year
# arth_raw_dat$year <- as.factor(strftime(arth_raw_dat$datetime, format = "%Y", tz="GMT"))


#Change NA in the columns Phylum, Order, Family, Genus, Species to 'sp.'
arth_raw_dat$Phylum[is.na(arth_raw_dat$Phylum)] <- "sp."
arth_raw_dat$Order[is.na(arth_raw_dat$Order)] <- "sp."
arth_raw_dat$Family[is.na(arth_raw_dat$Family)] <- "sp."
arth_raw_dat$Genus[is.na(arth_raw_dat$Genus)] <- "sp."
arth_raw_dat$Species[is.na(arth_raw_dat$Species)] <- "sp."

#restructure time
arth_raw_dat$Date <- as.POSIXct(arth_raw_dat$datetime, tz="GMT", format = "%Y-%m-%d")
arth_raw_dat$doy <- strftime(arth_raw_dat$datetime, format = "%j", tz="GMT")

# arth_raw_dat=arth_raw_dat%>%
#   mutate(year=year(Date), month=month(Date), dia=day(Date))%>%
#   mutate(across(16:23, as.integer))%>%
#   mutate(date=as.Date(paste(year, month, dia), "%Y %m %d"), DOY=yday(date))

#Exclude Acari (n=196033), Collembola (n=230643), Thysanoptera (n=494), (and Copepoda (n=6), Enchytraeida (n=36), Ostracoda (n=354), Psocoptera (n=1), Siphonaptera (n=3))
arth_raw_dat <- arth_raw_dat[arth_raw_dat$Order %in% c("Araneae", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera"),]
arth_raw_dat <- droplevels(arth_raw_dat)

#Exclude Art1 (window trap) and Art6 (no data for 1999 - 2018)
#Note that Art6 was only open in 1996, 1997, 1998 and 2019, 2020
#Note that Art4 only runs until (and including) 2016
#Note that Art7 only runs from (and including) 1999
arth_raw_dat <- arth_raw_dat[!arth_raw_dat$Plot.ID %in% c("Art1", "Art6"), ]
arth_raw_dat <- droplevels(arth_raw_dat)

#Sampling periods differ among years. To be able to make comparisons we restrict the dataset to a fixed sampling window (doy - doy)
arth_raw_dat= arth_raw_dat[!arth_raw_dat$doy<157, arth_raw_dat$doy>238,]
arth_raw_dat= arth_raw_dat[!arth_raw_dat$doy<157, arth_raw_dat$doy>238,]

doy_range <- cbind(with(df_Art_date, aggregate(doy, by=list(year), min)),  with(df_Art_date, aggregate(doy, by=list(year), max))$x)
colnames(doy_range) <- c("year", "doy_min", "doy_max")

#Sampling started c.a. 2 weeks later in 2020 due to covid restrictions, and 8 days later in 2018 due to excess snowcover. We
#therefore exclude both years from the analysis.
years <- sort(unique(df_Art_date$year))
years <- years[-which(years %in% c("2018", "2020"))]

#All remaining years have a sampling window that includes the period from doy 157 to doy 238
doy_min=157 #152 is first of June. We chose doy 157
doy_max=238 #243 is 31st of August. We chose doy 238
#Note: 2015 starts at doy 159, and 2018 at doy 165
#We thus miss the first few days of these two seasons, which is OK as arthropod numbers are very small early in the season

#Exclude all data outside that window
df_Art_date <- df_Art_date[df_Art_date$doy >= doy_min & df_Art_date$doy <= doy_max,]
df_Art_date <- droplevels(df_Art_date)

write.csv(arth_raw_dat, file="ZAC_arthdat_clean.csv", row.names = F)


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

#Exclude bombus if still in data and 2020? and an unidentified family
arth_raw_dat5=anti_join(arth_dat4,arth_dat5)%>%filter(!Genus=="Bombus", !year==2020,
                                                      !Family=="unidentified")

#exclude mites for 1996 because cups were not calculated separately
arth_dat6=arth_raw_dat5[which(arth_raw_dat5$Order=="Acari" &
                              arth_raw_dat5$year==1996 &
                              arth_raw_dat5$Plot.ID=="Art1"),]

arth_raw_dat6=anti_join(arth_raw_dat5,arth_dat6)

library(ggplot2)

ggplot(data=arth_raw_dat6, aes(x=DOY, y=TotalSpecimens, col=Family))+
  geom_point()+facet_wrap(~Order)+theme_classic()

