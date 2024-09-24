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
tmp_counts <- as.integer(arth_raw_dat[index_1996, "Plot_A"] / 4) #double precision
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

#keep the row in the dataframe if any of the traps contains counts of arthropos ANd info on trap-days
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

#Change columns with counts to integer:
arth_raw_dat[c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")] <-
  sapply(arth_raw_dat[c("Plot_A", "Plot_B", "Plot_C", "Plot_D", "Plot_E", "Plot_F", "Plot_G", "Plot_H")], as.integer)

#Add total number of individuals caught per day in all open traps:
arth_raw_dat <- arth_raw_dat %>%
  rowwise() %>%
  mutate(Total_Art = sum(Plot_A, Plot_B, Plot_C, Plot_D,
                         Plot_E, Plot_F, Plot_G, Plot_H, na.rm=TRUE))

#Specify per day the number of open traps where arthropods were counted:
arth_raw_dat <- arth_raw_dat %>%
  rowwise() %>%
  mutate(Total_Traps = sum(!is.na(Plot_A), !is.na(Plot_B), !is.na(Plot_C), !is.na(Plot_D),
                           !is.na(Plot_E), !is.na(Plot_F), !is.na(Plot_G), !is.na(Plot_H)))

#Calculate average number of individuals caught of each taxon per trap per day
arth_raw_dat$Mean_Art <- arth_raw_dat$Total_Art / arth_raw_dat$Total_Traps

# #Change NA in the columns Phylum, Order, Family, Genus, Species to 'sp.'
arth_raw_dat$Phylum[is.na(arth_raw_dat$Phylum)] <- "sp."
arth_raw_dat$Order[is.na(arth_raw_dat$Order)] <- "sp."
arth_raw_dat$Family[is.na(arth_raw_dat$Family)] <- "sp."
arth_raw_dat$Genus[is.na(arth_raw_dat$Genus)] <- "sp."
arth_raw_dat$Species[is.na(arth_raw_dat$Species)] <- "sp."

#restructure time
arth_raw_dat$Date <- as.POSIXct(arth_raw_dat$datetime, tz="GMT", format = "%Y-%m-%d")
arth_raw_dat$doy <- strftime(arth_raw_dat$datetime, format = "%j", tz="GMT")

arth_raw_dat=arth_raw_dat%>%
   mutate(year=year(Date), month=month(Date), dia=day(Date))%>%
   mutate(across(16:23, as.integer))%>%
   mutate(date=as.Date(paste(year, month, dia), "%Y %m %d", tz="GMT"))%>%
   select(-Date)

write.csv(arth_raw_dat, file="ZAC_arth_raw_int.csv")

##CLEANING DATA####
#for now, include all spp. but in the future maybe exclude
#Thysanoptera (n=494), and Copepoda (n=6), Enchytraeida (n=36), Ostracoda (n=354), Psocoptera (n=1), Siphonaptera (n=3))
#arth_raw_dat <- arth_raw_dat[arth_raw_dat$Order %in% c("Araneae", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Acari", "Collembola"),]
#arth_raw_dat <- droplevels(arth_raw_dat)

#Exclude Art1 (window trap) and Art6 (no data for 1999 - 2018)
#Note that Art6 was only open in 1996, 1997, 1998 and 2019, 2020
#Note that Art4 only runs until (and including) 2016
#Note that Art7 only runs from (and including) 1999
arth_raw_dat <- arth_raw_dat[!arth_raw_dat$Plot.ID %in% c("Art1", "Art6"), ]
arth_raw_dat <- droplevels(arth_raw_dat)

#Sampling started c.a. 2 weeks later in 2020 due to covid restrictions, and 8 days later in 2018 due to excess snowcover. We
#therefore exclude both years from the analysis.
arth_raw_dat= arth_raw_dat%>%filter(!year%in%c(2018, 2020))

