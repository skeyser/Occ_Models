#This script is intended to extract the norther breeding ranges of species in the BBS
#install.packages("pacman")
pacman::p_load("tidyverse", "reshape2", "vegan", "here")

#Reads in all the dataframes into a list 
temp = list.files(path = here("Data_BBS/States"), pattern="*.csv", full.names = T)
myfiles = lapply(temp, read.csv)

#Bind all states
master.df <- do.call(rbind, myfiles) 

#All 49 states are represented in this DF
length(unique(master.df$StateNum))

#Function reads all csv files in subdirectory and binds them in the tidy way 
#Throws error messages 
# tbl <- list.files(path = here("Data_BBS/States"), pattern="*.csv", full.names = T) %>%
#   map_df(~read_csv(.), col_types = cols(.default = "c"))

#Read in file for bird common names 
Species <- read.csv(here("Data_BBS/SpeciesList.csv"))

#Shape and organize species data for merging
Species <- Species[, c(2, 3, 5, 6, 7, 8)]
Species.character <- lapply(Species[2:6], as.character)
aou <- Species[,1]
Species.character <- as.data.frame(Species.character, stringsAsFactors = F)
Species <- cbind(aou, Species.character)
colnames(Species)[colnames(Species) == "aou"] <- "AOU"

#Merge DF by AOU
master.bbs <- merge(Species, master.df, by = "AOU")

#Read in Latitudes for routes 
rts <- read.csv(here("Data_BBS/routes.csv"))

#Make unique route identifier for both DF to merge 
master.bbs$Site.id <- paste0(master.bbs$StateNum, "", master.bbs$Route)
rts$Site.id <- paste0(rts$statenum, "", rts$Route)

bbs <- merge(master.bbs, rts, by = "Site.id")

#Pull out max for each AOU and Yr
bbs.max <- aggregate(Latitude ~ AOU + Year, data = bbs, FUN = "max")

#Bin Years for Each 5 yrs for avg range limit
bbs.1980 <- bbs.max #[bbs.max$Year >= 1980, ]

bbs.1980$yr.bin <- NA

#Typical assignment of year_bins 
#Need to find a more elegant way to do this 
bbs.1980$yr.bin[bbs.1980$Year >= 1966 & bbs.1980$Year < 1971] <- 1
bbs.1980$yr.bin[bbs.1980$Year >= 1971 & bbs.1980$Year < 1976] <- 2 
bbs.1980$yr.bin[bbs.1980$Year >= 1976 & bbs.1980$Year < 1981] <- 3 
bbs.1980$yr.bin[bbs.1980$Year >= 1981 & bbs.1980$Year < 1986] <- 4 
bbs.1980$yr.bin[bbs.1980$Year >= 1986 & bbs.1980$Year < 1991] <- 5 
bbs.1980$yr.bin[bbs.1980$Year >= 1991 & bbs.1980$Year < 1996] <- 6 
bbs.1980$yr.bin[bbs.1980$Year >= 1996 & bbs.1980$Year < 2001] <- 7
bbs.1980$yr.bin[bbs.1980$Year >= 2001 & bbs.1980$Year < 2006] <- 8 
bbs.1980$yr.bin[bbs.1980$Year >= 2006 & bbs.1980$Year < 2011] <- 9
bbs.1980$yr.bin[bbs.1980$Year >= 2016 & bbs.1980$Year <= 2017] <- 10

#Extract a DF that contains the range limits for each species between 1980 - 1985
bbs.1980 <- bbs.1980 %>% group_by(AOU, yr.bin) %>% summarise(Latitude = mean(Latitude)) %>% slice(which.min(yr.bin)) #%>% arrange(min.bin) %>% slice(n = 1)
max.range <- merge(bbs.1980, Species, by = "AOU")


#Remove Unidentified species 
max.clean <- max.range[!grepl("unid", max.range$English_Common_Name),]
range.merge <- max.clean[, c("Latitude", "English_Common_Name", "Genus", "Species")]

range.merge$sci_name <- paste0(range.merge$Genus, " ", range.merge$Species)

range.merge <- range.merge %>% select(-c(Genus, Species))

#write.csv(range.merge, here::here("Data_BBS/Generated DFs/MaxRange.csv"))
