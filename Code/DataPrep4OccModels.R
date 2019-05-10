########################################################
###########Data Cleaning for Occupancy Models###########
###############Created by Spencer R Keyser##############
########################################################

##Package Loading## 

library(pacman)
pacman::p_load(here, tidyverse, reshape2, ggplot2)

##End Package Loading##

##Load in Data for Species 

#Species Records for Gulf of Mexico
bbs_tx <- read.csv(here("Data_BBS/States_GoM/States_Complete/Texas_Complete.csv"))
bbs_al <- read.csv(here("Data_BBS/States_GoM/States_Complete/Alabama_Complete.csv"))
bbs_ms <- read.csv(here("Data_BBS/States_GoM/States_Complete/Mississ_Complete.csv"))
bbs_fl <- read.csv(here("Data_BBS/States_GoM/States_Complete/Florida_Complete.csv"))
bbs_la <- read.csv(here("Data_BBS/States_GoM/States_Complete/Louisia_Complete.csv"))

#All states put together
bbs <- rbind(bbs_tx, bbs_al, bbs_fl, bbs_ms, bbs_la)

#Species Identifiers DF
species <- read.csv(here("Data_BBS/States_GoM/SpeciesList.csv"), stringsAsFactors = F)

#Weather Data will include if the route meets criteria
weather <- read.csv(here("Data_BBS/weather.csv"), stringsAsFactors = F)

#Link Species data with Route Data
bbs <- merge(bbs, species, by = "AOU")

#Link BBS data with weather data
bbs <- merge(bbs, weather, by = "RouteDataID")
bbs <- bbs[, -22:-26]

#Remove all species records that are unidentified or hybrid
bbs_total <- bbs[!grepl("unid.", bbs$English_Common_Name),]
bbs_total <- bbs_total[!grepl("hybrid", bbs_total$English_Common_Name),]

#Change the column names from the "name.x" to "name"
names(bbs_total) <- gsub(x = names(bbs_total), pattern =  ".x", replacement = "")

#Load in the Routes file for selected routes
routes.99 <- read.csv(here("Data_BBS/States_GoM/bbsrts_1999_GoM_xy_mid.csv")) 
routes.66 <- read.csv(here("Data_BBS/States_GoM/bbsrts_1966_GoM_xy_mid.csv"))
#Load in data for all of the route info
routes <- read.csv(here("Data_BBS/routes.csv"), stringsAsFactors = F)

##Using the PRISM data where the envi data is from to get the used routes
routes.66$RTENAME <- as.character(routes.66$RTENAME)
routes.66 <- routes.66[routes.66$RTENAME == "SEMINOLE HLS" | routes.66$RTENAME == "ROSELAND" | routes.66$RTENAME == "BENNDALE" | routes.66$RTENAME == "MOON LAKE" | routes.66$RTENAME == "KINGSVILLE" | routes.66$RTENAME == "RAYMONDVILLE", ]


#Sum length for routes that are divided
#Decide what is the important length
rt_simple <- routes.99[, c("rte_length", "rteno", "RTENAME", "POINT_X", "POINT_Y")]
colnames(rt_simple) <- tolower(names(rt_simple))
colnames(rt_simple)[colnames(rt_simple) == "point_x"] <- "longitude"
colnames(rt_simple)[colnames(rt_simple) == "point_y"] <- "latitude"

rt66_simple <- routes.66[, c("RTELENG", "RTENO", "RTENAME", "POINT_X", "POINT_Y")]
colnames(rt66_simple) <- tolower(names(rt66_simple))
colnames(rt66_simple)[colnames(rt66_simple) == "rteleng"] <- "rte_length"
colnames(rt66_simple)[colnames(rt66_simple) == "point_x"] <- "longitude"
colnames(rt66_simple)[colnames(rt66_simple) == "point_y"] <- "latitude"


##Combine the two time periods using rbind
rts <- rbind(rt_simple, rt66_simple)
rts <- as.tibble(rts)

##Remove duplicate Routes 
rts_final <- rts[!duplicated(rts),]

#Merge Route Information for all Routes in our study 
rts_gom <- merge(rts_final, routes, by.x = "rtename", by.y = "RouteName")
rts_gom <- rts_gom[, -4:-5]
rts_gom <- rts_gom[!duplicated(rts_gom),]
rts_gom <- rts_gom[rts_gom$statenum == 83 | rts_gom$statenum == 25 |                   
                   rts_gom$statenum == 2 | rts_gom$statenum == 52 |
                   rts_gom$statenum == 42, ]

#Create a Route No. for BBs data for merging with routes
bbs_merge <- bbs_total
bbs_merge <- bbs_merge %>% mutate(StateNum = StateNum * 1000) %>%
             mutate(rteno = StateNum + Route)


#DF with all routes for GoM selected 
bbs_gom <- merge(bbs_merge, rts_gom, by ="rteno")

#Reduce this DF to only routes after 1980
bbs_gom <- bbs_gom[bbs_gom$Year >= 1980,]

#Make each 50 stop count P/A
bbs_gom$Count10 <- ifelse(bbs_gom$Count10 > 0, 1, 0)
bbs_gom$Count20 <- ifelse(bbs_gom$Count20 > 0, 1, 0)
bbs_gom$Count30 <- ifelse(bbs_gom$Count30 > 0, 1, 0)
bbs_gom$Count40 <- ifelse(bbs_gom$Count40 > 0, 1, 0)
bbs_gom$Count50 <- ifelse(bbs_gom$Count50 > 0, 1, 0)

#Need to now make a row for each Route + Segment Count Data 
birds <- bbs_gom[, c("rteno", "Year", "Genus", "Species", "Count10",
                     "Count20", "Count30", "Count40", "Count50")]

birds$sci_name <- paste0(birds$Genus, " ", birds$Species)

#Melt the data to get observations to link with segments
birds.melt <- melt(birds, id.vars = c("rteno", "Year", "sci_name"), measure.vars = c("Count10",
                                                                         "Count20",
                                                                         "Count30",
                                                                         "Count40",
                                                                         "Count50"))

#Change from factor to character
birds.melt[, 4] <- sapply(birds.melt[, 4], as.character)

birds.melt$unique_id <- paste0(birds.melt$rteno, "_", birds.melt$Year)

#Count Numbers are now Segments 
birds.melt[birds.melt$variable == "Count10", 4] <- 1
birds.melt[birds.melt$variable == "Count20", 4] <- 2
birds.melt[birds.melt$variable == "Count30", 4] <- 3
birds.melt[birds.melt$variable == "Count40", 4] <- 4
birds.melt[birds.melt$variable == "Count50", 4] <- 5

#Make ID for route + year + segment detection histories 
birds.melt$unique_id <- paste0(birds.melt$unique_id, "_", birds.melt$variable)
colnames(birds.melt)[colnames(birds.melt) == "value"] <- "Detected"
colnames(birds.melt)[colnames(birds.melt) == "variable"] <- "Segment"
birds.melt$rt_yr <- paste0(birds.melt$rteno, "_", birds.melt$Year)

#Take the Melt Data and Attach the other variables 
bbs_gom$rt_yr <- paste0(bbs_gom$rteno, "_", bbs_gom$Year)  

bbs_gom <- bbs_gom[, c("rteno", "rt_yr", "English_Common_Name", "StateNum", 
                       "ORDER", "Family", "Genus", "Species", "Month", "Day", 
                       "BCR", "ObsN")]
bbs_gom$sci_name <- paste0(bbs_gom$Genus, " ", bbs_gom$Species)

#Extract Info Not Related to Species for merging 
bbs_gom_rtinfo <- bbs_gom[, c("rteno", "rt_yr", "StateNum",
                              "Month", "Day", "BCR", "ObsN")]

#Remove duplicated rows before merging with species detection histories 
bbs_gom_rtinfo <- bbs_gom_rtinfo[!duplicated(bbs_gom_rtinfo),]

#Merge route info with species detection 
bbs_gom_merge <- merge(birds.melt, bbs_gom_rtinfo, by = "rt_yr")

#Need to clean up a few names of species 
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Setophaga coronata audoboni"] <- "Setophaga coronata"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Colaptes auratus cafer"] <- "Colaptes auratus"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Colaptes auratus auratus"] <- "Colaptes auratus"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Aphelocoma woodhouseii"] <- "Aphelocoma californica"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Aphelocoma wollweberi"] <- "Aphelocoma ultramarina"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Antrostomus arizonae"] <- "Antrostomus vociferus"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Anas platyrhynchos diazi"] <- "Anas platyrhynchos"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Ardea herodias occidentalis"] <- "Ardea herodias"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Hydroprogne caspia"] <- "Sterna caspia"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Antigone canadensis"] <- "Grus canadensis"
bbs_gom_merge$sci_name[bbs_gom_merge$sci_name == "Spinus psaltria"] <- "Carduelis psaltria"
bbs_gom_merge <- bbs_gom_merge[!bbs_gom_merge$sci_name == "Porphyrio porphyrio",]

#Upload the functional traits 
elton <- read.csv(file = here("Functional_Traits/Functional_Traits_ESA_Jetz_updated_csv.csv"), header = TRUE, stringsAsFactors = F) 
elton$English <- as.character(elton$English)
elton$Scientific <- as.character(elton$Scientific)

#Merge the functional data with the complete merged dataset 
bbs_gom_final <- merge(bbs_gom_merge, elton, by.x = "sci_name", by.y = "Scientific")
bbs_gom_final$Category <- as.character(bbs_gom_final$Category)

#Load in the Routes with more than 20% wetland 
wetland20 <- read.csv(here("Data_BBS/Generated DFs/RouteSegs_With_20percent.csv"), stringsAsFactors = F)

#Make all items capital letters for merging 
wetland20 <- data.frame(lapply(wetland20, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

#Make a "site.link" DF to connect name with Rteno for getting routes with only 20% into bbs_gom_final
site.link <- rts_final[, c("rteno", "rtename")]
site.link$rtename <- as.character(site.link$rtename)
site.link <- site.link[!duplicated(site.link),]

#Fix site.link names
site.link$rtename <- gsub(" ", "_", site.link$rtename)
site.link[site.link$rtename == "BETAIR", 2] <- "BELAIR"
site.link[site.link$rtename == "L._ATASCOSA_NWR", 2] <- "L__ATASCOSA_NWR"
site.link[site.link$rtename == "SUGARLOAF_KEY_2", 2] <- "SUGARLOAF_KEY"
site.link[site.link$rtename == "EGLIN_A.F.B.", 2] <- "EGLIN_A_F_B_"
site.link[site.link$rtename == "SEMINOLE_HLS", 2] <- "SEMINOLE_HILLS"
site.link[site.link$rtename == "ALABAMA_PORT", 2] <- "DAUPHIN_IS_2"


#Merge the site.link df with the Wetland 20
wetland.link <- merge(wetland20, site.link, by.x = "Route", by.y = "rtename")

#Remove the Name and Link With the RTENO
link.final20 <- wetland.link[, c("rteno", "year", "Segment")]

#Create Unique_id for filtering the gom_final df with this 
link.final20$unique_id <- paste0(link.final20$rteno, "_", link.final20$year, "_", link.final20$Segment)
link.final20 <- link.final20[!duplicated(link.final20$unique_id),]
rownames(link.final20) <- link.final20$unique_id

#Final DF that contains all of the Routes and Species for segments in all years that have 20% or more cover
final_gom <- bbs_gom_final[bbs_gom_final$unique_id %in% rownames(link.final20),]

#write.csv(final_gom, file = here("Data_BBS/Generated DFs/Final_GoM_DF.csv"))

final_gom <- read.csv(here("Data_BBS/Generated DFs/Final_GoM_DF.csv"))
length(unique(final_gom$Order))

       