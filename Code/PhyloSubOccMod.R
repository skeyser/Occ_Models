########################################################
###########Data Cleaning for Occupancy Models###########
###############Created by Spencer R Keyser##############
########################################################

rm(list = ls()) 

##Package Loading## 

library(pacman)
pacman::p_load(here, tidyverse, reshape2, ggplot2, data.table, lubridate, stringr)

#library( dplyr )
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
#library( tidyr )
##End Package Loading##

##Load in Data for Species 

#Species Records for Gulf of Mexico
bbs_tx <- read.csv(here::here("Data_BBS/States_GoM/States_Complete/Texas_Complete.csv"))
bbs_al <- read.csv(here::here("Data_BBS/States_GoM/States_Complete/Alabama_Complete.csv"))
bbs_ms <- read.csv(here::here("Data_BBS/States_GoM/States_Complete/Mississ_Complete.csv"))
bbs_fl <- read.csv(here::here("Data_BBS/States_GoM/States_Complete/Florida_Complete.csv"))
bbs_la <- read.csv(here::here("Data_BBS/States_GoM/States_Complete/Louisia_Complete.csv"))

#All states put together
bbs <- rbind(bbs_tx, bbs_al, bbs_fl, bbs_ms, bbs_la)

#Species Identifiers DF
species <- read.csv(here::here("Data_BBS/States_GoM/SpeciesList.csv"), stringsAsFactors = F)

#Weather Data will include if the route meets criteria
weather <- read.csv(here::here("Data_BBS/weather.csv"), stringsAsFactors = F)

#Link Species data with Route Data
bbs <- merge(bbs, species, by = "AOU")

#Link BBS data with weather data
bbs <- merge(bbs, weather, by = "RouteDataID")

#Remove rows that don't meet BBS standards
bbs <- bbs[!bbs$RunType == 0, ]

bbs <- bbs[, -22:-26]

#Remove all species records that are unidentified or hybrid
bbs_total <- bbs[!grepl("unid.", bbs$English_Common_Name),]
bbs_total <- bbs_total[!grepl("hybrid", bbs_total$English_Common_Name),]

#Change the column names from the "name.x" to "name"
names(bbs_total) <- gsub(x = names(bbs_total), pattern =  ".x", replacement = "")

#Load in the Routes file for selected routes
routes.99 <- read.csv(here::here("Data_BBS/States_GoM/bbsrts_1999_GoM_xy_mid.csv")) 
routes.66 <- read.csv(here::here("Data_BBS/States_GoM/bbsrts_1966_GoM_xy_mid.csv"))
#Load in data for all of the route info
routes <- read.csv(here::here("Data_BBS/routes.csv"), stringsAsFactors = F)

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
                       "BCR", "ObsN", "StartTime", "EndTime")]
bbs_gom$sci_name <- paste0(bbs_gom$Genus, " ", bbs_gom$Species)

#Extract Info Not Related to Species for merging 
bbs_gom_rtinfo <- bbs_gom[, c("rteno", "rt_yr", "StateNum",
                              "Month", "Day", "BCR", "ObsN", "StartTime", "EndTime")]

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
elton <- read.csv(file = here::here("Functional_Traits/Functional_Traits_ESA_Jetz_updated_csv.csv"), header = TRUE, stringsAsFactors = F) 
elton$English <- as.character(elton$English)
elton$Scientific <- as.character(elton$Scientific)

#Merge the functional data with the complete merged dataset 
bbs_gom_final <- merge(bbs_gom_merge, elton, by.x = "sci_name", by.y = "Scientific")
bbs_gom_final$site <- paste0(bbs_gom_final$rteno.x, "_", bbs_gom_final$Segment)
bbs_gom_final$Category <- as.character(bbs_gom_final$Category)

#Load in the Routes with more than 20% wetland 
wetland20 <- read.csv(here::here("Data_BBS/Generated DFs/CompleteSegments_With_20percent.csv"), stringsAsFactors = F)

#Make all items capital letters for merging 
wetland20 <- data.frame(lapply(wetland20, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))

wetland20.site <- unique(wetland20$site)
wetland20.site <- as.vector(wetland20.site)

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

wetland.link$rteno.segment <- paste0(wetland.link$rteno, "_", wetland.link$Segment) 
wetland20.site <- unique(wetland.link$rteno.segment)
wetland20.site <- as.vector(wetland20.site)

#Remove the Name and Link With the RTENO
link.final20 <- wetland.link[, c("rteno", "year", "Segment")]

#Create Unique_id for filtering the gom_final df with this 
link.final20$unique_id <- paste0(link.final20$rteno, "_", link.final20$year, "_", link.final20$Segment)
link.final20 <- link.final20[!duplicated(link.final20$unique_id),]
rownames(link.final20) <- link.final20$unique_id

#Final DF that contains all of the Routes and Species for segments in all years that have 20% or more cover
#DF only has species detections for each year of LULC data
final_gom <- bbs_gom_final[bbs_gom_final$unique_id %in% rownames(link.final20),]

#Final DF with all species for all years for sites w/ >=20%
final_sp_df <- bbs_gom_final[bbs_gom_final$site %in% wetland20.site,]

#write.csv(final_gom, file = here("Data_BBS/Generated DFs/Final_GoM_DF.csv"))
#write.csv(final_sp_df, file = here("Data_BBS/Generated DFs/Final_Sp_GoM_DF.csv"))


#final_gom <- read.csv(here("Data_BBS/Generated DFs/Final_GoM_DF.csv"))
length(unique(final_sp_df$Order))
length(unique(final_sp_df$sci_name))


#####################################################################
###################Creating Groups for species#######################
#####################################################################

#Split DF up by BCR
bcr31 <- final_sp_df[final_sp_df$BCR == 31,]
bcr26 <- final_sp_df[final_sp_df$BCR == 26,]
bcr27 <- final_sp_df[final_sp_df$BCR == 27,]
bcr37 <- final_sp_df[final_sp_df$BCR == 37,]
bcr36 <- final_sp_df[final_sp_df$BCR == 36,]

#Calculate # of species in BCR, Category, and Order 
BCR_breakdown <- setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = BCR]
#write.csv(BCR_breakdown, file = here::here("Data_BBS/Generated DFs/BCR_Count.csv"))

setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = Category]
order_breakdown <- setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = Order]
order_breakdown$Order[order_breakdown$Order == "Apodiformes"] <- "Apodiformes/Caprimulgiformes"
order_breakdown$Order[order_breakdown$Order == "Caprimulgiformes"] <- "Apodiformes/Caprimulgiformes"
order_breakdown$Order[order_breakdown$Order == "Falconiformes"] <- "Falconiformes/Psittaciformes"
order_breakdown$Order[order_breakdown$Order == "Psittaciformes"] <- "Falconiformes/Psittaciformes"
order_breakdown$Order[order_breakdown$Order == "Coraciiformes"] <- "Coraciiformes/Piciformes"
order_breakdown$Order[order_breakdown$Order == "Piciformes"] <- "Coraciiformes/Piciformes"
order_breakdown$Order[order_breakdown$Order == "Ciconiiformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes" 
order_breakdown$Order[order_breakdown$Order == "Pelecaniformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes"
order_breakdown$Order[order_breakdown$Order == "Suliformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes"
order_breakdown$Order[order_breakdown$Order == "Podicipediformes"] <- "Podicipediformes/Gruiformes"
order_breakdown$Order[order_breakdown$Order == "Gruiformes"] <- "Podicipediformes/Gruiformes"

order_breakdown <- aggregate(data = order_breakdown, count ~ Order, FUN = sum)
#write.csv(order_breakdown, file = here("Data_BBS/Generated DFs/PhyloGroups.csv"))

setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = Nocturnal]

passerine <- final_sp_df[final_sp_df$Order == "Passeriformes",]
passer <- setDT(passerine)[, .(count = uniqueN(sci_name)), by = Family_Latin]

#Rename phylo in the complete DF for grouping species
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Sturnidae"] <- "Sturnidae/Turdidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Turdidae"] <- "Sturnidae/Turdidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Alaudidae"] <- "Remizidae/Paridae/Alaudidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Paridae"] <- "Remizidae/Paridae/Alaudidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Remizidae"] <- "Remizidae/Paridae/Alaudidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Corvidae"] <- "Laniidae/Corvidae/Vireonidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Laniidae"] <- "Laniidae/Corvidae/Vireonidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Vireonidae"] <- "Laniidae/Corvidae/Vireonidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Fringillidae"] <- "Passeridae/Cardinalidae/Fringillidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Cardinalidae"] <- "Passeridae/Cardinalidae/Fringillidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Passeridae"] <- "Passeridae/Cardinalidae/Fringillidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Troglodytidae"] <- "Polioptilidae/Sittidae/Troglodytidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Sittidae"] <- "Polioptilidae/Sittidae/Troglodytidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Polioptilidae"] <- "Polioptilidae/Sittidae/Troglodytidae"

final_sp_df$Order[final_sp_df$Order == "Apodiformes"] <- "Apodiformes/Caprimulgiformes"
final_sp_df$Order[final_sp_df$Order == "Caprimulgiformes"] <- "Apodiformes/Caprimulgiformes"
final_sp_df$Order[final_sp_df$Order == "Falconiformes"] <- "Falconiformes/Psittaciformes"
final_sp_df$Order[final_sp_df$Order == "Psittaciformes"] <- "Falconiformes/Psittaciformes"
final_sp_df$Order[final_sp_df$Order == "Coraciiformes"] <- "Coraciiformes/Piciformes"
final_sp_df$Order[final_sp_df$Order == "Piciformes"] <- "Coraciiformes/Piciformes"
final_sp_df$Order[final_sp_df$Order == "Ciconiiformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes" 
final_sp_df$Order[final_sp_df$Order == "Pelecaniformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes"
final_sp_df$Order[final_sp_df$Order == "Suliformes"] <- "Ciconiiformes/Pelecaniformes/Suliformes"
final_sp_df$Order[final_sp_df$Order == "Podicipediformes"] <- "Podicipediformes/Gruiformes"
final_sp_df$Order[final_sp_df$Order == "Gruiformes"] <- "Podicipediformes/Gruiformes"

####Remove BCR 19#####
final_sp_df <- final_sp_df[!(final_sp_df$BCR == 19),]


########################################################
############From this point on every DF generated#######
####should begin from the final_sp_df for consistency###
########################################################

order.count <- setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = Order]
family.count <- setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = Family_Latin]
########################################################
############Creation of phylogenetic subgroup###########

final_sp_df <- final_sp_df[final_sp_df$Family_Latin == "Passeridae/Cardinalidae/Fringillidae",]

########################################################


final_sp_df <- final_sp_df[!final_sp_df$Detected == 0,]
final_sp_df <- plyr::rename(final_sp_df, c("Order" = "Phylo.V2", "Family_Latin" = "Phylo.V1"))

#Create Phylo Class 1 (Orders + Families for Passerines)
final_sp_df$Phylo.V1 <- ifelse(final_sp_df$Phylo.V2 != "Passeriformes", final_sp_df$Phylo.V2, final_sp_df$Phylo.V1)


final_sp_df <- transform(final_sp_df, spp.id = as.numeric(interaction(sci_name, drop = T)))
final_sp_df <- transform(final_sp_df, rteno.id = as.numeric(interaction(rteno.x, drop = T)))
final_sp_df <- transform(final_sp_df, year.id = as.numeric(interaction(Year, drop = T)))
final_sp_df <- transform(final_sp_df, site.id = as.numeric(interaction(site, drop = T)))
final_sp_df <- transform(final_sp_df, Phylo.V1.code = as.numeric(interaction(Phylo.V1, drop = T)))
final_sp_df <- transform(final_sp_df, Phylo.V2.code = as.numeric(interaction(Phylo.V2, drop = T)))
final_sp_df <- transform(final_sp_df, bcr.id = as.numeric(interaction(BCR, drop = T)))
final_sp_df <- transform(final_sp_df, ObsN.id = as.numeric(interaction(ObsN, drop = T)))

final_sp_df <- final_sp_df %>% select(sci_name, English, spp.id, Detected, Phylo.V1, Phylo.V1.code,
                                      Phylo.V2, Phylo.V2.code, rteno.x, rteno.id, rt_yr, Segment, 
                                      Year, year.id, Day, Month, BCR, bcr.id, site, site.id, 
                                      ObsN, ObsN.id, BodyMass, Nocturnal, Category, Inverts,
                                      Mammal.Bird, Herps, Fish, Unknown_Verts, Fruit,
                                      Nect, Seed, Other_plant, ForStrat_watbelowsurf, 
                                      ForStrat_wataroundsurf, ForStrat_ground, 
                                      ForStrat_understory, ForStrat_midhigh, 
                                      ForStrat_canopy, ForStrat_aerial, 
                                      PelagicSpecialist, StartTime, EndTime)  
#Create"group" identifier for each family group or order

#Write csv for final_sp_df
#write.csv(final_sp_df, file = here::here("Data_BBS/Generated DFs/Final_Sp_Df_DetectionsOnly.csv"))

##########################


#Compute total # of detections per Family and Order
obs.phylo.fam <- aggregate(Detected ~ Family_Latin + Order, final_sp_df, FUN = sum)
obs.phylo.fam <- obs.phylo.fam[obs.phylo.fam$Order == "Passeriformes",]

obs.phylo.order <- aggregate(Detected ~ Order, final_sp_df, FUN = sum)
obs.phylo.order$Family_Latin <- NA

#Rearrange and bind order data and family data
obs.phylo.fam <- obs.phylo.fam %>% select(Order, Family_Latin, Detected)
obs.phylo.order <- obs.phylo.order %>% select(Order, Family_Latin, Detected)
obs.phylo <- rbind(obs.phylo.fam, obs.phylo.order)




#############################################################
##########Create augmented species x site matrix#############
#############################################################

#Take all of the sites where wetland >=20% from entire dataset
#total.sp.mat <- bbs_gom_final[bbs_gom_final$site %in% wetland20.site,]
total.sp.mat <- final_sp_df
total.sp.mat$unique_id <- paste0(total.sp.mat$site, "_", total.sp.mat$Year)
raw.sp.mat <- total.sp.mat[, c("unique_id", "sci_name", "Detected")]
raw.sp.mat <- dcast(raw.sp.mat, unique_id ~ sci_name, fun.aggregate = sum, value.var = "Detected")
#rownames(raw.sp.mat) <- raw.sp.mat$unique_id
#raw.sp.mat <- subset(raw.sp.mat, select = -unique_id)
raw.sp.mat[raw.sp.mat > 1] <- 1

#Calculate # of species in each Category and Order by BCR
setDT(bcr31)[, .(count = uniqueN(sci_name)), by = Category]
setDT(bcr31)[, .(count = uniqueN(sci_name)), by = Order]
setDT(bcr26)[, .(count = uniqueN(sci_name)), by = Category]
setDT(bcr26)[, .(count = uniqueN(sci_name)), by = Order]
setDT(bcr27)[, .(count = uniqueN(sci_name)), by = Category]
setDT(bcr27)[, .(count = uniqueN(sci_name)), by = Order]
setDT(bcr37)[, .(count = uniqueN(sci_name)), by = Category]
setDT(bcr37)[, .(count = uniqueN(sci_name)), by = Order]
setDT(bcr36)[, .(count = uniqueN(sci_name)), by = Category]
setDT(bcr36)[, .(count = uniqueN(sci_name)), by = Order]

#Extract a list of species for phylogenetic subset
#spp.occ <- as.data.frame(unique(final_sp_df$sci_name))

###################################
####Extracting for JAGS coding#####
###################################

#Pull out BCRs of species detections
bcr.detections.uncast <- final_sp_df %>% select(sci_name, spp.id, bcr.id, Detected)

#BCR by Species Matrix for later
bcr.detections <- dcast(bcr.detections.uncast, spp.id ~ bcr.id, fun.aggregate = sum, value.var = "Detected")
bcr.detections <- data.frame(bcr.detections, row.names = 1)
bcr.detections[bcr.detections > 1] <- 1
colnames(bcr.detections)[colnames(bcr.detections) == "X1"] <- "BCR26"
colnames(bcr.detections)[colnames(bcr.detections) == "X2"] <- "BCR27"
colnames(bcr.detections)[colnames(bcr.detections) == "X3"] <- "BCR31"
colnames(bcr.detections)[colnames(bcr.detections) == "X4"] <- "BCR36"
colnames(bcr.detections)[colnames(bcr.detections) == "X5"] <- "BCR37"
bcr.detections$spp.id <- rownames(bcr.detections) 
rownames(bcr.detections) <- NULL
#bcr.detections <- bcr.detections %>% mutate(Species.code = 1:n())
bcr.detections <- bcr.detections %>% select(spp.id, BCR26, BCR27, BCR31, BCR36, BCR37)
bcr.detections$spp.id <- as.numeric(bcr.detections$spp.id)
bcr.jags <- bcr.detections[-1:-1]
bcr.jags <- as.matrix(bcr.jags)

#Write csv
#write.table(bcr.jags, file = here::here("Data_BBS/Generated DFs/bcr.detection.raw.csv"), sep = ",", row.names = F, col.names = F)

#Write BCR 
#write.csv(bcr.detections, file = here::here("Data_BBS/Generated DFs/BCR.Detections.csv"), row.names = F)

#Put BCR detections with species for the species DF
spp.occ <- final_sp_df %>% select(sci_name, spp.id, Phylo.V1, Phylo.V1.code, 
                                  Phylo.V2, Phylo.V2.code, BodyMass)
spp.occ <- spp.occ[!duplicated(spp.occ),]
spp.occ <- merge(spp.occ, bcr.detections, by = "spp.id")
spp.occ <- spp.occ %>% select(sci_name, spp.id, Phylo.V1, 
                              Phylo.V1.code, Phylo.V2, 
                              Phylo.V2.code, BodyMass)
#view
head( spp.occ)

###Write Species DF csv
#write.csv(spp.occ, file = here::here("Data_BBS/Generated DFs/Spp.Occ.csv"), row.names = F)

###################################
###################################

#BCR matrix 
bcr.occ <- as.data.frame(unique(final_sp_df$BCR))
bcr.occ <- bcr.occ %>% mutate(bcr.code = 1:n()) %>%
  select(bcr.code, everything())
colnames(bcr.occ)[colnames(bcr.occ) == "unique(final_sp_df$BCR)"] <- "BCR"

bcr.occ <- bcr.detections
rownames(bcr.occ) <- bcr.occ$spp.id
bcr.occ <- bcr.occ[, -1]
colnames(bcr.occ) <- c(1, 2, 3, 4, 5)

#Write csv
#write.csv(bcr.occ, file = here("Data_BBS/Generated DFs/BCR.occ.csv"))

################################################
#####Creating Time and Segment JAGs DFs#########
################################################

###Create the time component for each segment###
time <- final_sp_df[, c("rteno.x", "rteno.id", 
                        "StartTime", "BCR", "bcr.id",
                        "EndTime", "Day", "Month", 
                        "ObsN", "ObsN.id", "site", 
                        "site.id", "Year", "year.id", "Segment")]
time$rt_yr <- paste0(time$rteno.x, "_", time$Year)
time <- time[!duplicated(time), ]
#segments <- final_sp_df[, c("rt_yr", "Segment")]

#Getting Segments and Time of Day Data Together 
time.occ <- time

#Time of Day for each segment calculated
#All route segments assumed to be surveyed at consistent interval
time.occ$StartTime <- as.character(time.occ$StartTime)
time.occ$EndTime <- as.character(time.occ$EndTime)

#Split military time up to hours and minutes
time.occ$StrtH <- str_sub(time.occ$StartTime, 1, -3)
time.occ$StrtMin <- str_sub(time.occ$StartTime, -2)
time.occ$EndH <- str_sub(time.occ$EndTime, 1, -3)
time.occ$EndMin <- str_sub(time.occ$EndTime, -2)

time.occ[, c("StrtH", "StrtMin", "EndH", "EndMin", "Segment")] <- lapply(time.occ[, c("StrtH", "StrtMin", "EndH", "EndMin", "Segment")], as.numeric)

#Convert to time elapsed since midnight
#Calculate duration of entire survey 
#Calculate duration of each segment
#Multiply duration  of each segment by segment #
#Add back to start time 
time.occ <- time.occ %>% mutate(StrtElapsed = (StrtH * 60) + StrtMin) %>%
  mutate(EndElapsed = (EndH * 60) + EndMin) %>% mutate(Duration = EndElapsed - StrtElapsed) %>%
  mutate(SegDur = Duration / 5) %>% mutate(TOD = (SegDur * Segment) + StrtElapsed)


#Calculate J-Date 
time.occ$Date <- paste0(time.occ$Month, "/", time.occ$Day, "/", time.occ$Year)
time.occ$Date <- as.Date(time.occ$Date, "%m/%d/%Y")
time.occ$OrdinalDate <- yday(time.occ$Date)

site.occ.df <- time.occ %>% select(rteno.x, rteno.id, site, 
                                   site.id, BCR, bcr.id, 
                                   Year, year.id, ObsN, 
                                   ObsN.id, TOD, 
                                   OrdinalDate)

#view
head( site.occ.df )

#Assigning codes for new observers

obs.change <- site.occ.df[, c("ObsN", "Year", "site", "rteno.x")]
obs.change <- obs.change %>% separate(site, into = c("rteno.x", "segment"), by = "_")
obs.change$site <- paste0(obs.change$rteno.x, "_", obs.change$segment)
obs.change$ObsID <- paste0(obs.change$ObsN, "_", obs.change$rteno.x, "_", obs.change$Year)
obs.change$rt_yr <- paste0(obs.change$rteno.x, "_", obs.change$Year, "_", obs.change$segment)
obs.change$Change <- NA
obs.change$ObsRt <- paste0(obs.change$rteno, "_", obs.change$ObsN)
obs.change <- transform(obs.change, RtObs.id = as.numeric(interaction(ObsRt, drop = T)))
obs.n.vect <- obs.change$RtObs.id
obs.n.vect.s <- obs.n.vect
obs.change <- obs.change[order(obs.change$rt_yr),]
#obs.change <- obs.change[!duplicated(obs.change),]
obs.change <- obs.change %>% mutate(obs.id.lag = lag(RtObs.id))

chang.fun <- function(x, y){
  x = x
  y = y
  ifelse(x == y, 0,1)
}

obs.change$Change <- chang.fun(obs.change$RtObs.id, obs.change$obs.id.lag)
obs.change$rt_yr <- paste0(obs.change$rteno.x, "_", obs.change$Year)
obs.change[1, 8] <- 1
obs.change <- obs.change %>% group_by(rt_yr) %>% mutate(ChangeD = first(Change)) %>%
  ungroup(obs.change)

obs.change <- obs.change %>% select(ObsID, ChangeD)
site.occ.df$ObsID <- paste0(site.occ.df$ObsN, "_", site.occ.df$rteno.x, "_", site.occ.df$Year)
site.occ.df <- merge(site.occ.df, obs.change, by = "ObsID")
site.occ.df <- site.occ.df[, -1]
site.occ.df <- site.occ.df[!duplicated(site.occ.df),]

years.occ <- time.occ %>% select(Year, year.id)
years.occ <-years.occ[!duplicated(years.occ),]

#write.csv
#write.csv(site.occ.df, file = here::here("Data_BBS/Generated DFs/Site.Occ.csv"), row.names = F)
#write.csv(years.occ, file = here::here("Data_BBS/Generated DFs/Years.Occ.csv"), row.names = F)


###########################################################################
################################################################################
##### save workspace and dataframes ###########################################
#save workspace
#save.image(file = here("R Workspace/DataPrep4OccModel5_15_2019.RData"))
#################### end of script 1###########################################

expit <- function(p){
  return(exp(p) / (exp(p) + 1))
}

logit <- function(p, e = 0.0005){
  p = ifelse(p == 1, p- e, p)
  p = ifelse(p == 0, p + e, p)
  return(log(p / (1 - p)))
}

######Standardizing functions######
######Generic Function Gelman's Suggestion######

standardise <- function(xmat, stdev = 1, marg = c(1, 2, 3)){
  mean.xmat = mean(as.vector(xmat), na.rm = T)
  sd.xmat = sd(as.vector(xmat), na.rm = T)
  std.xmat = apply(xmat, marg, function(x){
    (x - mean.xmat) / (stdev * sd.xmat)})
  retun(std.xmat)
}

scale2 <- function(x, na.rm = F) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


#####Standardizing Vectors#####
standardise.vector <- function(x, stdev = 1){(x - mean(x, na.rm = T))/
    (stdev * sd(x, na.rm = T))}

#####Ending Standardization Functions#####

###Load in DFs###

#final_sp_df <- read.csv(file = here::here("Data_BBS/Generated DFs/Final_Sp_Df_DetectionsOnly.csv"), header = T)
#site.occ.df <- read.csv(file = here::here("Data_BBS/Generated DFs/Site.Occ.csv"), header = T)
#spp.occ <- read.csv(file = here::here("Data_BBS/Generated DFs/Spp.Occ.csv"), header = T)
#Years.Occ <- read.csv(file = here::here("Data_BBS/Generated DFs/Years.Occ.csv"), header = T)
#bcr.occ <- read.csv(file = here::here("Data_BBS/Generated DFs/bcr.detection.raw.csv"), header = F)


#Pulling out one phylo grouping for test run
setDT(spp.occ)[, .(count = uniqueN(spp.id)), by = Phylo.V1]


###Remove column names from BCR
colnames(bcr.occ) <- NULL
bcr.occ <- as.matrix(bcr.occ)
###Fix two species body masses###
spp.occ$BodyMass[spp.occ$sci_name == "Sitta carolinensis"] <- 20
spp.occ$BodyMass[spp.occ$sci_name == "Auriparus flaviceps"] <- 8

#Standardize variables for analysis

##Long way to scale variables
#site.occ.df <- site.occ.df %>% mutate(TOD.scaled = scale(TOD, center = T, scale = T))
#site.occ.df$TOD.scaled <- as.numeric(site.occ.df$TOD.scaled)
#site.occ.df <- site.occ.df %>% mutate(Ordinal.scaled = scale(OrdinalDate, center = T, scale = T))
#site.occ.df$Ordinal.scaled <- as.numeric(site.occ.df$Ordinal.scaled)
#site.occ.df$ChangeD.scaled <- site.occ.df$ChangeD

site.occ.scaled <- site.occ.df %>% mutate_at(c("TOD", "OrdinalDate"), scale2)

site.occ.scaled$ChangeD.scaled <- site.occ.scaled$ChangeD

site.occ.scaled$ChangeD.scaled[site.occ.scaled$ChangeD.scaled == 0] <- -1

spp.occ <- spp.occ %>% mutate(Mass.scaled = scale(BodyMass, center = T, scale = T))
spp.occ$Mass.scaled <- as.numeric(spp.occ$Mass.scaled)

#Finished Standardizing and Centering

####Create Matrices
site.occ.ma <- site.occ.scaled %>% select(site, site.id, Year,
                                          year.id, TOD, OrdinalDate,
                                          ChangeD.scaled) %>%
  arrange(site.id)

##Site x Year Matrix with scaled.tod values and NAs for missing
#Time of day matrix
TOD.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum, 
                value.var = "TOD", fill = NA_real_, drop = F)
TOD.ma <- TOD.ma[, -1]
TOD.ma <- as.matrix(TOD.ma)
colnames(TOD.ma) <- NULL
TOD.ma[is.na(TOD.ma)] <- 0

#Ordinal Date Matrix
Ord.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum,
                value.var = "OrdinalDate", fill = NA_real_, drop = F)
Ord.ma <- Ord.ma[, -1]
Ord.ma <- as.matrix(Ord.ma)
colnames(Ord.ma) <- NULL
Ord.ma[is.na(Ord.ma)] <- 0

#1st year observer matrix (1 new obs, -1 same observer)
Obs.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum,
                value.var = "ChangeD.scaled", fill = NA_real_, drop = F)
Obs.ma <- Obs.ma[, -1]
Obs.ma <- as.matrix(Obs.ma)
colnames(Obs.ma) <- NULL
Obs.ma[is.na(Obs.ma)] <- 0

#Route, Segment, and BCR DF (jdf)
jdf <- site.occ.df[, c("rteno.x", "rteno.id", "site",
                       "site.id", "BCR", "bcr.id")]
jdf <- jdf[!duplicated(jdf),]

jdf <- jdf %>% arrange(site.id)


#write.csv(jdf, file = here::here("Data_BBS/Generated DFs/jdf.csv"))
#write.csv(TOD.ma, file = here::here("Data_BBS/Generated DFs/TOD_ma.csv"))
#write.csv(Ord.ma, file = here::here("Data_BBS/Generated DFs/Ord_ma.csv"))
#write.csv(Obs.ma, file = here::here("Data_BBS/Generated DFs/Obs_ma.csv"))
#write.csv(spp.occ, file = here::here("Data_BBS/Generated DFs/spp_occ.csv"))
#write.csv(bcr.occ, file = here::here("Data_BBS/Generated DFs/bcr_occ.csv"))

#####Import species specific data#####
#Creation of the complete ydf matrix
#alternative approach ## 
#head( spp.df )
head(  final_sp_df )
head( site.occ.df )
glimpse( site.occ.df )
head( spp.occ)

#check for duplicates
site.occ.df[ duplicated( site.occ.df ), ] #none!

#create reduced dataframe from final_sp_df
spp.df <- final_sp_df %>% dplyr::select( sci_name, site, Year, Detected )
#view
head( spp.df ); dim( spp.df )
#remove zero detections:
spp.df <- spp.df %>% dplyr::filter( Detected == 1 )
#view
head( spp.df ); dim( spp.df )
glimpse( spp.df )
#append species data
#first relabel first column
colnames( spp.df )[ 1 ] <- "Species"
colnames( spp.occ )[ 1 ] <- "Species"
spp.df <- left_join( spp.df, spp.occ, by = "Species" )
#check
#view
head( spp.df ); dim( spp.df ) #it didn't add columns..hooray!
#check for duplicates:
spp.df[ which( (spp.df$spp.id == 1) & ( spp.df$site == '2141_4') ),  ]
#append siteXyear info #ensure we turn all.y=F so that sitesXyear unsurveyed are not
# added
spp.df <- left_join( spp.df, site.occ.df, by = c( "site", "Year" ), all.x=T, all.y=F )
#view
tail( spp.df ); dim( spp.df ) #row numbers stayed the same
#did it add rows
sum( spp.df$Detected ) 
#check for duplicates
spp.df[ which( (spp.df$spp.id == 1) & ( spp.df$site.id == 12) ),  ]


#remove duplicates for species for now #you need to not do this with the new
#final df. #####
spp.df <- spp.df %>% group_by( spp.id, site.id, year.id ) %>%
  slice( 1 ) #only keeps one record for each speciesXsiteXyear
tail( spp.df ); dim( spp.df )
#removed ~200,000 records

# setdimensions ###
#total number of species
S <- max(spp.df$spp.id)
#total number of segments
J <- max(spp.df$site.id)
#total number of routes
M <- max(spp.df$rteno.id)
#total number of sampling years
K <- max(spp.df$year.id)
#number of species in each BCR vector
B <- colSums(bcr.occ)
#number of phylogenetic groups 
G <- length(unique(spp.occ$Phylo.V1.code))


### working out missing sampling years for a given segment:
#we will work from site.occ.df which has all details of which segments were surveyed 
# which year:
#view
tail( site.occ.df ); dim( site.occ.df )
site.occ.df[ which( site.occ.df$site.id == 1), ]
spp.df[ which( spp.df$site.id == 1), ]
glimpse( site.occ.df )
#use it to create complete dataframe
JKdf <- site.occ.df %>% dplyr::select( site.id, year.id )
#check for duplicates
JKdf[ duplicated( JKdf ), ] #some duplicates present!!!!
#remove
#JKdf <- JKdf[ !duplicated( JKdf ), ]
#convert year code to factor
JKdf$year.id <- as.factor( JKdf$year.id )
#check
levels( JKdf$year.id )
#convert site code to factor
JKdf$site.id <- as.factor( JKdf$site.id )
#check
levels( JKdf$site.id )
#add surveyed column
JKdf$surveyed <- 1 
#view
head( JKdf ); dim( JKdf )
#add missing combinations
JKdf <- JKdf %>% tidyr::complete( site.id, year.id )
#view
head( JKdf ); dim( JKdf )
#dimensions should equal: 
J*K
#turn to wide format
JKmat <- tidyr::spread( JKdf, key = year.id, value = surveyed )
#it has to have J rows and K columns
head( JKmat); dim( JKmat )
#convert to matrix:
JKmat <- as.matrix( JKmat[ ,2:dim( JKmat )[2] ] )
# #remove column names from matrix
colnames( JKmat ) <- NULL
#view
head( JKmat); dim( JKmat )
#define how many segments were surveyed each year:
surveyedJ <- colSums( JKmat, na.rm = TRUE )


#create indicator for calculating summary statistics
JKsurv <- JKmat
JKsurv[ is.na(JKsurv) ] <- 0 
#check
JKmat[1:5,1:10]
JKsurv[1:5,1:10]


####### now we create the ydf ##### 
#create observations dataframe
ydf <- array( 0, dim = c(S, J, K) )
ydf[ 1, , ] * JKmat
# assigned 1 when species was detected on given year and route
for( i in 1:dim( spp.df )[1] ){
  ydf[ as.numeric( spp.df[i,'spp.id'] ), as.numeric(spp.df[i,'site.id']),
       as.numeric( spp.df[i,'year.id'] ) ] <- 1 #as.numeric( spp.df[i, 'Detected'] )
}

# #check that it worked for species 1:
table( ydf[4,,] )
table( spp.df$Detected[ which( spp.df$spp.id == 4) ] )  #77 
spp.df[ which( spp.df$spp.id == 4), ]
ydf[1,8,32]

#add missing values 
for( i in 1:S ){
  ydf[ i, , ] <- ydf[ i, , ] * JKmat
}

ydf[6, 26, ]



#####Workspace for testing subsetted groupings######
rm(list = setdiff(ls(), c("ydf", "jdf", "bcr.occ", 
                          "spp.occ", "TOD.ma", "Ord.ma",
                          "Obs.ma", "JKmat")))
save.image(here::here("R Workspace/CardinalsFinchesSubgroup.RData"))
