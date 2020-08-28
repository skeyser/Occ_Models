#######################################################################
#############Script for calculating diversity metrics##################
#######and pairing avian data with the LULC and climate data###########
#######################for model running###############################
##################Created by: Spencer R Keyser#########################
#######################################################################

#Start from line 748 to get at what the issue is for the sites
#Today it seemed ok? Need to double check 
#Should hve climate data for 84 rts and 274 sites


##Clean Global Environment  
rm(list = ls())

#Package Loading
#install.packages("pacman")
#library("pacman")
pacman::p_load("here", "tidyverse", "reshape2", "vegan", "data.table", "cowplot", "lme4", "sjPlot", 
               "sjstats", "car", "betapart")

#Packages to call in the function without loading
#pacman::p_load(MASS, effects, nlme, lsr, lmerTest, MUMIn)


# #Load in cleaned spp x site DF
# final_sp_df <- read.csv(here::here("Data_BBS/Generated DFs/Final_Sp_Df_DetectionsOnly.csv"))
# bbs.combined <- read.csv(here::here("Data_BBS/Generated DFs/BBS_Total_DF.csv"))

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

#Pull out wetland sites lat/long
site.coords <- rts_gom[rts_gom$rteno %in% wetland20.site,]

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


#Rename phylo in the complete DF for grouping species
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Sturnidae"] <- "Sturnidae/Turdidae/Mimidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Turdidae"] <- "Sturnidae/Turdidae/Mimidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Mimidae"] <- "Sturnidae/Turdidae/Mimidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Alaudidae"] <- "Remizidae/Paridae/Alaudidae/Hirundinidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Paridae"] <- "Remizidae/Paridae/Alaudidae/Hirundinidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Remizidae"] <- "Remizidae/Paridae/Alaudidae/Hirundinidae"
final_sp_df$Family_Latin[final_sp_df$Family_Latin == "Hirundinidae"] <- "Remizidae/Paridae/Alaudidae/Hirundinidae"
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
final_sp_df$Order[final_sp_df$Order == "Cuculiformes"] <- "Cuculiformes/Columbiformes"
final_sp_df$Order[final_sp_df$Order == "Columbiformes"] <- "Cuculiformes/Columbiformes"
final_sp_df$Order[final_sp_df$Order == "Anseriformes"] <- "Anseriformes/Galliformes"
final_sp_df$Order[final_sp_df$Order == "Galliformes"] <- "Anseriformes/Galliformes"


####Remove BCR 19#####
final_sp_df <- final_sp_df[!(final_sp_df$BCR == 19),]

########################################################
############From this point on every DF generated#######
####should begin from the final_sp_df for consistency###
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

final_sp_df <- final_sp_df %>% dplyr::select(sci_name, English, spp.id, Detected, Phylo.V1, Phylo.V1.code,
                                             Phylo.V2, Phylo.V2.code, rteno.x, rteno.id, rt_yr, Segment, 
                                             Year, year.id, Day, Month, BCR, bcr.id, site, site.id, 
                                             ObsN, ObsN.id, BodyMass, Nocturnal, Category, Inverts,
                                             Mammal.Bird, Herps, Fish, Unknown_Verts, Fruit,
                                             Nect, Seed, Other_plant, ForStrat_watbelowsurf, 
                                             ForStrat_wataroundsurf, ForStrat_ground, 
                                             ForStrat_understory, ForStrat_midhigh, 
                                             ForStrat_canopy, ForStrat_aerial, 
                                             PelagicSpecialist, StartTime, EndTime)  


#Clean the workspace except for the files needed
#rm(list = setdiff(ls(), c()))



#Pull Out the Lat/Long
site.list <- unique(final_sp_df$rteno.x)
rt.xy <- rts_gom[rts_gom$rteno %in% site.list,]
rt.xy <- rt.xy[rt.xy$BCR != 19, ]
rt.xy <- rt.xy %>% dplyr::select(rtename, rteno, countrynum, Latitude,
                                 Longitude)

#Create a spp by site matrix
#Keep only the columns we need
#Create site x year id
final_sp_df$unique_id <- paste0(final_sp_df$site, "_", final_sp_df$Year)

bbs_simple <- final_sp_df[,c("sci_name", "unique_id", "Year", "site")]
bbs_simple$Detected <- 1

bbs_simple <- merge(bbs_simple, Bins, by = "Year") #Bring in the bins 
bbs_simple$unique_id <- paste0(bbs_simple$site, "_", bbs_simple$Year) #Moving back to year to year

#bbs_plot <- bbs_total[, c("unique_id", "statenum")]
#bbs_plot <- as.data.frame(bbs_plot)

#Cast Data for calculating the sums for site diversity 
bbs_cast <- dcast(data = bbs_simple, formula = unique_id ~ sci_name, fun.aggregate = sum)

#Remove column 1 after setting rownames
rownames(bbs_cast) <- bbs_cast$unique_id
bbs_cast <- bbs_cast[,-1]
bbs_cast[bbs_cast > 1] <- 1

#Calculate site diversity by summming across sites 
Site_div <- rowSums(bbs_cast)
Site_div <- as.data.frame(Site_div)
Site_div$unique_id <- rownames(Site_div)
rownames(Site_div) <- NULL


#check names are the same so we can just cbind them
bbs_total <- final_sp_df[!duplicated(final_sp_df),] #Remove any duplicated rows

bbs_total <- bbs_total[bbs_total$Year >= 1980,]
bbs_total$Yr_bin <- 2
bbs_total$Yr_bin[bbs_total$Year >= 1985 & bbs_total$Year <= 1989] <- 3
bbs_total$Yr_bin[bbs_total$Year >= 1990 & bbs_total$Year <= 1994] <- 4
bbs_total$Yr_bin[bbs_total$Year >= 1995 & bbs_total$Year <= 1999] <- 5
bbs_total$Yr_bin[bbs_total$Year >= 2000 & bbs_total$Year <= 2004] <- 6
bbs_total$Yr_bin[bbs_total$Year >= 2005 & bbs_total$Year <= 2009] <- 7
bbs_total$Yr_bin[bbs_total$Year >= 2010 & bbs_total$Year <= 2014] <- 8
bbs_total$Yr_bin[bbs_total$Year >= 2015 & bbs_total$Year <= 2018] <- 9

bbs_total$rt_yr <- paste0(bbs_total$rteno.x, "_", bbs_total$Year) #Make this all based on the yr_bin 

bbs_total <- bbs_total %>% dplyr::select(rteno.x, rt_yr, Segment, Year,
                                         BCR, site, unique_id) #Selection variables 

bbs_total$unique_id <- paste0(bbs_total$site, "_", bbs_total$Year) 

bbs_total <- bbs_total[!duplicated(bbs_total),] #Remove any duplicated rows 
bbs_total <- arrange(bbs_total, desc(unique_id)) #Order bbs_total
Site_div <- arrange(Site_div, desc(unique_id)) #Order Site_dive
bbs_total <- cbind(bbs_total, Site_div) #cbind two dfs
bbs_total <- bbs_total[, -9]

bbs_total <- merge(bbs_total, rt.xy, by.x = "rteno.x", by.y = "rteno")

#Creating DF with only Wetland Species for Calculating Beta Div
#Load in the Army Corps Data
#Bring in AOU for merging 
aou <- species[, c("AOU", "Genus", "Species")]
aou$species <- paste0(aou$Genus, " ", aou$Species)
aou <- aou[, c("AOU", "species")]

#Fixing some annoying name issues
aou[aou$species == "Antigone canadensis", 2] <- "Grus canadensis"
aou[aou$species == "Spinus psaltria", 2] <- "Carduelis psaltria"
aou[aou$species == "Hydroprogne caspia", 2] <- "Sterna caspia" 

#Read in the nesting DF
nest <- read.csv(here::here("Functional_Traits/Army_Corps_Bird_Data.csv"), stringsAsFactors = F)

#Create a list of birds from the complete DF
bbs.spp <- data.frame(species = unique(final_sp_df$sci_name))
bbs.spp$species <- as.character(bbs.spp$species)

#Make list of birds from the nesting DF
nest.spp <- nest[, c("AOU", "SCIENTIFIC_NAME")]
nest.spp$SCIENTIFIC_NAME <- as.character(nest.spp$SCIENTIFIC_NAME)

#Merge the species list with the AOU for future merging
bbs.spp <- merge(bbs.spp, aou, by = "species")

#Northern Flicker wasn't matching up 
bbs.spp[bbs.spp$AOU == 4123, 2] <- 4120


bbs_nest <- merge(bbs.spp, nest, by = "AOU")
bbs_nest <- bbs_nest %>% dplyr::select(-c("SCIENTIFIC_NAME"))

#Get habitats in with the full spp_df
bbs_nest <- merge(final_sp_df, bbs_nest, by.x = "sci_name", by.y = "species")

#Create vector of habitats we want
habs <- c("Coastal", "Swamp/Marsh", "Shoreline")

#Filter the DF with the specific habitats 
bbs_nest <- bbs_nest %>% filter(grepl(paste(habs, collapse = "|"), HABITAT))


#Mean SR for each route
bbs_total$site <- as.character(bbs_total$site)
Site.means <- bbs_total %>% group_by(site) %>%
  summarise(mean.alpha = mean(Site_div)) %>%
  ungroup(bbs_total)

Site.means.bin <- bbs_total %>% group_by(site, Yr_bin) %>%
  summarise(alpha.bin = mean(Site_div)) %>%
  ungroup(bbs_total)

Site.means.bin$unique_id <- paste0(Site.means.bin$site, "_", Site.means.bin$Yr_bin)
Site.means.bin <- Site.means.bin[, c("unique_id", "alpha.bin")]


#Bring in the MCMC DFs
beta.mcmc <- read.csv(file = here::here("Data_BBS/Generated DFs/beta_means.csv"))
beta.mcmc <- beta.mcmc[, -1]
beta.mcmc <- beta.mcmc %>% rename(unique_id = Sites)
beta.mcmc$unique_id <- as.character(beta.mcmc$unique_id)

sr.mcmc <- read.csv(file = here::here("Data_BBS/Generated DFs/sr_means.csv"))
sr.mcmc <- sr.mcmc[, -1]
sr.mcmc <- sr.mcmc %>% rename(unique_id = Sites)
sr.mcmc$unique_id <- as.character(sr.mcmc$unique_id)


#Make bins DF
Bins <- data.frame(matrix(ncol = 2, nrow = 40))
l <- c("Year", "Yr_bin")
colnames(Bins) <- l
Bins$Year <- 1980:2019
Bins$Yr_bin <- rep(2:9, times = 1, each = 5)

#Read in Occ Data 50 cutoff
occ <- read.csv(here::here("Data_BBS/Generated DFs/occ50.csv"), stringsAsFactors = F)

# #Look at structure of the occ model
# occ.str <- occ %>%
#   select(site, Yr_bin, Year) %>%
#   group_by(site, Year) %>%
#   distinct() %>%
#   add_count()

#Merge it with the year bins 
#occ <- merge(occ, Bins, by = "Year")
occ$unique_id <- paste0(occ$site, "_", occ$Year)
#occ$unique_id <- paste0(occ$site, "_", occ$Year)

#Making a occ DF with only wetland birds
occ.spp <- data.frame(species = unique(occ$Species))
occ.spp$species <- as.character(occ.spp$species)

#Merge the species list with the AOU for future merging
occ.spp <- merge(occ.spp, aou, by = "species")

#Northern Flicker wasn't matching up 
bbs.spp[bbs.spp$AOU == 4123, 2] <- 4120

#Bring in the AOU for cleaner merging 
#Names of birds constantly changing 
occ_nest <- merge(occ.spp, nest, by = "AOU")
occ_nest <- occ_nest %>% dplyr::select(-c("SCIENTIFIC_NAME"))

#Get habitats in with the full spp_df
occ_nest <- merge(occ, occ_nest, by.x = "Species", by.y = "species")
occ.wet <- occ_nest %>% filter(grepl(paste(habs, collapse = "|"), HABITAT))

#Cast spp x site_yr.bin
bbs_cast2 <- dcast(occ, unique_id ~ Species, value.var = "Occupancy", fun.aggregate = sum)
rownames(bbs_cast2) <- bbs_cast2$unique_id
bbs_cast2 <- bbs_cast2[, -1]

#Cast spp x site_year for wetland spp
wet_cast <- dcast(occ.wet, unique_id ~ Species, value.var = "Occupancy", fun.aggregate = sum)
rownames(wet_cast) <- wet_cast$unique_id
wet_cast <- wet_cast[, -1]

#Make all detections 1
bbs_cast2[bbs_cast2 >= 1] <- 1
wet_cast[wet_cast >= 1] <- 1

#Calculate SR for 50 cutoff
occSR <- data.frame(unique_id = rownames(bbs_cast2), SR50 = NA, SR.wet = NA)
occSR <- occSR[order(occSR$unique_id), ]
bbs_cast2 <- bbs_cast2[order(rownames(bbs_cast2)), ]
wet_cast <- wet_cast[order(rownames(wet_cast)), ]
occSR$SR50 <- rowSums(bbs_cast2)
occSR$SR.wet <- rowSums(wet_cast)
occSR$unique_id <- as.character(occSR$unique_id)

#Initializing DFs for loops 
n.sites <- length(unique(bbs_total$site))
site.list <- as.character(unique(bbs_total$site))

slopes_sites <- data.frame(site.list, slope = NA, slope.mcmc = NA, slope50 = NA, slope.wet = NA)


sr.mcmc <- sr.mcmc %>%
  separate(unique_id, into = c("route", "segment", "Year")) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(Bins, by = "Year") %>%
  mutate(unique_id = paste0(route, "_", segment, "_", Year)) %>%
  group_by(unique_id) %>%
  summarise(Mean_SR = mean(Mean_SR),
            Sd_SR = mean(Sd_SR))
  


bbs_div_means <- bbs_total
bbs_div_means <- bbs_div_means %>% left_join(sr.mcmc, by = "unique_id") %>% left_join(occSR, by = "unique_id")
bbs_div_means <- bbs_div_means %>% rename(SR_MCMC = Mean_SR, SD_MCMC = Sd_SR)
colnames(bbs_div_means)[colnames(bbs_div_means) == "Site_div.y"] <- "Site_div" 
#bbs_div_means <- bbs_div_means %>% group_by(site, Yr_bin)%>% summarise(mean.div = mean(Site_div))

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_total_temp <- bbs_div_means[bbs_div_means$site == site.temp,]
  if (nrow(site_total_temp) > 1){
    lm.temp <- lm(site_total_temp$Site_div ~ site_total_temp$Year)
    lm.temp2 <- lm(site_total_temp$SR_MCMC ~ site_total_temp$Year)
    lm.temp3 <- lm(site_total_temp$SR50 ~ site_total_temp$Year)
    lm.temp4 <- lm(site_total_temp$SR.wet ~ site_total_temp$Year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slope.temp2 <- summary(lm.temp2)$coefficients[2,1]
    slope.temp3 <- summary(lm.temp3)$coefficients[2,1]
    slope.temp4 <- summary(lm.temp4)$coefficients[2,1]
    slopes_sites[i, 2] <- slope.temp
    slopes_sites[i, 3] <- slope.temp2
    slopes_sites[i, 4] <- slope.temp3
    slopes_sites[i, 5] <- slope.temp4
  }
}

hist(slopes_sites$slope, breaks = 50)
hist(slopes_sites$slope.mcmc, breaks = 50)
hist(slopes_sites$slope50, breaks = 50)
hist(slopes_sites$slope.wet, breaks = 50)

slopes_sites$site <- as.character(slopes_sites$site.list)
slopes_sites <- slopes_sites %>% dplyr::select(site, slope, slope.mcmc, slope50, slope.wet)

#Mean alpha at the site level 
bbs_div_means <- bbs_div_means[, unique(colnames(bbs_div_means))]
site_means <- bbs_div_means %>% dplyr::select(Site_div, SR_MCMC, SR50, SR.wet, site)
site_means <- site_means %>% group_by(site) %>% summarise(MeanDiv = mean(Site_div), MeanDivMCMC = mean(SR_MCMC),
                                                          MeanDiv50 = mean(SR50), MeanDivWet = mean(SR.wet))

#Alpha diversity at the site x yr level
site_alpha <- bbs_div_means %>% dplyr::select(Site_div, SR50, SR_MCMC, SR.wet, SD_MCMC, unique_id)

#Mean alpha and alpha slopes at the site level 
site_means2 <- merge(site_means, slopes_sites, by = "site")


mean(site_means2$slope, na.rm = T)
median(site_means2$slope, na.rm = T)
mean(site_means2$slope.mcmc, na.rm = T)
median(site_means2$slope.mcmc, na.rm = T)
mean(site_means2$slope50, na.rm = T)
median(site_means2$slope50, na.rm = T)

############################################################
#############Beta diversity preparation#####################
############################################################

#Binning the species specific dataset for community bins
bbs_simple <- bbs_simple[bbs_simple$Year >= 1980,]
bbs_simple$Yr_bin <- 2
bbs_simple$Yr_bin[bbs_simple$Year >= 1985 & bbs_simple$Year <= 1989] <- 3
bbs_simple$Yr_bin[bbs_simple$Year >= 1990 & bbs_simple$Year <= 1994] <- 4
bbs_simple$Yr_bin[bbs_simple$Year >= 1995 & bbs_simple$Year <= 1999] <- 5
bbs_simple$Yr_bin[bbs_simple$Year >= 2000 & bbs_simple$Year <= 2004] <- 6
bbs_simple$Yr_bin[bbs_simple$Year >= 2005 & bbs_simple$Year <= 2009] <- 7
bbs_simple$Yr_bin[bbs_simple$Year >= 2010 & bbs_simple$Year <= 2014] <- 8
bbs_simple$Yr_bin[bbs_simple$Year >= 2015 & bbs_simple$Year <= 2018] <- 9

#Calculate community means per yr_bin
bbs_simple$unique_id <- paste0(bbs_simple$site, "_", bbs_simple$Year)

#Cast spp x site_yr.bin
bbs_cast <- dcast(bbs_simple, unique_id ~ sci_name, value.var = "Detected", fun.aggregate = sum)
rownames(bbs_cast) <- bbs_cast$unique_id
bbs_cast <- bbs_cast[, -1]

#Make all detections 1
bbs_cast[bbs_cast >= 1] <- 1

#mean temporal beta diversity by site
mean.betas <- data.frame(site.list, mean.beta = NA)

#Create matrix for storing Betas 
beta.matrix <- data.frame(unique_id = bbs_simple$unique_id, beta.turn = NA, beta.nest = NA, beta.jac = NA, 
                          beta.wet.turn = NA, beta.wet.nest = NA, beta.wet.jac = NA,
                          beta50.turn = NA, beta50.nest = NA, beta50.jac = NA)

beta.matrix$unique_id <- as.character(beta.matrix$unique_id)

#Remove duplites
beta.matrix <- beta.matrix[!duplicated(beta.matrix),]

#Loop through the sites so that we calculate beta based on year bins
#Change the loop to do year to year beta 

#Loop below was created for the year bin to year bin calculation. NOT USED!

# for (n in 1:n.sites){
#   site.temp <- site.list[n]
#   #Find all the years for which that site was surveys
#   yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Yr_bin"])
#   yrs.temp <- yrs.temp[yrs.temp >= 1]
#   #only consider sites that were observed in more than 5 years 
#   if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin 
#     #Create a baseline year from years temp
#     #Sorting
#     yrs.temp <- yrs.temp[order(yrs.temp)]
#     #find the first 5 years 
#     for (y in 1:length(yrs.temp)){
#       first.yr <- yrs.temp[y]
#       if(first.yr == max(yrs.temp)){
#         next
#         }
#       
#       second.yr <- yrs.temp[y + 1]
#       print(paste0("Change from ", first.yr, " to ", second.yr))
#       site.yrs.first <- paste0(site.temp, "_", first.yr)
#       site.yrs.second <- paste0(site.temp, "_", second.yr)
#       bbs_cast_first <- bbs_cast[rownames(bbs_cast) %in% site.yrs.first,]
#       bbs_cast_second <- bbs_cast[rownames(bbs_cast) %in% site.yrs.second,]
#       bbs_cast_tmp <- rbind(bbs_cast_first, bbs_cast_second)
#       beta.temp <- betapart::beta.multi(bbs_cast_tmp, index.family = "jaccard")
#       #beta.matrix$beta[beta.matrix$unique_id == other.site.temp] <- beta.temp
#       beta.matrix$beta.turn[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[1]) #store it in df
#       beta.matrix$beta.nest[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[2])
#       beta.matrix$beta.jac[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[3])
#       
#     }}}


for (n in 1:n.sites){
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Year"])
  yrs.temp <- yrs.temp[yrs.temp >= 1]
  #only consider sites that were observed in more than 5 years
  if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin
    #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1])] #pulls first year bin
    print(paste0("Starting Year for ", site.temp, " is ", first.yrs))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years
    bbs_cast_first <- bbs_cast[rownames(bbs_cast) %in% site.yrs.first,]
    first.comm.temp <- bbs_cast_first
    #get remaining years
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1])]
    #Now loop the remaining years
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y] #pull out another yr
      other.site.temp <- paste0(site.temp, "_", other.yr.temp) #make the unique id
      bbs_cast_other <- bbs_cast[rownames(bbs_cast) == other.site.temp,] #pull the comm data
      bbs_cast_other <- rbind(bbs_cast_other, first.comm.temp) #put the current yr and baseline
      #beta.temp <- vegdist(sqrt(sqrt(bbs_cast_other)), method = "jaccard") #calc jaccard disim
      beta.temp <- betapart::beta.multi(bbs_cast_other, index.family = "jaccard")
      #beta.matrix$beta[beta.matrix$unique_id == other.site.temp] <- beta.temp
      beta.matrix$beta.turn[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[1]) #store it in df
      beta.matrix$beta.nest[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[2])
      beta.matrix$beta.jac[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[3])
    }
  }}

#Loop through the occ50 sites so that we calculate beta based on year bins
#Change the loop to do year to year beta 
# for (n in 1:n.sites){
#   site.temp <- site.list[n]
#   #Find all the years for which that site was surveys
#   yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Year"])
#   yrs.temp <- yrs.temp[yrs.temp >= 1]
#   #only consider sites that were observed in more than 5 years 
#   if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin 
#     #Create a baseline year from years temp
#     #Sorting
#     yrs.temp <- yrs.temp[order(yrs.temp)]
#     #find the first 5 years 
#     for (y in 1:length(yrs.temp)){
#       first.yr <- yrs.temp[y]
#       if(first.yr == max(yrs.temp)){
#         next
#       }
#       
#       second.yr <- yrs.temp[y + 1]
#       print(paste0("Change from ", first.yr, " to ", second.yr))
#       site.yrs.first <- paste0(site.temp, "_", first.yr)
#       site.yrs.second <- paste0(site.temp, "_", second.yr)
#       bbs_cast_first <- bbs_cast2[rownames(bbs_cast2) %in% site.yrs.first,]
#       bbs_cast_second <- bbs_cast2[rownames(bbs_cast2) %in% site.yrs.second,]
#       bbs_cast_tmp <- rbind(bbs_cast_first, bbs_cast_second)
#       beta.temp <- betapart::beta.multi(bbs_cast_tmp, index.family = "jaccard")
#       #beta.matrix$beta[beta.matrix$unique_id == other.site.temp] <- beta.temp
#       beta.matrix$beta50.turn[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[1]) #store it in df
#       beta.matrix$beta50.nest[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[2])
#       beta.matrix$beta50.jac[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[3])
#       
#     }}}

# #Loop for Occ50 Matrix
for (n in 1:n.sites){
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(occ[which(occ$site == site.temp),"Year"])
  yrs.temp <- yrs.temp[yrs.temp >= 1]
  #only consider sites that were observed in more than 5 years
  if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin
    #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1])] #pulls first year bin
    print(paste0("Starting Year for ", site.temp, " is ", first.yrs))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years
    bbs_cast_first <- bbs_cast2[rownames(bbs_cast2) %in% site.yrs.first,]
    first.comm.temp <- bbs_cast_first
    #get remaining years
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1])]
    #Now loop the remaining years
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y] #pull out another yr
      other.site.temp <- paste0(site.temp, "_", other.yr.temp) #make the unique id
      bbs_cast_other <- bbs_cast2[rownames(bbs_cast2) == other.site.temp,] #pull the comm data
      bbs_cast_other <- rbind(bbs_cast_other, first.comm.temp) #put the current yr and baseline
      beta.temp <- betapart::beta.multi(bbs_cast_other, index.family = "jaccard")
      beta.matrix$beta50.turn[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[1]) #store it in df
      beta.matrix$beta50.nest[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[2])
      beta.matrix$beta50.jac[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[3])
      #beta.temp <- vegdist(sqrt(sqrt(bbs_cast_other)), method = "jaccard") #calc jaccard disim
      #beta.matrix$beta50[beta.matrix$unique_id == other.site.temp] <- beta.temp #store it in df
    }
  }}
# 

#Loop through the wetland only sites so that we calculate beta based on year bins
#Change the loop to do year to year beta 
# for (n in 1:n.sites){
#   site.temp <- site.list[n]
#   #Find all the years for which that site was surveys
#   yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Yr_bin"])
#   yrs.temp <- yrs.temp[yrs.temp >= 1]
#   #only consider sites that were observed in more than 5 years 
#   if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin 
#     #Create a baseline year from years temp
#     #Sorting
#     yrs.temp <- yrs.temp[order(yrs.temp)]
#     #find the first 5 years 
#     for (y in 1:length(yrs.temp)){
#       first.yr <- yrs.temp[y]
#       if(first.yr == max(yrs.temp)){
#         next
#       }
#       
#       second.yr <- yrs.temp[y + 1]
#       print(paste0("Change from ", first.yr, " to ", second.yr))
#       site.yrs.first <- paste0(site.temp, "_", first.yr)
#       site.yrs.second <- paste0(site.temp, "_", second.yr)
#       bbs_cast_first <- wet_cast[rownames(wet_cast) %in% site.yrs.first,]
#       bbs_cast_second <- wet_cast[rownames(wet_cast) %in% site.yrs.second,]
#       bbs_cast_tmp <- rbind(bbs_cast_first, bbs_cast_second)
#       beta.temp <- betapart::beta.multi(bbs_cast_tmp, index.family = "jaccard")
#       #beta.matrix$beta[beta.matrix$unique_id == other.site.temp] <- beta.temp
#       beta.matrix$beta.wet.turn[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[1]) #store it in df
#       beta.matrix$beta.wet.nest[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[2])
#       beta.matrix$beta.wet.jac[beta.matrix$unique_id == site.yrs.second] <- as.numeric(beta.temp[3])
#       
#     }}}
# 
# #Loop for wetland spp
for (n in 1:n.sites){
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Year"])
  yrs.temp <- yrs.temp[yrs.temp >= 1]
  #only consider sites that were observed in more than 5 years
  if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin
    #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1])] #pulls first year bin
    print(paste0("Starting Year for ", site.temp, " is ", first.yrs))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years
    wet_cast_first <- wet_cast[rownames(wet_cast) %in% site.yrs.first,]
    first.comm.temp <- wet_cast_first
    #get remaining years
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1])]
    #Now loop the remaining years
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y] #pull out another yr
      other.site.temp <- paste0(site.temp, "_", other.yr.temp) #make the unique id
      wet_cast_other <- wet_cast[rownames(wet_cast) == other.site.temp,] #pull the comm data
      wet_cast_other <- rbind(wet_cast_other, first.comm.temp) #put the current yr and baseline
      beta.temp <- betapart::beta.multi(wet_cast_other, index.family = "jaccard")
      beta.matrix$beta.wet.turn[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[1]) #store it in df
      beta.matrix$beta.wet.nest[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[2])
      beta.matrix$beta.wet.jac[beta.matrix$unique_id == other.site.temp] <- as.numeric(beta.temp[3])
      #beta.temp <- vegdist(sqrt(sqrt(wet_cast_other)), method = "jaccard") #calc jaccard disim
      #beta.matrix$beta.wet[beta.matrix$unique_id == other.site.temp] <- beta.temp #store it in df
    }
  }}

alpha_merge <- site_alpha %>%
  dplyr::select(unique_id, Site_div, SR50, SR.wet) %>%
  separate(unique_id, c("Route", "Segment", "Year"), sep = "_")


beta.mat <- beta.matrix %>% 
  separate(unique_id, c("Route", "Segment", "Year"), sep = "_") %>%
  mutate(Site = paste0(Route, "_", Segment)) %>%
  #group_by(Site) %>%
  #mutate(FirstYr = min(Year),
  #       LastYr = max(Year)) %>%
  complete(Site, nesting(Year))

beta.yrs <- beta.matrix %>%
  separate(unique_id, c("Route", "Segment", "Year"), sep = "_") %>%
  mutate(Site = paste0(Route, "_", Segment)) %>%
  group_by(Site) %>%
  mutate(FirstYr = min(Year),
         LastYr = max(Year)) %>%
  distinct(Site, FirstYr, LastYr) %>%
  left_join(beta.mat) %>%
  filter(Year >= FirstYr & Year <= LastYr) %>%
  fill(Route, Segment) %>%
  left_join(alpha_merge, by = c("Route", "Segment", "Year"))
  
beta.matrix <- beta.yrs 



test <- beta.yrs %>%
  group_by(Site) %>%
  #add_count() %>%
  slice(2:n) %>%
  filter(is.na(beta.turn)) %>%
  add_count() %>%
  select(Site, n) %>%
  distinct() %>%
  ggplot(aes(x = n)) + 
  geom_histogram()
  

beta.acf <- beta.matrix %>% 
  select(unique_id, beta.jac) %>%
  separate(unique_id, into = c("Route", "Seg", "Year"), sep = "_", remove = F) %>%
  mutate(Site = paste0(Route, "_", Seg)) %>%
  mutate_at("Year", as.integer) %>%
  group_by(Site) %>%
  mutate(start_yr = min(Year),
         end_yr = max(Year)) %>%
  complete(Year = start_yr:end_yr) %>%
  ungroup(Site) %>%
  mutate(unique_id = paste0(Site, "_", Year)) %>%
  select(-Site, -Seg, -Year, -Route, -start_yr, -end_yr) %>%
  separate(unique_id, into = c("Route", "Seg", "Year"), sep = "_")

grouped_acf_values <- beta.acf %>%
  
  select(Site, beta.jac) %>%
  filter(!is.na(beta.jac)) %>%
  group_by(Site) %>%
  add_count(Site) %>%
  filter(n >= 2) %>%
  ungroup() %>%
  select(Site, beta.jac) %>%
  tidyr::nest(-Site) %>%
  dplyr::mutate(acf_results = purrr::map(data, ~ acf(.x$beta.jac,
                                                     type = "correlation", plot = T)),
                acf_values = purrr::map(acf_results, ~ drop(.x$acf))) %>%
  tidyr::unnest(acf_values) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(lag = seq(0, n() - 1))

head(grouped_acf_values)

#Quick plot trends through time 
beta.plot <- beta.matrix
beta.plot[is.na(beta.plot)] <- 0
beta.plot <- beta.plot %>%
  separate(unique_id, into = c("Rt", "Seg", "Yr_bin"), remove = F) %>%
  mutate(Site = paste0(Rt, "_", Seg))

ggplot(data = beta.plot, aes(x = Yr_bin, y = beta50.jac, group = Site, color = Site)) + 
  geom_point() + 
  geom_line() + 
  theme(legend.position = "none")

beta.mcmc <- beta.mcmc %>%
  separate(unique_id, into = c("route", "segment", "Year")) %>%
  mutate_if(is.character, as.numeric) %>%
  left_join(Bins, by = "Year") %>%
  mutate(unique_id = paste0(route, "_", segment, "_", Yr_bin)) %>%
  group_by(unique_id) %>%
  summarise(mean.beta = mean(mean.beta),
            beta.sd = mean(beta.sd))

beta.total <- left_join(beta.matrix, beta.mcmc, by = "unique_id")
beta.total <- beta.total %>% rename(beta.mcmc = mean.beta, beta.mcmc.sd = beta.sd)


#Remove years and duplicates in the bbs_total DF so we can link yr_bin calculations
#bbs_total$unique_id <- paste0(bbs_total$site, "_", bbs_total$Yr_bin)
#bbs_total <- bbs_total %>% select(-c("Year", "rt_yr", "Site_div"))
#bbs_total <- bbs_total[!duplicated(bbs_total), ]

#Add the betas into the bbs_total df 
#At this point the bbs_total DF has the necessary info to merge with 
#climate and LULC data for analyses at the binned level
bbs_total <- merge(bbs_total, beta.total, by = "unique_id")
bbs_total <- merge(bbs_total, site_alpha, by = "unique_id")

ggplot(data = bbs_total, aes(x = Yr_bin, y = beta.mcmc, group = site, color = site)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme(legend.position = "none")



#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, betaturn.slope = NA, betanest.slope = NA, betajac.slope = NA, 
                                beta50turn.slope = NA, beta50nest.slope = NA, beta50jac.slope = NA, 
                                betawetturn.slope = NA, betawetnest.slope = NA, betawetjac.slope = NA,
                                beta.mcmc.slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  bbs_total_temp <- bbs_total[bbs_total$site == site.temp,]
  bbs_total_temp <- bbs_total_temp[complete.cases(bbs_total_temp),]
  if (nrow(bbs_total_temp) > 1){
    lm.temp <- lm(bbs_total_temp$beta.turn ~ bbs_total_temp$Yr_bin)
    lm.temp2 <- lm(bbs_total_temp$beta.nest ~ bbs_total_temp$Yr_bin)
    lm.temp3 <- lm(bbs_total_temp$beta.jac ~ bbs_total_temp$Yr_bin)
    lm.temp4 <- lm(bbs_total_temp$beta50.turn ~ bbs_total_temp$Yr_bin)
    lm.temp5 <- lm(bbs_total_temp$beta50.nest ~ bbs_total_temp$Yr_bin)
    lm.temp6 <- lm(bbs_total_temp$beta50.jac ~ bbs_total_temp$Yr_bin)
    lm.temp7 <- lm(bbs_total_temp$beta.wet.turn ~ bbs_total_temp$Yr_bin)
    lm.temp8 <- lm(bbs_total_temp$beta.wet.nest ~ bbs_total_temp$Yr_bin)
    lm.temp9 <- lm(bbs_total_temp$beta.wet.jac ~ bbs_total_temp$Yr_bin)
    lm.temp10 <- lm(bbs_total_temp$beta.mcmc ~ bbs_total_temp$Yr_bin)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slope.temp2 <- summary(lm.temp2)$coefficients[2,1]
    slope.temp3 <- summary(lm.temp3)$coefficients[2,1]
    slope.temp4 <- summary(lm.temp4)$coefficients[2,1]
    slope.temp5 <- summary(lm.temp5)$coefficients[2,1]
    slope.temp6 <- summary(lm.temp6)$coefficients[2,1]
    slope.temp7 <- summary(lm.temp7)$coefficients[2,1]
    slope.temp8 <- summary(lm.temp8)$coefficients[2,1]
    slope.temp9 <- summary(lm.temp9)$coefficients[2,1]
    slope.temp10 <- summary(lm.temp10)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
    slopes_sites_beta[i,3] <- slope.temp2
    slopes_sites_beta[i,4] <- slope.temp3
    slopes_sites_beta[i,5] <- slope.temp4
    slopes_sites_beta[i,6] <- slope.temp5
    slopes_sites_beta[i,7] <- slope.temp6
    slopes_sites_beta[i,8] <- slope.temp7
    slopes_sites_beta[i,9] <- slope.temp8
    slopes_sites_beta[i,10] <- slope.temp9
    slopes_sites_beta[i,11] <- slope.temp10
  }
}

hist(slopes_sites_beta$betajac.slope)
t.test(slopes_sites_beta$betajac.slope, mu = 0)

hist(slopes_sites_beta$beta.mcmc.slope)
t.test(slopes_sites_beta$beta.mcmc.slope, mu = 0)

hist(slopes_sites_beta$beta50jac.slope)
t.test(slopes_sites_beta$beta50jac.slope, mu = 0)

hist(slopes_sites_beta$betawetjac.slope)
t.test(slopes_sites_beta$betawetjac.slope, mu = 0)


#Create a DF that calculates beta and alpha diversity in the year_bin formats
#for models 
bbs_bin <- bbs_total %>% dplyr::select(rteno.x, site, BCR, Site_div.x, rtename, Yr_bin, beta.turn, beta.nest, beta.jac, 
                                       beta.wet.turn, beta.wet.nest, beta.wet.jac, beta50.turn, beta50.nest, beta50.jac,
                                       beta.mcmc, SR_MCMC, SR50, SR.wet) %>%
  group_by(site, Yr_bin) %>% summarize(
    beta.turn = mean(beta.turn), beta.nest = mean(beta.nest), beta.jac = mean(beta.jac),
    beta.wet.turn = mean(beta.wet.turn), beta.wet.nest = mean(beta.wet.nest), beta.wet.jac = mean(beta.wet.jac),
    beta50.turn = mean(beta50.turn), beta50.nest = mean(beta50.nest), beta50.jac = mean(beta50.jac),
    beta.mcmc = mean(beta.mcmc), alpha = mean(Site_div.x), alpha.mcmc = mean(SR_MCMC), alpha50 = mean(SR50), alpha.wet = mean(SR.wet), 
  )


bbs_bin$unique_id <- paste0(bbs_bin$site, "_", bbs_bin$Yr_bin)

#Plots for Beta Diversity
bbs_total <- bbs_total[!bbs_total$BCR == 36, ]
bbs_total$BCR_name <- NA
bbs_total$BCR_name[bbs_total$BCR == "31"] <- "Peninsular Florida"
bbs_total$BCR_name[bbs_total$BCR == "37"] <- "Gulf Coast Prairie"
bbs_total$BCR_name[bbs_total$BCR == "27"] <- "Southeastern Coastal Plain"
bbs_total$BCR_name[bbs_total$BCR == "26"] <- "Mississippi Alluvial Valley"

gghist_beta <- ggplot(data = slopes_sites_beta, aes(slopes_sites_beta$betawetnest.slope)) + 
  #geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  geom_density(alpha = 0.4, fill = "#F98C0AFF") + #042333b2(Magma) #450A69FF(Viridis) #F98C0AFF (plasma)
  labs(title = "") +
  labs(x = expression(paste("Slopes of Wetland Bird ", beta, "-diversity"[italic("Nest")])), y = "# of Sites") +
  theme_cowplot(font_size = 15, line_size = 1.2) +
  coord_flip()

#site_data_merge$Year <- site_data_merge$count_yr + 1900

plot_beta <- (ggplot(bbs_total, aes(x = Yr_bin, y = beta.wet.turn, group = BCR, 
                                    colour = factor(BCR, 
                                                    labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = loess, se = T, aes(x = Yr_bin, y = beta.wet.nest, group = BCR, 
                                                        colour = factor(BCR,
                                                                        labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                #scale_color_brewer(palette = )
                xlab("Year") +
                ylab(expression(paste("Wetland Bird ", beta[italic("Nest")]))) +
                labs(colour = "Bird Conservation Region") +
                #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 15),
                      axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(vjust = -1, size = 15),
                      axis.title.y = element_text(vjust = 1.5, size = 15),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=15))) 

plot_beta <- plot_beta + scale_color_viridis(option = "plasma", discrete = T)

ggdraw() + 
  draw_plot(plot_beta + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.72, y = .030, width = .2, height = .75, scale = 1) 

#Species Richness Trends GoM
#Add a column to Summary_Stats to indicate 5 year bin to look at long term trends
#start at 1975 as this was the year standard protocols were implemented
# bbs_total <- bbs_total[bbs_total$Year >= 1980,]
# bbs_total$Yr_bin <- 1
# bbs_total$Yr_bin[bbs_total$Year >= 1985 & bbs_total$Year <= 1989] <- 2
# bbs_total$Yr_bin[bbs_total$Year >= 1990 & bbs_total$Year <= 1994] <- 3
# bbs_total$Yr_bin[bbs_total$Year >= 1995 & bbs_total$Year <= 1999] <- 4
# bbs_total$Yr_bin[bbs_total$Year >= 2000 & bbs_total$Year <= 2004] <- 5
# bbs_total$Yr_bin[bbs_total$Year >= 2005 & bbs_total$Year <= 2009] <- 6
# bbs_total$Yr_bin[bbs_total$Year >= 2010 & bbs_total$Year <= 2014] <- 7
# bbs_total$Yr_bin[bbs_total$Year >= 2015 & bbs_total$Year <= 2018] <- 8



#Mean_alpha_bin <- aggregate(data = bbs_total, Site_div ~ Yr_bin + site, FUN = "mean")
#Mean_alpha_bin <- separate(Mean_alpha_bin, site, c("state", "rteno", "yr"), sep = "_")

plot(bbs_total$alpha.bin ~ bbs_total$Yr_bin)
lm_alpha <- lm(alpha.bin ~ Yr_bin + site, data = bbs_total)
summary(lm_alpha)
anova(lm_alpha)
library(lme4)
lmer_alpha <- lmer(alpha.bin ~ Yr_bin + (1|site), data = bbs_total)
summary(lmer_alpha)
anova(lmer_alpha)
sjPlot::tab_model(lmer_alpha)
dev.off()

#Slopes for mean 1980 rarefied alphas 
n.sites <- length(unique(bbs_total$site))
site.list <- as.character(unique(bbs_total$site))

slopes_sites <- data.frame(site.list, slope = NA, slope50 = NA, slope.wet = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  bbs_temp <- bbs_total[bbs_total$site == site.temp & bbs_total$Year >= 1980, ]
  bbs_temp <- bbs_temp[!duplicated(bbs_temp$unique_id),]
  if (nrow(bbs_temp) > 1){
    lm.temp <- lm(bbs_temp$Site_div.x ~ bbs_temp$Year)
    lm.temp2 <- lm(bbs_temp$SR50 ~ bbs_temp$Year)
    lm.temp3 <- lm(bbs_temp$SR.wet ~ bbs_temp$Year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slope.temp2 <- summary(lm.temp2)$coefficients[2,1]
    slope.temp3 <- summary(lm.temp3)$coefficients[2,1]
    slopes_sites[i, 2] <- slope.temp
    slopes_sites[i, 3] <- slope.temp2
    slopes_sites[i, 4] <- slope.temp3
  }
}

plot(slopes_sites$slope ~ slopes_sites$site.list)
abline(h = 0)

t.test(slopes_sites$slope, mu = 0)

hist(slopes_sites$slope, breaks = 8, xlim = c(-1.5, 2.0), ylim = c(0, 50), main = "", 
     xlab = expression(paste("Slopes of ", alpha, "-diversity")), ylab = "# Sites")

gghist <- ggplot(data = slopes_sites, aes(slopes_sites$slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 30, binwidth = 0.25) +
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", alpha, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()


plot_alpha <- (ggplot(bbs_total, aes(x = Year, y = SR50, group = BCR_name, 
                                     colour = factor(BCR_name))) + 
                                                     #labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = loess, se = T, aes(x = Year, y = Site_div, group = BCR_name)) + 
                                                         #colour = factor(BCR, 
                                                        #                 labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 xlab("Year") +
                 ylab(expression(paste("Detection-correct ", alpha, "-diversity"))) +
                 labs(colour = "Bird Conservation Region") +
                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black", size =1.2),
                       axis.text.x = element_text(size = 14),
                       axis.text.y = element_text(size = 14),
                       axis.title.x = element_text(vjust = -1, size = 14),
                       axis.title.y = element_text(vjust = 1.5, size = 14),
                       axis.ticks = element_line(size = 1.2),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.margin = unit(c(1,1,2,2), "lines"),
                       text = element_text(size=14))) +
                scale_color_viridis(discrete = T)


plot_alpha
ggdraw() + 
  draw_plot(plot_alpha + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist, x = 0.75, y = .025, width = .2, height = .75, scale = 1) 
#draw_plot_label(c("A", "B"), c(0, 1.5), c(.75, 0.75), size = 12)


 
#######################################################################################################
###########################Climate at the Segment Level Centroid#######################################
#######################################################################################################

#Read in segment level climate data
prism.base <- read.csv(file = here::here("Data_Envi/PRISM Data/Segment_Level/PRISMsegs_Baseline.csv"), stringsAsFactors = F)
prism.seg1 <- read.csv(file = here::here("Data_Envi/PRISM Data/Segment_Level/PRISMsegs_1980_1994.csv"), stringsAsFactors = F)
prism.seg2 <- read.csv(file = here::here("Data_Envi/PRISM Data/Segment_Level/PRISMsegs_1995_2009.csv"), stringsAsFactors = F)
prism.seg3 <- read.csv(file = here::here("Data_Envi/PRISM Data/Segment_Level/PRISMsegs_2010_2019.csv"), stringsAsFactors = F)

#Bind data, clear NA's from spacing in csv, and separate dates
prism.seg <- rbind(prism.base, prism.seg1, prism.seg2, prism.seg3)
#prism.seg <- rbind(prism.seg1, prism.seg2, prism.seg3)
prism.seg <- prism.seg[complete.cases(prism.seg), ]
prism.seg <- separate(prism.seg, col = Date, into = c("Year", "Month"), by = "-", remove = T)

#Rename variables, wierd CSV notation
colnames(prism.seg)[colnames(prism.seg) == "ppt..inches."] <- "precip"
colnames(prism.seg)[colnames(prism.seg) == "tmax..degrees.F."] <- "tmax"
colnames(prism.seg)[colnames(prism.seg) == "tmin..degrees.F."] <- "tmin"
colnames(prism.seg)[colnames(prism.seg) == "tmean..degrees.F."] <- "tmean"

prism.seg <- prism.seg[!prism.seg$Year < 1975, ]

#Bin set up
Bins <- data.frame(matrix(ncol = 2, nrow = 45))
l <- c("Year", "Yr_bin")
colnames(Bins) <- l
Bins$Year <- 1975:2019
Bins$Yr_bin <- rep(1:9, times = 1, each = 5)

# Bins2 <- data.frame(matrix(ncol = 2, nrow = 10))
# l2 <- c("Year", "Yr_bin")
# colnames(Bins2) <- l2
# Bins2$Year <- 1970:1979
# Bins2$Yr_bin <- 1

#Bins <- rbind(Bins2, Bins)

#Pull the yr_bins into the data
prism.seg <- merge(prism.seg, Bins, by = "Year")

#Bin months for season (use only spring and summer)
prism.seg$Month <- as.numeric(prism.seg$Month)
prism.seg$season <- "Winter"
prism.seg$season[prism.seg$Month >= 04 & prism.seg$Month <= 06] <- "Spring"
prism.seg$season[prism.seg$Month >= 07 & prism.seg$Month <= 09] <- "Summer"
prism.seg$season[prism.seg$Month >= 10 & prism.seg$Month <= 12] <- "Fall"
prism.seg$season[prism.seg$Month >= 01 & prism.seg$Month <= 03] <- "Winter"
prism.seg$season[prism.seg$Month >= 05 & prism.seg$Month <= 08] <- "Breeding"
prism.seg$RainSeason <- "Dry"
prism.seg$RainSeason[prism.seg$Month >= 11 & prism.seg$Month <= 04] <- "Dry"
prism.seg$RainSeason[prism.seg$Month >= 05 & prism.seg$Month <= 10] <- "Wet"
#prism.seg$RainSeason[prism.seg$Month >= 05 & prism.seg$Month <= 08] <- "Breeding"



#Reduce data for just Spring and Summer seasons
prism.seg.sp <- prism.seg[prism.seg$season == "Spring", ] 
prism.seg.s <- prism.seg[prism.seg$season == "Summer", ]
prism.seg.f <- prism.seg[prism.seg$season == "Fall", ]
prism.seg.w <- prism.seg[prism.seg$season == "Winter", ]
prism.seg.bt <- prism.seg[prism.seg$season == "Breeding", ]
prism.seg.dry <- prism.seg[prism.seg$RainSeason == "Dry", ]
prism.seg.wet <- prism.seg[prism.seg$RainSeason == "Wet", ]
#prism.seg.bp <- prism.seg[prism.seg$RainSeason == "Breeding", ]


#Put indicator of season in colnames for calling 
colnames(prism.seg.sp)[7:10] <- paste0(colnames(prism.seg.sp)[7:10], "_", unique(prism.seg.sp$season)) 
colnames(prism.seg.s)[7:10] <- paste0(colnames(prism.seg.s)[7:10], "_", unique(prism.seg.s$season)) 
colnames(prism.seg.f)[7:10] <- paste0(colnames(prism.seg.f)[7:10], "_", unique(prism.seg.f$season)) 
colnames(prism.seg.w)[7:10] <- paste0(colnames(prism.seg.w)[7:10], "_", unique(prism.seg.w$season)) 
colnames(prism.seg.bt)[7:10] <- paste0(colnames(prism.seg.bt)[7:10], "_", unique(prism.seg.bt$season)) 
colnames(prism.seg.dry)[7:10] <- paste0(colnames(prism.seg.dry)[7:10], "_", unique(prism.seg.dry$RainSeason)) 
colnames(prism.seg.wet)[7:10] <- paste0(colnames(prism.seg.wet)[7:10], "_", unique(prism.seg.wet$RainSeason)) 
#colnames(prism.seg.bp)[7:10] <- paste0(colnames(prism.seg.bp)[7:10], "_", unique(prism.seg.bp$RainSeason)) 

#Calculate Yr_bin & Seasonal Means for site 
prism.seg.sp <- prism.seg.sp %>% group_by(Site, Year) %>% summarise(tmean_Spring = mean(tmean_Spring),
                                                                    tmin_Spring = mean(tmin_Spring), 
                                                                    tmax_Spring = mean(tmax_Spring),
                                                                    precip_Spring = mean(precip_Spring)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()

prism.seg.s <- prism.seg.s %>% group_by(Site, Year) %>% summarise(tmean_Summer = mean(tmean_Summer),
                                                                  tmin_Summer = mean(tmin_Summer), 
                                                                  tmax_Summer = mean(tmax_Summer),
                                                                  precip_Summer = mean(precip_Summer)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()

prism.seg.f <- prism.seg.f %>% group_by(Site, Year) %>% summarise(tmean_Fall = mean(tmean_Fall),
                                                                  tmin_Fall = mean(tmin_Fall), 
                                                                  tmax_Fall = mean(tmax_Fall),
                                                                  precip_Fall = mean(precip_Fall)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()

prism.seg.w <- prism.seg.w %>% group_by(Site, Year) %>% summarise(tmean_Winter = mean(tmean_Winter),
                                                                  tmin_Winter = mean(tmin_Winter), 
                                                                  tmax_Winter = mean(tmax_Winter),
                                                                  precip_Winter = mean(precip_Winter)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()

prism.seg.bt <- prism.seg.bt %>% group_by(Site, Year) %>% summarise(tmean_Breeding = mean(tmean_Breeding),
                                                                  tmin_Breeding = mean(tmin_Breeding), 
                                                                  tmax_Breeding = mean(tmax_Breeding),
                                                                  precip_Breeding = mean(precip_Breeding)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()


prism.seg.wet <- prism.seg.wet %>% group_by(Site, Year) %>% summarise(tmean_Wet = mean(tmean_Wet),
                                                                      tmin_Wet = mean(tmin_Wet),
                                                                      tmax_Wet = mean(tmax_Wet),
                                                                      precip_Wet = mean(precip_Wet)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()

prism.seg.dry <- prism.seg.dry %>% group_by(Site, Year) %>% summarise(tmean_Dry = mean(tmean_Dry),
                                                                      tmin_Dry = mean(tmin_Dry),
                                                                      tmax_Dry = mean(tmax_Dry),
                                                                      precip_Dry = mean(precip_Dry)) %>%
  unite("unique_id", c("Site", "Year")) %>% ungroup()


#Bring all the seasons together 
prism.seg.means <- prism.seg.sp %>% 
  left_join(prism.seg.s, by = "unique_id") %>% 
  left_join(prism.seg.f, by = "unique_id") %>%
  left_join(prism.seg.w, by = "unique_id") %>% 
  left_join(prism.seg.bt, by = "unique_id") %>%
  left_join(prism.seg.dry, by = "unique_id") %>%
  left_join(prism.seg.wet, by = "unique_id")

prism.seg.means <- prism.seg.means %>% separate(col = unique_id, into = c("Rteno", "Segment", "Year"))

prism.seg.means$Site <- paste0(prism.seg.means$Rteno, "_", prism.seg.means$Segment) 

#Calculate the total mean column for each site 
prism.seg.means <- prism.seg.means %>% mutate(tmean = (tmean_Spring + tmean_Summer + tmean_Fall + tmean_Winter) / 4,
                                              tmax = (tmax_Spring + tmax_Summer + tmax_Fall + tmax_Winter) / 4,
                                              tmin = (tmin_Spring + tmin_Summer + tmin_Fall + tmin_Winter) / 4,
                                              pmean = (precip_Spring + precip_Summer + precip_Fall + precip_Winter) / 4)
                                              #tmean.bird = (tmean_Spring + tmean_Summer) / 2,
                                              #tmin.bird = (tmin_Spring + tmin_Summer) / 2,
                                              #tmax.bird = (tmax_Spring + tmax_Summer) / 2,
                                              #pmean.bird = (precip_Spring + precip_Summer) / 2)

#Metric Conversion
seg.climate <- prism.seg.means %>% mutate(tmean.c = (tmean - 32) / 1.8, 
                                          tmax.c = (tmax - 32) / 1.8,
                                          tmin.c = (tmin - 32) / 1.8,
                                          tmean.breed.c = (tmean_Breeding - 32) / 1.8,
                                          tmin.breed.c = (tmin_Breeding - 32) / 1.8,
                                          tmax.breed.c = (tmax_Breeding - 32) / 1.8,
                                          pmean.c = (pmean * 2.54),
                                          tmean.spring.c = (tmean_Spring - 32) / 1.8,
                                          tmean.summer.c = (tmean_Summer - 32) / 1.8,
                                          tmean.fall.c = (tmean_Fall - 32) / 1.8,
                                          tmean.winter.c = (tmean_Winter - 32) / 1.8,
                                          tmin.spring.c = (tmin_Spring - 32) / 1.8,
                                          tmin.summer.c = (tmin_Summer - 32) / 1.8,
                                          tmin.fall.c = (tmin_Fall - 32) / 1.8,
                                          tmin.winter.c = (tmin_Winter - 32) / 1.8,
                                          tmax.spring.c = (tmax_Spring - 32) / 1.8,
                                          tmax.summer.c = (tmax_Summer - 32) / 1.8,
                                          tmax.fall.c = (tmax_Fall - 32) / 1.8,
                                          tmax.winter.c = (tmax_Winter - 32) / 1.8,
                                          tmean.wet.c = (tmean_Wet - 32) / 1.8,
                                          tmax.wet.c = (tmax_Wet - 32) / 1.8,
                                          tmin.wet.c = (tmin_Wet - 32) / 1.8,
                                          tmean.dry.c = (tmean_Dry - 32) / 1.8,
                                          tmax.dry.c = (tmax_Dry - 32) / 1.8,
                                          tmin.dry.c = (tmin_Dry - 32) / 1.8,
                                          precip.spring.c = (precip_Spring * 2.54),
                                          precip.summer.c = (precip_Summer * 2.54),
                                          precip.fall.c = (precip_Fall * 2.54),
                                          precip.winter.c = (precip_Winter * 2.54),
                                          precip.wet.c = (precip_Wet * 2.54),
                                          precip.dry.c = (precip_Dry * 2.54),
                                          precip.breed.c = (precip_Breeding * 2.54),)

#Old command when I wanted to keep all of the climate covariates
# seg.climate <- seg.climate[, c("Site", "Rteno", "Year", "tmean.c", "tmax.c", "tmin.c", "pmean.c",
#                                "tmean.bird.c", "tmax.bird.c", "tmin.bird.c", "pmean.bird.c", "tmean_Spring.c",
#                                "tmin_Spring.c", "tmax_Spring.c", "precip_Spring.c", "tmean_Summer.c", "tmax_Summer.c",
#                                "tmin_Summer.c", "precip_Summer.c", "tmean_Fall.c", "tmin_Fall.c", "tmax_Fall.c", "precip_Fall.c",
#                                "tmean_Winter.c", "tmin_Winter.c", "tmax_Winter.c", "precip_Winter.c", "precip_Dry.c", "precip_Wet.c",
#                                "tmax_Wet.c", "tmin_Wet.c", "tmax_Dry.c", "tmin_Dry.c", "tmean_Dry.c", "tmean_Wet.c")]

seg.climate <- seg.climate[, c("Site", "Rteno", "Year", 
                               "tmean.breed.c", "tmax.breed.c", 
                               "tmin.breed.c", "precip.breed.c")]
                               


# lag.calc <- function(x) {
#   x - lag(x)
# }

diff.fun <- function(x) {
  x - first(x)
}

#clim.vars <- names(seg.climate)[4:35]
clim.vars <- names(seg.climate)[4:7]

bbs.info <- beta.matrix %>%
  dplyr::select(Site, Year) %>%
  mutate(unique_id = paste0(Site, "_", Year)) %>%
  dplyr::select(unique_id) %>%
  mutate(as.character(unique_id))

bbs.info <- bbs.info$unique_id  

survey.range <- beta.matrix %>%
  dplyr::select(Site, Year, FirstYr, LastYr) %>%
  #mutate(unique_id = paste0(Site, "_", Year)) %>%
  ungroup() %>%
  dplyr::select(Site, FirstYr, LastYr) %>%
  distinct()

seg.clim <- seg.climate %>%
  mutate(unique_id = paste0(Site, "_", Year)) %>%
  filter(unique_id %in% bbs.info) 

seg.climate <- seg.climate %>%
  right_join(survey.range, by = "Site") %>%
  mutate_at(c("FirstYr", "LastYr", "Year"), ~as.numeric(.))
  

first.yrs <- seg.climate %>% group_by(Site) %>% arrange(Year) %>%
          mutate(Year = as.numeric(Year),
                 index.yr.min = FirstYr - Year,
                 index.yr.max = LastYr - Year,
                 Year = as.character(Year)) %>%
          filter(between(index.yr.min, 0, 4)) %>%
          mutate_at(c("FirstYr", "LastYr"), ~as.character(.)) %>%
          dplyr::select(-index.yr.min, -index.yr.max) %>%
          summarize_if(is.numeric, mean) %>% ungroup()

first.yrs <- first.yrs %>%
  left_join(survey.range) %>%
  mutate(Year = FirstYr) %>% 
  separate(Site, into = c("Rteno", "Segment"), remove = F) %>%
  dplyr::select(-Segment) %>%
  mutate(unique_id = paste0(Site, "_", Year))
  


seg.anom <- seg.climate %>%
  mutate(Year = as.character(Year),
         FirstYr = as.character(FirstYr),
         LastYr = as.character(LastYr)) %>%
  mutate(unique_id = paste0(Site, "_", Year)) %>%
  filter(unique_id %in% bbs.info) %>%
  group_by(Site) %>%
  arrange(Site, Year) %>%
  slice(-1) %>%
  ungroup() %>%
  bind_rows(first.yrs) %>%
  group_by(Site) %>%
  arrange(Site, Year) %>%
  mutate_at(vars(-Site, -Rteno, -Year, -unique_id, -FirstYr, -LastYr), function(x){x - first(x)}) %>%
  rename_at(vars(clim.vars), ~paste0(.,"_anom")) %>%
  ungroup() %>%
  dplyr::select(-Rteno, -Year, -Site) 

seg.clim.fin <- seg.clim %>%
  left_join(seg.anom, by = "unique_id") %>%
  mutate(unique_id = paste0(Site, "_", Year))


Bins <- data.frame(matrix(ncol = 2, nrow = 40))
l <- c("Year", "Yr_bin")
colnames(Bins) <- l
Bins$Year <- 1980:2019
Bins$Yr_bin <- rep(1:8, times = 1, each = 5)
Bins$Year <- as.character(Bins$Year)

bbs.clim <- beta.matrix %>%
  mutate(unique_id = paste0(Site, "_", Year)) %>%
  right_join(seg.clim.fin) %>%
  right_join(Bins, by = "Year")

unique_id <- paste0(bbs.clim$Site, "_", bbs.clim$Yr_bin)

#Put betas with temperature data 
#DF contains unique_id, site, yr_bin, raw climate, and anomalies
#lmer.df <- merge(temp.avgs, bbs_bin, by = "unique_id")


#Bring in the LULC Data 
bbs_lulc <- read.csv(here::here("Data_BBS/Generated DFs/BBS_LULC.csv"))
colnames(bbs_lulc)[colnames(bbs_lulc) == "year"] <- "Year"

bbs_ma <- read.csv(here::here("Data_Envi/BBS_Mangrove_50stop/MangroveData.csv"), stringsAsFactors = F)

bbs_lulc <- bbs_lulc %>% dplyr::select(-c(site, X, unique_id))

#Fix years for LULC
bbs_lulc[bbs_lulc$Year == 1986, 2] <- 1985
bbs_lulc[bbs_lulc$Year == 1996, 2] <- 1995
bbs_lulc[bbs_lulc$Year == 2001, 2] <- 2000
bbs_lulc[bbs_lulc$Year == 2006, 2] <- 2005
bbs_lulc[bbs_lulc$Year == 2016, 2] <- 2015

bbs_lulc <- merge(Bins, bbs_lulc, by = "Year")
bbs_lulc$site <- paste0(bbs_lulc$Route, "_", bbs_lulc$Segment)
bbs_lulc$unique_id <- paste0(bbs_lulc$rteno.x, "_", bbs_lulc$Yr_bin)


 
mangrove.list <- unique(bbs_ma$site)
lulc.list <- unique(bbs_lulc$site)
lulc.list <- as.character(lulc.list)

#Mangrove Data
bbs_ma <- merge(bbs_ma, Bins, by.x = "year", by.y = "Year")

#Make link df for mangrove
seg.link <- wetland.link[, c("site", "rteno.segment")]
seg.link$site <- as.character(seg.link$site)
seg.link <- seg.link[!duplicated(seg.link$rteno.segment), ]
seg.sites <- rep(unique(seg.link$rteno.segment), 1, each = 8)
seg.year <- rep(seq(from = 1980, to = 2015, by = 5), time = 380, each = 1)
seg.rts <- rep(unique(seg.link$site), 1, each = 8)

seg.link <- as.data.frame(cbind(seg.rts, seg.sites, seg.year))
seg.link$unique_id <- paste0(seg.link$seg.rts, "_", seg.link$seg.year)

#Link up mangroves wth rt#
bbs_ma <- merge(bbs_ma, seg.link, by = "unique_id")

#Pul out the values that exist in the BBS_lulc DF
bbs_ma <- bbs_ma[bbs_ma$site %in% lulc.list, ]
#bbs_ma <- merge(bbs_ma, Bins, by.x = "year", by.y = "Year")
bbs_ma$unique_id <- paste0(bbs_ma$seg.sites, "_", bbs_ma$seg.year)

bbs_lulc$Emergent_Wetlands <- bbs_lulc$Emergent_Wetlands + 0.000001
bbs_lulc$Woody_Wetlands <- bbs_lulc$Woody_Wetlands + 0.000001


#bbs_lulc$mangrove <- NA
bbs_ma <- bbs_ma[bbs_ma$mangrove > 0, ]
bbs_ma$unique_id <- paste0(bbs_ma$seg.sites, "_", bbs_ma$Yr_bin)
#Write a loop to place values in the mangrove DF into the BBS_LULC DF
# for (o in 1:length(bbs_ma$unique_id)){
#   mangrove.tmp <- bbs_ma[o, ]
#   mangrove.val <- mangrove.tmp[5]
#   mangrove.val <- mangrove.val[1,1]
#   bbs_lulc[bbs_lulc$unique_id == mangrove.tmp$unique_id, "mangrove"] <- mangrove.val
# }

bbs_lulc <- bbs_ma %>% dplyr::select(c("unique_id", "mangrove")) %>% full_join(bbs_lulc, by = "unique_id") 

bbs_lulc[is.na(bbs_lulc$mangrove), 2] <- 0
#bbs_lulc <- bbs_lulc %>% group_by(site) %>% mutate(pct.man = mangrove / total_cover_nb, diff.from.first.man = (pct.man - first(pct.man)))#,
#scale.pman = scale(pct.man), scale.pdman = scale(diff.from.first.man))

#Merge with the site centroid DF
#bbs_full <- bbs_lulc %>% left_join(temp.avgs, by = "unique_id") %>% ungroup()
#bbs_lulc <- merge(Bins, bbs_lulc, by.x = "Year", by.y = "year")
bbs_lulc$unique_id <- paste0(bbs_lulc$rteno.x, "_", bbs_lulc$Yr_bin)

bbs_lulc <- bbs_lulc %>%
  filter(bbs_lulc$unique_id %in% unique_id) %>%
  group_by(rteno.x) %>%
  arrange(rteno.x) %>%
  mutate(total_cover_nb = Urban + Ag + Grassland + Forest + Woody_Wetlands + Emergent_Wetlands + Bare + Water,
         pct.nat = (Grassland + Forest) / total_cover_nb, pct.human = (Urban + Ag) / total_cover_nb,
         pct.wetland = (Woody_Wetlands + Emergent_Wetlands) / total_cover_nb, pct.ag = Ag / total_cover_nb,
         pct.grass = Grassland / total_cover_nb, pct.ww = Woody_Wetlands / total_cover_nb, pct.ew = Emergent_Wetlands / total_cover_nb, 
         pct.ur = Urban / total_cover_nb, pct.for = Forest / total_cover_nb, pct.wat = Water / total_cover_nb, pct.bare = Bare / total_cover_nb, 
         pct.man = (mangrove / total_cover_nb),WW_NoMan = Woody_Wetlands - mangrove, pct.wwnm = WW_NoMan / total_cover_nb, ratio.ww = (Woody_Wetlands / Emergent_Wetlands),
         ratio.ew = (Emergent_Wetlands / Woody_Wetlands), ratio.man = (mangrove / Emergent_Wetlands), ratio.mnww = (mangrove / WW_NoMan))%>%
  mutate(diff.ww = (pct.ww - first(pct.ww))) %>%
  mutate(ratio.ww = (Woody_Wetlands / Emergent_Wetlands)) %>%
  mutate(ratio.man = (mangrove / Emergent_Wetlands)) %>%
  mutate(ratio.manww = (mangrove / WW_NoMan)) %>%
  mutate(diff.ew = (pct.ew - first(pct.ew))) %>%
  mutate(ratio.ew = (Emergent_Wetlands / Woody_Wetlands)) %>%
  mutate(diff.ratio.wet = (ratio.ww - first(ratio.ww))) %>%
  mutate(diff.man = (pct.man - first(pct.man))) %>%
  mutate(diff.wwnm = (pct.wwnm - first(pct.wwnm))) %>%
  mutate(diff.nat = (pct.nat - first(pct.nat))) %>%
  mutate(diff.human = (pct.human - first(pct.human))) %>%
  mutate(diff.ur = (pct.ur - first(pct.ur))) %>%
  mutate(diff.ag = (pct.ag - first(pct.ag))) %>%
  mutate(diff.wet = (pct.wetland - first(pct.wetland))) %>%
  mutate(diff.for = (pct.for - first(pct.for))) %>%
  mutate(diff.bare = (pct.bare - first(pct.bare))) %>%
  mutate(diff.wat = (pct.wat - first(pct.wat))) %>%
  mutate(diff.grass = (pct.grass - first(pct.grass))) %>%
  rename(Site = rteno.x) %>%
  ungroup() %>%
  dplyr::select(-rteno, -site, -Year, -Segment, -rteno, -unique_id, -Route)

Bins$Year <- as.character(Bins$Year)

#Take the remaining missing info for analyses
bbs.info <- bbs_total %>%
  dplyr::select(Site = site, BCR, BCR_name, Latitude, Longitude) %>%
  distinct()

#Merge it all together
bbs.comp <- bbs.clim %>%
  inner_join(bbs_lulc, by = c("Site", "Yr_bin")) %>%
  arrange(Site, Year) %>%
  left_join(bbs.info, by = "Site") %>%
  dplyr::select(Site, Year, Yr_bin, Route, BCR, BCR_name, Latitude, Longitude,
                beta = beta.wet.jac, turn = beta.wet.turn, nest = beta.wet.nest, alpha = SR.wet,
                beta.full = beta50.jac, turn.full = beta50.turn, nest.full = beta50.nest, alpha.full = SR50,
                tmean = tmean.breed.c, tmax = tmax.breed.c, tmin = tmin.breed.c, ppt = precip.breed.c, 
                tmean.anom = tmean.breed.c_anom, tmin.anom = tmin.breed.c_anom, tmax.anom = tmax.breed.c_anom,
                ppt.anom = precip.breed.c_anom, 
                pct.wetland, pct.nat, pct.human,
                pct.man, pct.ww, pct.wwnm, pct.ew,
                pct.ag, pct.ur, pct.for, pct.wat, pct.bare, pct.grass,
                diff.wet, diff.man, diff.ww, diff.wwnm, diff.ew, diff.nat, 
                diff.human, diff.ur, diff.ag, diff.for, diff.wat, diff.grass,
                diff.bare, ratio.ww, ratio.ew, ratio.man, ratio.mnww, ratio.man)


#write.csv(bbs.comp, here::here("Data_BBS/Generated DFs/INLAdfFinalRevised.csv"))


#############################################################################################
#############################CMRL for the Occupancy Data#####################################
#############################################################################################

#Load in the occupancy DF
occ50 <- read.csv(file = here::here("Data_BBS/Generated DFs/occ50.csv"), stringsAsFactors = F)
occ50 <- occ50[, -1:-2]

#Load in the CMRL Max Ranges
spp.range <- read.csv(here::here("Data_BBS/Generated DFs/MaxRange.csv"), stringsAsFactors = T)

#Clean some names up for matching
spp.range$sci_name <- as.character(spp.range$sci_name)
spp.range$sci_name[spp.range$sci_name == "Antigone canadensis"] <- "Grus canadensis"
spp.range$sci_name[spp.range$sci_name == "Hydroprogne caspia"] <- "Sterna caspia"
spp.range$sci_name[spp.range$sci_name == "Spinus psaltria"] <- "Carduelis psaltria"
spp.range$sci_name[spp.range$sci_name == "Colaptes auratus auratus"] <- "Colaptes auratus"

#Merge occ data with CMRL
occ50 <- merge(occ50, spp.range, by.x = "Species", by.y = "sci_name")
present <- occ50[occ50$Occupancy == 1,]

#Calculate the CMRL by Year
cmrl.occ <- aggregate(data = present, Latitude ~ site + Year, FUN = mean)
lm.cmrl <- lmer(data = cmrl.occ, Latitude ~ Year + (1|site), REML = F)
summary(lm.cmrl)
Anova(lm.cmrl)
tab_model(lm.cmrl)

cmrl.occ <- merge(cmrl.occ, Bins, by = "Year")

#Get the change in CMRL 
cmrl.occ <- cmrl.occ %>% group_by(site) %>% arrange(Year) %>% #summarise(CMRL = mean(Latitude)) %>% 
  mutate(change.cmrl = Latitude - first(Latitude),
         change.cmrl.km = change.cmrl * 111)
#Find the mean by site
cmrl.mean <- cmrl.occ %>% group_by(site) %>% summarise(mean.change = mean(change.cmrl.km))

#Plot some results
ggplot(data = cmrl.occ, aes(x = Year, y = Latitude, group = site)) + 
  geom_smooth(data = cmrl.occ, aes(x = Year, y = Latitude, group = site), method = lm, se = F)

#Pull out the last for each group
cmrl.last <- cmrl.occ %>% rename(CMRL = Latitude, Site = site) %>% arrange(Yr_bin) %>% slice(n()) %>% dplyr::select(-c("Year", "Yr_bin"))

#Put the cmrl with the last mod DF
bbs_last <- merge(bbs_last, cmrl.last, by = "Site")

#write.csv(bbs_last, here::here("Data_BBS/Generated DFs/DF4Analysis_Seg_Last_New.csv"))

#bbs_last <- bbs_last %>% mutate(scale.cmrl = scale(change.cmrl.km))

cmrl.occ$unique_id <- paste0(cmrl.occ$site, "_", cmrl.occ$Yr_bin)
bbs_full <- merge(bbs_full, cmrl.occ, by = "unique_id")

#save.image(here::here("R Workspace/Diversity_Script_Full.RData"))
