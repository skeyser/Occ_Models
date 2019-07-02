#######################################################################
#############Script for calculating diversity metrics##################
#######and pairing avian data with the LULC and climate data###########
#######################for model running###############################
##################Created by: Spencer R Keyser#########################
#######################################################################


##Clean Global Environment  
rm(list = ls())

#Package Loading
#install.packages("pacman")
library("pacman")
pacman::p_load("tidyverse", "reshape2", "vegan", "data.table", "cowplot")
# pacman::p_load("sjPlot", "sjstats", "nlme","reshape2", "effects","vegan", 
#                "lmerTest", "tidyverse", "lsr", "cowplot", "here", "purrr", 
#                "lme4", "car", "MASS", "MuMIn")

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


#Clean the workspace except for the files needed
#rm(list = setdiff(ls(), c()))



#Pull Out the Lat/Long
site.list <- unique(final_sp_df$rteno.x)
rt.xy <- rts_gom[rts_gom$rteno %in% site.list,]
rt.xy <- rt.xy[rt.xy$BCR != 19, ]
rt.xy <- rt.xy %>% select(rtename, rteno, countrynum, Latitude,
                          Longitude)

#Create a spp by site matrix
#Keep only the columns we need
#Create site x year id
final_sp_df$unique_id <- paste0(final_sp_df$site, "_", final_sp_df$Year)

bbs_simple <- final_sp_df[,c("sci_name", "unique_id")]
bbs_simple$Detected <- 1
#bbs_plot <- bbs_total[, c("unique_id", "statenum")]
#bbs_plot <- as.data.frame(bbs_plot)

#Cast Data
bbs_cast <- dcast(data = bbs_simple, formula = unique_id ~ sci_name, fun.aggregate = sum)

rownames(bbs_cast) <- bbs_cast$unique_id
bbs_cast <- bbs_cast[,-1]

Site_div <- rowSums(bbs_cast)
Site_div <- as.data.frame(Site_div)
Site_div$unique_id <- rownames(Site_div)
rownames(Site_div) <- NULL


#check names are the same so we can just cbind them
bbs_total <- final_sp_df[!duplicated(final_sp_df),]
bbs_total <- bbs_total %>% dplyr::select(rteno.x, rt_yr, Segment, Year,
                                  BCR, site, unique_id)
bbs_total <- bbs_total[!duplicated(bbs_total),]
bbs_total <- arrange(bbs_total, desc(unique_id))
Site_div <- arrange(Site_div, desc(unique_id))
bbs_total <- cbind(bbs_total, Site_div)
bbs_total <- bbs_total[, -9]

bbs_total <- merge(bbs_total, rt.xy, by.x = "rteno.x", by.y = "rteno")


#Mean SR for each route
bbs_total$site <- as.character(bbs_total$site)
Site.means <- bbs_total %>% group_by(site) %>%
              summarise(mean.alpha = mean(Site_div))


n.sites <- length(unique(bbs_total$site))
site.list <- as.character(unique(bbs_total$site))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_total_temp <- bbs_total[bbs_total$site == site.temp,]
  if (nrow(site_total_temp) > 1){
    lm.temp <- lm(site_total_temp$Site_div ~ site_total_temp$Year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites[i,2] <- slope.temp
  }
}

hist(slopes_sites$slope, breaks = 50)

slopes_sites$site <- as.character(slopes_sites$site.list)
slopes_sites <- slopes_sites %>% dplyr::select(site, slope)

site_means2 <- merge(Site.means, slopes_sites, by = "site")

mean(site_means2$slope, na.rm = T)
median(site_means2$slope, na.rm = T)


#mean temporal beta diversity by site
mean.betas <- data.frame(site.list, mean.beta = NA)

#Create matrix for storing Betas 
beta.matrix <- data.frame(unique_id = bbs_total$unique_id, beta = NA)
beta.matrix$unique_id <- as.character(beta.matrix$unique_id)

#Loop for calculating Beta Diversity from baseline 5 years
#If I average the Beta diversity after will this generate 
#statistical biases???
#Should we loop it from year bin to year bin?
for (n in 1:n.sites){
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(bbs_total[which(bbs_total$site == site.temp),"Year"])
  yrs.temp <- yrs.temp[yrs.temp >= 1980]
  #only consider sites that were observed in more than 5 years 
  if (length(yrs.temp) > 5){
    #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years 
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1] + 4)]
    print(paste0("Number of years for site ", site.temp, " is ", length(first.yrs)))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years 
    bbs_cast_first <- bbs_cast[rownames(bbs_cast) %in% site.yrs.first,]
    first.comm.temp <- colMeans(bbs_cast_first)
    #get remaining years 
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1] + 4)]
    #Now loop the remaining years 
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y]
      other.site.temp <- paste0(site.temp, "_", other.yr.temp)
      bbs_cast_other <- bbs_cast[rownames(bbs_cast) == other.site.temp,]
      bbs_cast_other <- rbind(bbs_cast_other, first.comm.temp)
      beta.temp <- vegdist(sqrt(sqrt(bbs_cast_other)), method = "jaccard")
      beta.matrix$beta[beta.matrix$unique_id == other.site.temp] <- beta.temp
    }
  }}


#Put betas back in site daTA MERGE
bbs_total <- merge(bbs_total, beta.matrix, by = "unique_id")

#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, beta.slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  bbs_total_temp <- bbs_total[bbs_total$site == site.temp,]
  bbs_total_temp <- bbs_total_temp[complete.cases(bbs_total_temp),]
  if (nrow(bbs_total_temp) > 1){
    lm.temp <- lm(bbs_total_temp$beta ~ bbs_total_temp$Year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
  }
}

hist(slopes_sites_beta$beta.slope)
t.test(slopes_sites_beta$beta.slope, mu = 0)

#####Cleaned to this point#####
###############################



#Plots for Beta Diversity

gghist_beta <- ggplot(data = slopes_sites_beta, aes(slopes_sites_beta$beta.slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", beta, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()

#site_data_merge$Year <- site_data_merge$count_yr + 1900


plot_beta <- (ggplot(bbs_total, aes(x = Year, y = beta, group = site, 
                                          colour = factor(BCR, 
                                                          labels = c("BCR 26", "BCR 27", "BCR 31", "BCR 36", "BCR 37")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = lm, se = FALSE, aes(x = Year, y = beta, group = site, 
                                                         colour = factor(BCR,
                                                                         labels = c("BCR 26", "BCR 27", "BCR 31", "BCR 36", "BCR 37")))) +
                xlab("Year") +
                ylab(expression(paste(beta, "-diversity"))) +
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
                      text = element_text(size=14)))

ggdraw() + 
  draw_plot(plot_beta + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.75, y = .025, width = .2, height = .75, scale = 1) 

#Species Richness Trends GoM
#Add a column to Summary_Stats to indicate 5 year bin to look at long term trends
#start at 1975 as this was the year standard protocols were implemented
bbs_total <- bbs_total[bbs_total$Year >= 1980,]
bbs_total$Yr_bin <- 1
bbs_total$Yr_bin[bbs_total$Year >= 1985 & bbs_total$Year <= 1989] <- 2
bbs_total$Yr_bin[bbs_total$Year >= 1990 & bbs_total$Year <= 1994] <- 3
bbs_total$Yr_bin[bbs_total$Year >= 1995 & bbs_total$Year <= 1999] <- 4
bbs_total$Yr_bin[bbs_total$Year >= 2000 & bbs_total$Year <= 2004] <- 5
bbs_total$Yr_bin[bbs_total$Year >= 2005 & bbs_total$Year <= 2009] <- 6
bbs_total$Yr_bin[bbs_total$Year >= 2010 & bbs_total$Year <= 2014] <- 7
bbs_total$Yr_bin[bbs_total$Year >= 2015 & bbs_total$Year <= 2018] <- 8



Mean_alpha_bin <- aggregate(data = bbs_total, Site_div ~ Yr_bin + site, FUN = "mean")
#Mean_alpha_bin <- separate(Mean_alpha_bin, site, c("state", "rteno", "yr"), sep = "_")

plot(Mean_alpha_bin$Site_div ~ Mean_alpha_bin$Yr_bin)
lm_alpha <- lm(Site_div ~ Yr_bin + site, data = Mean_alpha_bin)
summary(lm_alpha)
anova(lm_alpha)
library(lme4)
lmer_alpha <- lmer(Site_div ~ Yr_bin + (1|site), data = Mean_alpha_bin)
summary(lmer_alpha)
anova(lmer_alpha)
sjPlot::tab_model(lmer_alpha)
dev.off()

#Slopes for mean 1980 rarefied alphas 
n.sites <- length(unique(bbs_total$site))
site.list <- as.character(unique(bbs_total$site))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  bbs_temp <- bbs_total[bbs_total$site == site.temp & bbs_total$Year >= 1980, ]
  bbs_temp <- bbs_temp[!duplicated(bbs_temp$unique_id),]
  if (nrow(bbs_temp) > 1){
    lm.temp <- lm(bbs_temp$Site_div ~ bbs_temp$Year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites[i,2] <- slope.temp
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


plot_alpha <- (ggplot(bbs_total, aes(x = Year, y = Site_div, group = site, 
                                           colour = factor(BCR, 
                                                           labels = c("BCR 26", "BCR 27", "BCR 31", "BCR 36", "BCR 37")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = glm,se = FALSE, aes(x = Year, y = Site_div, group = site, 
                                                          colour = factor(BCR, 
                                                                          labels = c("BCR 26", "BCR 27", "BCR 31", "BCR 36", "BCR 37")))) +
                 xlab("Year") +
                 ylab(expression(paste("Rarefied ", alpha, "-diversity"))) +
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
                       text = element_text(size=14)))


plot_alpha
ggdraw() + 
  draw_plot(plot_alpha + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 
#draw_plot_label(c("A", "B"), c(0, 1.5), c(.75, 0.75), size = 12)

######################################################
#############Plots for thesis proposal################
######################################################
# my_mean.alpha = aggregate(Mean_alpha_bin$rarefied.alpha , by=list(Mean_alpha_bin$Yr_bin) , mean); colnames(my_mean.alpha)=c("Yr_Bin" , "mean")
# se <- function(x) sqrt(var(x)/length(x))
# my_sd_alpha=aggregate(Mean_alpha_bin$rarefied.alpha , by=list(Mean_alpha_bin$Yr_bin) , se); colnames(my_sd_alpha)=c("Yr_Bin" , "sd")
# my_info.rarefied=merge(my_mean.alpha , my_sd_alpha , by.x=1 , by.y=1)
# 
# 
# box4 <- ggplot(my_info.rarefied) + 
#   # geom_point(aes(x = Yr_Bin, y = CMRL) , colour=rgb(0.8,0.7,0.1,0.4) , size=5) + 
#   geom_point(data = my_info.rarefied, aes(x=Yr_Bin , y = mean) , colour = rgb(0.6,0.5,0.4,0.7) , size = 8) +
#   geom_errorbar(data = my_info.rarefied, aes(x = Yr_Bin, y = sd, ymin = mean - sd, ymax = mean + sd), colour = rgb(0.4,0.8,0.2,0.4) , width = 0.7 , size=1.5) +
#   scale_y_continuous(name = "Rarefied Alpha-diversity", breaks = seq(37, 47, 1)) +
#   scale_x_continuous(name = "Year Bin", breaks = seq(1, 8, 1), 
#                      labels = c("1980 - 1984", "1985 - 1989", "1990 - 1994", "1995 - 1999", 
#                                 "2000 - 2004", "2005 - 2009", "2010 - 2014", "2015 - 2019")) +
#   coord_cartesian(ylim = c(35, 50)) 
# 
# box4
######################################################




#Environmental Data Loading and Cleaning 
#PRISM Data Script 
prism.1 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_198001_199412.csv"))
prism.2 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_199501_200912.csv"))
prism.3 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_201001_201712.csv"))
prism.4 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_198001_199412.csv"))
prism.5 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_199501_200912.csv"))
prism.6 <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_201001_201712.csv"))
prism.1.monthly <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1966_GoM_mids_monthly_normals.csv")) 
prism.2.monthly <- read.csv(here::here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_monthly_normals.csv"))

#Bind the datasets together 
prism.1966 <- do.call("rbind", list(prism.1, prism.2, prism.3))
prism.1966$gis_layer <- 1966
prism.1966 <- prism.1966[prism.1966$Name == "SEMINOLE HLS" | prism.1966$Name == "ROSELAND" | prism.1966$Name == "BENNDALE" | prism.1966$Name == "MOON LAKE" | prism.1966 == "KINGSVILLE" | prism.1966 == "RAYMONDVILLE", ]
prism.1966 <- prism.1966[complete.cases(prism.1966),]
prism.1999 <- do.call("rbind", list(prism.4, prism.5, prism.6))
prism.1999$gis_layer <- 1999
prism.total <- rbind(prism.1966, prism.1999)
prism.total <- prism.total[complete.cases(prism.total),]
prism.total <- prism.total[!duplicated(prism.total),]
prism.total <- separate(prism.total, 5, c("Year", "Month"), sep = "-", remove = F)
prism.total <- prism.total %>% unite(time, c("Year", "Month"), sep = "", remove = F)

#Rename variables, wierd CSV notation
colnames(prism.total)[colnames(prism.total) == "ppt..inches."] <- "precip"
colnames(prism.total)[colnames(prism.total) == "tmax..degrees.F."] <- "tmax"
colnames(prism.total)[colnames(prism.total) == "tmin..degrees.F."] <- "tmin"
colnames(prism.total)[colnames(prism.total) == "tmean..degrees.F."] <- "tmean"

#Year bins for analysis 
prism.total$Year <- as.numeric(prism.total$Year)
prism.total$yr_bin <- 1
prism.total$yr_bin[prism.total$Year >= 1980 & prism.total$Year <= 1984] <- 1
prism.total$yr_bin[prism.total$Year >= 1985 & prism.total$Year <= 1989] <- 2
prism.total$yr_bin[prism.total$Year >= 1990 & prism.total$Year <= 1994] <- 3
prism.total$yr_bin[prism.total$Year >= 1995 & prism.total$Year <= 1999] <- 4
prism.total$yr_bin[prism.total$Year >= 2000 & prism.total$Year <= 2004] <- 5
prism.total$yr_bin[prism.total$Year >= 2005 & prism.total$Year <= 2009] <- 6
prism.total$yr_bin[prism.total$Year >= 2010 & prism.total$Year <= 2014] <- 7
prism.total$yr_bin[prism.total$Year >= 2015 & prism.total$Year <= 2018] <- 8

#Bin months for seasons
prism.total$Month <- as.numeric(prism.total$Month)
prism.total$season <- "Winter"
prism.total$season[prism.total$Month >= 03 & prism.total$Month <= 05] <- "Spring"
prism.total$season[prism.total$Month >= 06 & prism.total$Month <= 08] <- "Summer"
prism.total$season[prism.total$Month >= 09 & prism.total$Month <= 11] <- "Fall"
prism.total$season[prism.total$Month >= 12 & prism.total$Month <= 02] <- "Winter"

#Reduce data for just Spring and Summer seasons
prism.total <- prism.total[prism.total$season == "Spring" | prism.total$season == "Summer", ]

prism.simple <- prism.total[, -4:-6]
colnames(prism.simple)[colnames(prism.simple) == "Name"] <- "Site"
prism.simple$Site <- as.character(prism.simple$Site)

#Calculate Monthly Averages in Temperature for bins
#Might be able to integrate precip in here together
avgs <- aggregate(tmean ~ Site +  Month + yr_bin, prism.simple, FUN = "mean")
avg.max <- aggregate(tmax ~ Site + Month + yr_bin, prism.simple, FUN = "mean")
avg.max <- avg.max[, "tmax"]
avg.max <- as.data.frame(avg.max)
avg.min <- aggregate(tmin ~ Site + Month + yr_bin, prism.simple, FUN = "mean")
avg.min <- avg.min[, "tmin"]
avg.min <- as.data.frame(avg.min)
precipitation <- aggregate(precip ~ Site + Month + yr_bin, prism.simple, FUN = "mean")
precipitation <- precipitation[, "precip"]
precipitation <- as.data.frame(precipitation)
avgs <- do.call("cbind", list(avgs, avg.max, avg.min, precipitation))
avgs$site.mo.yr <- paste0(avgs$Site, "_", avgs$Month, "_", avgs$yr_bin)



#Binning the beta matrix data 
beta.matrix$unique_id <- as.character(beta.matrix$unique_id)
betas <- beta.matrix %>% separate(unique_id, c("rteno", "segment", "year"), sep = "_") %>% unite(site, c("rteno", "segment"), sep = "_")
betas$year <- as.numeric(betas$year)
betas$Yr_bin <- 1
betas$Yr_bin[betas$year >= 1985 & betas$year <= 1989] <- 2
betas$Yr_bin[betas$year >= 1990 & betas$year <= 1994] <- 3
betas$Yr_bin[betas$year >= 1995 & betas$year <= 1999] <- 4
betas$Yr_bin[betas$year >= 2000 & betas$year <= 2004] <- 5
betas$Yr_bin[betas$year >= 2005 & betas$year <= 2009] <- 6
betas$Yr_bin[betas$year >= 2010 & betas$year <= 2014] <- 7
betas$Yr_bin[betas$year >= 2015 & betas$year <= 2018] <- 8

#Getting the Betas into the environmental Data 
#Pull the minumum year-bin for the linking up the 
#climate data to the first year bin with a beta value
min.betas <- aggregate(Yr_bin ~ site, betas, FUN = "min")

betas.mod <- betas

#Calculate Mean Beta for each year-bin 
betas.mod <- aggregate(beta ~ site + Yr_bin, betas, FUN = "mean")
bbs_site_temp <- bbs_total[, c("rteno.x", "rtename", "site", "rt_yr")]
bbs_site_temp$rtename <- as.character(bbs_site_temp$rtename)
bbs_site_temp <- bbs_site_temp[!duplicated(bbs_site_temp$site), ]
colnames(bbs_site_temp)[colnames(bbs_site_temp) == "rtename"] <- "Site"

#Create a site.link dataframe 
#site.merge <- bbs_site_temp[, c("site", "", "rteno")]
site.merge <- bbs_site_temp
site.merge <- site.merge[!duplicated(site.merge),]

site.merge$Site <- gsub(" ", "_", site.merge$Site)
#site.merge[site.merge$rtename == "BETAIR", 2] <- "BELAIR"
#site.merge[site.merge$rtename == "L._ATASCOSA_NWR", 2] <- "L__ATASCOSA_NWR"
site.merge[site.merge$Site == "SUGARLOAF_KEY_2", 2] <- "SUGARLOAF_KE"
site.merge[site.merge$Site == "EGLIN_A.F.B.", 2] <- "EGLIN_A_F_B_"
site.merge[site.merge$Site == "SEMINOLE_HLS", 2] <- "SEMINOLE_HILLS"
site.merge[site.merge$Site == "ALABAMA_PORT", 2] <- "DAUPHIN_IS_2"





min.betas <- merge(min.betas, site.merge, by = "site")

colnames(min.betas)[colnames(min.betas) == "Yr_bin"] <- "min.yr.bin"

#Merge the climate averages with the site data for ref
avgs$Site <- gsub(" ", "_", avgs$Site)
#avgs[avgs$Site == "BETAIR", 1] <- "BELAIR"
#avgs[avgs$Site == "L._ATASCOSA_NWR", 1] <- "L__ATASCOSA_NWR"
avgs[avgs$Site == "SUGARLOAF_KE", 1] <- "SUGARLOAF_KE"
avgs[avgs$Site == "EGLIN_A.F.B.", 1] <- "EGLIN_A_F_B_"
avgs[avgs$Site == "SEMINOLE_HLS", 1] <- "SEMINOLE_HILLS"
avgs[avgs$Site == "ALABAMA_PORT", 1] <- "DAUPHIN_IS_2"
avgs[avgs$Site == "FORKED_ISLAN", 1] <- "FORKED_ISLAND"

avgs <- merge(avgs, site.merge, by = "Site")
avgs$unique_id <- paste0(avgs$site, "_", avgs$yr_bin)
avgs <- merge(avgs, min.betas, by = "site")
avgs <- avgs[avgs$yr_bin >= avgs$min.yr.bin,]

avgs <- avgs %>% select(site, Site.x, Month, yr_bin, tmean, avg.max, avg.min,
                        precipitation, site.mo.yr, rteno.x.x, min.yr.bin)

avgs$site_bin <- paste0(avgs$site, "_", avgs$yr_bin)
colnames(avgs)[colnames(avgs) == "Site.x"] <- "Rt_name"
colnames(avgs)[colnames(avgs) == "rteno.x.x"] <- "rteno"


#Metric Conversion
avgs <- avgs %>% mutate(tmean.c = (tmean - 32) / 1.8, 
                        tmax.c = (avg.max - 32) / 1.8,
                        tmin.c = (avg.min - 32) / 1.8,
                        precipitation.c = (precipitation * 2.54))

#Calculated anomalies relative to the base year
avgs <- avgs %>% group_by(site, Month) %>% 
  arrange(yr_bin, .by_group = T) %>% 
  mutate(anomalies = tmean.c - first(tmean.c)) %>% 
  mutate(max.anomalies = tmax.c - first(tmax.c)) %>% 
  mutate(min.anomalies = tmin.c - first(tmin.c)) %>%
  mutate(precip.anomalies = precipitation.c - first(precipitation.c))

#Same as above for anomalies, need to merge the data for betas & anomalies 
#using the join field of State_Rteno

temp.agg <- avgs %>% ungroup(avgs) %>% select(site, yr_bin, tmean.c, tmax.c, tmin.c, precipitation.c,
                            anomalies, max.anomalies, min.anomalies, precip.anomalies)


#Aggregate all the anomalies by year_bin for each site
temp.agg <- aggregate(. ~ site + yr_bin, temp.agg, FUN = "mean")
temp.avgs <- temp.agg
betas.mod$unique_id <- paste0(betas.mod$site, "_", betas.mod$Yr_bin)
temp.avgs$unique_id <- paste0(temp.avgs$site, "_", temp.avgs$yr_bin)

#Put betas with temperature data 
lmer.df <- merge(temp.avgs, betas.mod, by = "unique_id")
test <- lmer.df

ggplotRegression(tempmod)

ggplot(temp.avgs, aes(x = yr_bin, y = anomalies, group = site)) + labs(x = "Year Bin", ylab = "Temperature Anomalies") + geom_line()
ggplot(temp.avgs, aes(x = yr_bin, y = precip.anomalies, group = site)) + labs(x = "Year Bin", ylab = "Temperature Anomalies") + geom_line()

tempmod <- lm(temp.avgs$anomalies ~ temp.avgs$yr_bin)
summary(tempmod)


##Create df to link diversity, climate, and LULC data together## 

true.sites <- final_sp_df[, c("rteno.x", "site")]
true.sites <- true.sites[!duplicated(true.sites),]
true.sites <- true.sites$site
true.sites <- wetland.link[wetland.link$rteno.segment %in% true.sites,]

#write.csv(true.sites, here::here("Data_BBS/Generated DFs/LULC_Final_Sites_DF.csv"))

######################################################################
##############Script cleaned to here 7/1/2019#########################
######################################################################







#Average the Richness at Sites in the year bins for use as a covariate
alpha.1980 <- bbs_total
alpha.1980 <- alpha.1980[alpha.1980$Year >= 1980,]
alpha.1980$yr_bin <- 1
alpha.1980$yr_bin[alpha.1980$Year >= 1985 & alpha.1980$Year <= 1989] <- 2
alpha.1980$yr_bin[alpha.1980$Year >= 1990 & alpha.1980$Year <= 1994] <- 3
alpha.1980$yr_bin[alpha.1980$Year >= 1995 & alpha.1980$Year <= 1999] <- 4
alpha.1980$yr_bin[alpha.1980$Year >= 2000 & alpha.1980$Year <= 2004] <- 5
alpha.1980$yr_bin[alpha.1980$Year >= 2005 & alpha.1980$Year <= 2009] <- 6
alpha.1980$yr_bin[alpha.1980$Year >= 2010 & alpha.1980$Year <= 2014] <- 7
alpha.1980$yr_bin[alpha.1980$Year >= 2015 & alpha.1980$Year <= 2020] <- 8
#alpha.1980 <- separate(alpha.1980, unique_id, c("State", "rteno"), sep = "_")
rownames(alpha.1980) <- c()
colnames(alpha.1980)[colnames(alpha.1980) == "abbrev"] <- "site"
alpha.1980 <- alpha.1980[, c("site", "yr_bin", "Year", "Total_Abundance", "Total_SR", "rarefied.alpha")]
alpha.agg <- aggregate(. ~ site + yr_bin, alpha.1980, FUN = "mean")
alpha.agg$unique_id <- paste0(alpha.agg$site, "_", alpha.agg$yr_bin)

test <- merge(test, alpha.agg, "unique_id")
test <- test %>% group_by(site) %>% arrange(Yr_bin) %>% 
  mutate(p.change.alpha = (rarefied.alpha - first(rarefied.alpha)) / first(rarefied.alpha) * 100) %>%
  ungroup(test)
colnames(test)[colnames(test) == "site.x"] <- "Site"
test <- test[, c("unique_id", "Site", "Yr_bin", "Year", "rarefied.alpha", 
                 "p.change.alpha", "anomalies", "max.anomalies", "min.anomalies", "precip.anomalies", 
                 "beta")]

#Remove the first time bin for each site, that way we don't have a lot of 
#untrue 0s in the tests, this is only for p.change
test1 <- test %>%
  group_by(Site) %>%
  slice(2:n())

test <- separate(test, Site, into = c("State", "Rteno"), sep = "_")
test$unique.id <- paste0(test$Rteno, "_", test$Yr_bin)
#1og-likelihood 
mod.alpha <- lmer(p.change.alpha ~ precip.anomalies + (1|Site), data = test1, REML = F)
summary(mod.alpha)
Anova(mod.alpha)
eta_sq(mod.alpha)
t.test(test1$p.change.alpha, mu = 0)
r.squaredGLMM(mod.alpha)
plot(test1$precip.anomalies, test1$p.change.alpha)
abline(lm(test1$p.change.alpha ~ test1$precip.anomalies))

mod.beta <- lmer(beta ~ anomalies + rarefied.alpha + precip.anomalies + (1|Rteno), data = test)
summary(mod.beta)
tab_model(mod.beta)



##########################################################################
###################Code Cleaned Up To Here 4/24/2019######################
##########################################################################




lmer.df <- separate(lmer.df, col = site.x, into = c("State", "Rteno"),sep = "_")
lmer.df$unique_id <- paste0(lmer.df$Rteno, "_", lmer.df$yr_bin)
lmer2 <- lmer.df %>% group_by(Rteno) %>% first(lmer.df, order_by = "yr_bin")

mod.complete <- merge(lmer.df, lulc.agg, by.x = "unique_id", by.y = "unique.id")

mod.lulc <- lmer(data = lmer.df, beta ~ anomalies + (1|Rteno))
summary(mod.lulc)
tab_model(mod.lulc)

lmer.df <- 
  
  a.test <- test1[, c("site", "Yr_bin", "p.change.alpha", "rarefied.alpha", "anomalies", "precip.anomalies", "min.anomalies", "max.anomalies")]
a.test <- a.test[a.test$p.change.alpha != 0, ]
s <- scale(a.test$p.change.alpha)
GHQ <- lmer(p.change.alpha ~ anomalies + precip.anomalies + (1 | Site), data = a.test) # Set nAGQ to # of desired iterations
summary(GHQ)

##Get residuals for partial regression plots 
ggplot (data = prplot, aes(x = precip.anomalies, y = y_partial)) + geom_smooth(method = "lm")
devtools::install_github("hohenstein/remef")
library(remef)
beta.t.anom <- remef(mod.2, fix = "precip.anomalies", ran = "all")
beta.t.anom <- as.data.frame(beta.t.anom)
prplot <- test
prplot <- cbind(prplot, beta.t.anom)
precip.partial <- remef(mod2, fix = "anomalies", ran = "all")
precip.partial <- as.data.frame(precip.partial)
prplot <- cbind(prplot, precip.partial)
alpha.partial <- 
  ggplot (data = prplot, aes(x = anomalies, y = beta.t.anom)) + geom_smooth(method = "lm")

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), x = "Year Bins", y = ("Temperature Anomalies (Degree C)"))}

temperature <- lm(y_partial ~ anomalies, data = prplot)
ggplotRegression(temperature) + labs(x.axis = "Rarefied Alpha")

precip <- lm(precip.partial ~ precip.anomalies, data = prplot)
ggplotRegression(precip)

##If the model include rarefied alpha need to fix for all three effects
mod.fix <- lmer(beta ~ anomalies + rarefied.alpha + precip.anomalies + (1|Rteno), data = test, REML = F)
tab_model(mod.fix)
r.squaredGLMM(mod.fix)
summary(mod.fix)
ano <- Anova(mod.fix)
fix <- anova_stats(car::Anova(mod.fix, type = 2))
eta_sq(mod.fix)
predictors <- test[, c("anomalies", "max.anomalies", "min.anomalies", "precip.anomalies", "rarefied.alpha")]
cor.test <- round(cor(predictors, use = "pair"), 2)
print(cor.test)
r.squaredGLMM(mod.fix)
est <- effect("anomalies", partial.residuals = T, mod.fix)
plot(est)
summary(est)
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}
r2.corr.mer(mod.fix)

test$Years <- 1
test$Years[test$Yr_bin == 2] <- "1985-1989"
test$Years[test$Yr_bin == 3] <- "1990-1994"
test$Years[test$Yr_bin == 4] <- "1995-1999"
test$Years[test$Yr_bin == 5] <- "2000-2004"
test$Years[test$Yr_bin == 6] <- "2005-2009"
test$Years[test$Yr_bin == 7] <- "2010-2014"
test$Years[test$Yr_bin == 8] <- "2015-2018"
test$Years <- as.character(test$Years)

test <- separate(test, unique_id, c("State", "RTENO", "waste"), sep = "_")

beta.plot <- (ggplot(test, aes(x = Years, y = beta, group = site.x, 
                               colour = factor(BCR, 
                                               labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Tamaulipan Brushlands", "Gulf Coastal Prairie")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_line(aes(x = Years, y = beta, group = site.x, 
                              colour = factor(BCR, 
                                              labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Tamaulipan Brushlands", "Gulf Coastal Prairie"))), linetype = 1, size = 1) +
                xlab("Years") +
                ylab(expression(paste(beta, "-Diversity"))) +
                labs(colour = "Bird Conservation Region") + 
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(vjust = -1, size = 14),
                      axis.title.y = element_text(vjust = 1.5, size = 14),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=14)) +
                scale_color_manual(values = cbbPalette))
beta.plot 
cbbPalette <- c("#99FFCC","#99FFFF", "#66CCCC", "#339999", "#006666")
cbbPalette <- c("#000000", "#000000", "#000000", "#000000", "#000000")
#test$scaled.anomalies <- scale(test$anomalies, center = F, scale = apply(test$anomalies, 2, sd, na.rm = T))

beta.t.anom <- remef(mod.fix, fix = c("precip.anomalies", "rarefied.alpha"), ran = "all")
beta.t.anom <- as.data.frame(beta.t.anom)
prplot <- test
prplot <- cbind(prplot, beta.t.anom)
precip.partial <- remef(mod.fix, fix = c("anomalies", "rarefied.alpha"), ran = "all")
precip.partial <- as.data.frame(precip.partial)
prplot <- cbind(prplot, precip.partial)
alpha.partial <- remef(mod.fix, fix = c("anomalies", "precip.anomalies"), ran = "all")
prplot <- cbind(prplot, alpha.partial)

temp.fix <- lm(beta.t.anom ~ scaled.anomalies, data = prplot)
ggplotRegression(temp.fix) 

precip.fix <- lm(precip.partial ~ precip.anomalies, data = prplot)
ggplotRegression(precip.fix)

alpha.fix <- lm(alpha.partial ~ rarefied.alpha, data = prplot)
ggplotRegression(alpha.fix)



#Stats and Model Testing 
#Testing to structure of the data
#Nomal - Best so far
lmer.df$beta.t <- lmer.df$beta + 1
qqp(lmer.df$beta.t, "norm")

#Log-normal distribution - No 
qqp(lmer.df$beta.t, "lnorm")
qqp(test$p.change.alpha, "lnorm")

#Negative binomial distribution - No - not discrete data
nbinom <- fitdistr(lmer.df$beta.t, "Negative Binomial")
qqp(lmer.df$beta.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#poisson dist - No
poisson <- fitdistr(a.test$p.change.alpha, "Poisson")
qqp(a.test$p.change.alpha, "pois", lambda = poisson$estimate)

#gamma - Similar to the normal distribution...between the two
#I wouldrather us ethe normal dist for simplicity 
gamma <- fitdistr(a.test$p.change.alpha, "gamma")
qqp(lmer.df$beta.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

mod2 <- lmer(beta ~ anomalies + precip.anomalies + (1|site.x), data = lmer.df, REML = F)
summary(mod2)
Anova(mod2)
r.squaredGLMM(mod2)

etaSquared(mod2, type = 2, anova = F)
predictors <- lmer.df[, c("anomalies", "max.anomalies", "min.anomalies", "precip.anomalies")]
cor.test <- round(cor(predictors, use = "pair"), 2)
print(cor.test)

hist(lmer.df$anomalies)

plot(lmer.df$yr_bin, lmer.df$anomalies)
t.test(lmer.df$anomalies, mu = 0)
t.test(lmer.df$max.anomalies, mu = 0)
t.test(lmer.df$precip.anomalies, mu = 0)

r.squaredGLMM(mod2)
#min.anomalies .88 correlation with anomalies

plot(lmer.df$anomalies, lmer.df$beta)

#Ploting all the lines 
ggplot(data = lmer.df, aes(x = yr_bin, y = beta))+
  geom_smooth()+
  theme(legend.position = "none")


lmer.df.2015 <- lmer.df[lmer.df$yr_bin == 8,]
mod.2015 <- lm(beta ~ anomalies + precip.anomalies, data = lmer.df.2015)
summary(mod.2015)
Anova(mod.2015)

##Code for getting a VIF 
vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- Yr_Bin(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  Yr_Bin(v) <- nam
  v }
vif.lme(mod2)

#Figures 
richness <- aggregate(rarefied.alpha ~ abbrev, site_data_merge, FUN = "mean")
spat <- site_data_merge[, c("abbrev", "latitude", "longitude")]
spat <- spat[!duplicated(spat),]
richness <- merge(richness, spat, by = "abbrev")
ggplot(data = richness, aes(x = latitude, y = rarefied.alpha)) + geom_smooth() + theme(legend.position = "none")
r <- lm(rarefied.alpha ~ latitude, data = richness)
summary(r)
ggplot(data = richness, aes(x = longitude, y = rarefied.alpha)) + geom_smooth() + theme(legend.position = "none")
plot(richness$latitude, richness$rarefied.alpha)
abline(lm(richness$rarefied.alpha ~ richness$latitude))
#This is likely a relic of the first attempt being messy and incohesiveness of two scripts 
# # prcip.agg <- aggregate(anomalies.pcp ~ Site + yr_bin, avgs.pcp, FUN = "mean")
# # prcip.avgs <- prcip.agg
# 
# bbs_site_temp <- merge(temp.avgs, bbs_site_temp, by = "Site", all = F)
# bbs_site_temp$unique_id <- paste0(bbs_site_temp$Site, "_", bbs_site_temp$yr_bin)
# # prcip.avgs$unique_id <- paste0(prcip.avgs$Site, "_", prcip.avgs$yr_bin)
# # bbs_site_temp <- merge(prcip.avgs, bbs_site_temp, by = "unique_id", all = F)
# 
# bbs_site_temp <- bbs_site_temp[, c(-5:-6)]
# colnames(bbs_site_temp)[colnames(bbs_site_temp) == "Site.x"] <- "Site"
# colnames(bbs_site_temp)[colnames(bbs_site_temp) == "yr_bin.x"] <- "yr_bin"
# bbs_site_temp <- separate(bbs_site_temp, unique_id, c("name", "waste"), sep = "_")
# bbs_site_temp <- bbs_site_temp[, c(-1:-2)]
# bbs_site_temp$unique_id <- paste0(bbs_site_temp$state, "_", bbs_site_temp$rteno, "_", bbs_site_temp$yr_bin)
# #bbs_site_temp <- bbs_site_temp[, c(3, 4, 5, 6, 8, 9, 10, 13)]
# colnames(betas.mod)[colnames(betas.mod) == "site"] <- "sites"
# betas.mod$sites <- as.character(betas.mod$sites)
# betas.mod$unique_id <- paste0(betas.mod$sites, "_", betas.mod$Yr_bide

#Partial Regression Plots
res.beta <- residuals(lm(beta ~ anomalies + pct.ur, data = test1))
res.temp <- residuals(lm(anomalies ~ pct.ur + rarefied.alpha, data = test1))
res <- cbind(as.data.frame(res.beta), as.data.frame(res.temp))
res.beta.u <- residuals(lm(beta ~ pct.ur + rarefied.alpha, data = test1))
res.ur <- residuals(lm(pct.ur ~ anomalies + rarefied.alpha, data = test1))
res2 <- cbind(as.data.frame(res.beta.u), as.data.frame(res.ur))
res.beta.a <- residuals(lm(beta ~ anomalies + pct.ur, data = test1))
res.alpha <- residuals(lm(rarefied.alpha ~ anomalies + pct.ur, data = test1))
res3 <- cbind(as.data.frame(res.beta.a), as.data.frame(res.alpha))
plot(res.temp, res.beta)
k <- lm(res.plot$res.beta ~ res.plot$res.temp)
p <- lm(res.plot$res.beta.u ~ res.plot$res.ur)
g <- lm(res.plot$res.beta.a ~ res.plot$res.alpha)
ggplotRegression(k)
ggplotRegression(p)
ggplotRegression(g)

res.plot <- cbind(prplot, res)
res.plot <- cbind(res.plot, res2)
res.plot <- cbind(res.plot, res3)

etaSquared(mod.fix, type = 2, anova = F)

t <- lm(beta ~ anomalies, data = test)
ggplotRegression(t)
p <- lm(beta ~ precip.anomalies, data = test)
ggplotRegression(p)
a <- lm(beta ~ rarefied.alpha, data = test)
ggplotRegression(a)
# plot(res$res.temp, res$res.beta)
# abline(lm(res$res.beta ~ res$res.temp))
# l <- lm(res$res.beta ~ res$res.temp)
# summary(l)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = Yr_Bin(fit$model)[2], y = Yr_Bin(fit$model)[1])) + 
    geom_point() +
    geom_smooth(method = "lm", col = "red") +
    labs(x = expression(paste("Precipitation Anomaly ", (cm))), y = expression(paste(beta, " -Diversity")))}
# title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
#                  # "Intercept =",signif(fit$coef[[1]],5 ),
#                  " Slope =",signif(fit$coef[[2]], 5),
#                  " P =",signif(summary(fit)$coef[2,4], 5)), x = expression(paste("Rarefied ", alpha, "-Diversity")), y = expression(paste(beta, " -Diversity")))}


##Extracting BCR for figures
bcr <- read.csv(here::here("Data_BBS/States_GoM/bbsrts_bcr.csv"))
bcr <- bcr[, c("RTENO", "BCR")]
bcr$RTENO <- as.numeric(bcr$RTENO)
bcr$BCR <- as.numeric(bcr$BCR)
test$RTENO <- as.numeric(test$RTENO)
bcr <- bcr[!duplicated(bcr$RTENO), ]
test <- merge(test, bcr, by = "RTENO")


#Table of output 
tab_model(mod.fix, p.val = "wald", dv.labels = "Beta Diversity", pred.labels = c("Intercept", "Temperature Anomalies", "Precipitation Anomalies", "Alpha Diversity"), title = "Linear Mixed-Effects Model Output", digits = 4, digits.p = 4, p.style = "asterisk", file = here::here("Figures and Tables"))

#Getting Dataset to make a map for route locations 
climate.map <- separate(slopes, site.list, c("State", "rteno"), sep = "_")
climate.map <- merge(climate.mao, rts, by = "rteno")

write.csv(climate.map, file = here::here("climate.map.csv"))



#######################################################################
######Link Spacies Name with Range and looking at tropicalization######
#######################################################################

bbs_total$english_common_name <- as.character(bbs_total$english_common_name)
bbs_total <- merge(bbs_total, range.merge, by.x = "english_common_name", by.y = "English_Common_Name")
bbs.agg <- aggregate(data = bbs_total, Latitude ~ abbrev + year, FUN = mean)
bbs.agg <- bbs.agg[bbs.agg$year >= 1980, ]
bbs.agg$Yr_bin <- 1
bbs.agg$Yr_bin[bbs.agg$year >= 1985 & bbs.agg$year <= 1989] <- 2
bbs.agg$Yr_bin[bbs.agg$year >= 1990 & bbs.agg$year <= 1994] <- 3
bbs.agg$Yr_bin[bbs.agg$year >= 1995 & bbs.agg$year <= 1999] <- 4
bbs.agg$Yr_bin[bbs.agg$year >= 2000 & bbs.agg$year <= 2004] <- 5
bbs.agg$Yr_bin[bbs.agg$year >= 2005 & bbs.agg$year <= 2009] <- 6
bbs.agg$Yr_bin[bbs.agg$year >= 2010 & bbs.agg$year <= 2014] <- 7
bbs.agg$Yr_bin[bbs.agg$year >= 2015 & bbs.agg$year <= 2018] <- 8
bbs.agg$unique_id <- paste0(bbs.agg$abbrev, "_", bbs.agg$Yr_bin)
betas$unique_id <- paste0(betas$site, "_", betas$Yr_bin) 
betas <- betas[, c(-1, -2, -4)]
bbs.agg <- merge(bbs.agg, betas, by = "unique_id")

colnames(bbs.agg)[colnames(bbs.agg) == "Latitude"] <- "CMRL"

##Making some figures for tropicalization data 


line = "#0000CC"
fill = "#6699CC"
box1 <- ggplot(bbs.agg, aes(x = Yr_bin, y = CMRL, group = Yr_bin)) +
  geom_boxplot(alpha = 0.3, notch = T, fill = fill, colour = line) + 
  scale_y_continuous(name = "Community Mean Range Limit", breaks = seq(46, 70, 2)) +
  scale_x_continuous(name = "Year Bin (5 years from baseline)") + 
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Dark2") + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange") +
  stat_summary(fun.y = mean,
               geom = "line")
box1


range_t <- bbs.agg %>% group_by(abbrev) %>% arrange(Yr_bin) %>%
  filter(row_number()==1 | row_number()==n()) 

box2 <- ggplot(df.mer, aes(x = Yr_bin, y = CMRL, group = Yr_bin, fill = region)) +
  geom_boxplot(aes(x = Yr_bin, y = CMRL, group = Yr_bin, fill = region)) + 
  scale_y_continuous(name = "Community Mean Range Limit", breaks = seq(46, 70, 2)) +
  scale_x_continuous(name = "Year Bin (5 years from baseline)") +
  theme(legend.position = "none") 


box2

#Calculation of mean and sd of each group ?
my_mean=aggregate(bbs.agg$CMRL , by=list(bbs.agg$Yr_bin) , mean); colnames(my_mean)=c("Yr_Bin" , "mean")
se <- function(x) sqrt(var(x)/length(x))
my_sd=aggregate(bbs.agg$CMRL , by=list(bbs.agg$Yr_bin) , se); colnames(my_sd)=c("Yr_Bin" , "sd")
my_info=merge(my_mean , my_sd , by.x=1 , by.y=1)


#Plotting beta and cmrl
plot5 <- ggplot(bbs.agg, aes(x = beta, y = CMRL, group = abbrev)) +
  geom_smooth(method = "lm", se = F)
plot5
# Make the plot
coefs <- coef(lm(df.mer$CMRL ~ df.mer$Yr_bin))
box3 <- ggplot(my_info) + 
  # geom_point(aes(x = Yr_Bin, y = CMRL) , colour=rgb(0.8,0.7,0.1,0.4) , size=5) + 
  geom_point(data = my_info, aes(x=Yr_Bin , y = mean) , colour = rgb(0.6,0.5,0.4,0.7) , size = 8) +
  geom_errorbar(data = my_info, aes(x = Yr_Bin, y = sd, ymin = mean - sd, ymax = mean + sd), colour = rgb(0.4,0.8,0.2,0.4) , width = 0.7 , size=1.5) +
  scale_y_continuous(name = "Community Mean Range Limit", breaks = seq(59, 62, 0.5)) +
  scale_x_continuous(name = "Year Bin", breaks = seq(1, 8, 1), 
                     labels = c("1980 - 1984", "1985 - 1989", "1990 - 1994", "1995 - 1999", 
                                "2000 - 2004", "2005 - 2009", "2010 - 2014", "2015 - 2019")) +
  coord_cartesian(ylim = c(59, 62)) 

box3 + annotate("text", x = 6, y = 61.8, size = 6, label = "Average difference between year bins corresponds 
                to a decrease in community latitude of 1.6 km/y")


mean.change <- my_info %>% mutate(lagged = lag(mean)) %>% mutate(dif = mean - lagged) %>% 
  mutate(avg.dif = mean(dif, na.rm = T)) 

mod5 <- lm(range_t$CMRL ~ range_t$Yr_bin)
summary(mod5)
Anova(mod5)

mod.range <- lm(bbs.agg$CMRL ~ bbs.agg$beta)
summary(mod.range)
Anova(mod.range)


bbs.testing <-  aggregate(data = bbs.agg, CMRL ~ abbrev + Yr_bin, FUN = mean)
bbs.testing <- separate(bbs.testing, abbrev, into = c("State", "rteno"), sep = "_")
bbs.testing$abbrev <- paste0(bbs.testing$State, "_", bbs.testing$rteno)
link <- link.me.up
link <- link %>% mutate_if(is.integer, as.character) 

df.mer <- merge(bbs.testing, link, by = "rteno")
df.mer$region <- "North"
df.mer$region[df.mer$latitude <= 28.27796] <- "South"

plot_cmrl2 <- ggplot(data = df.mer, aes(x = Yr_bin, y = CMRL, group = region)) +
  geom_smooth(method = "lm", se = F)
plot_cmrl2

#Plotting 
plot_cmrl <- ggplot(data = df.mer, aes(x = Yr_bin, y = CMRL, colour = region)) +
  geom_boxplot(aes(x = Yr_bin, y = CMRL, colour = region))
plot_cmrl

bbs.agg <- bbs.agg[complete.cases(bbs.agg),]

range.exp <- lmer(bbs.agg$beta ~ bbs.agg$CMRL + (1|bbs.agg$abbrev))
summary(range.exp)
tab_model(range.exp)
Anova(range.exp)
plot(bbs.agg$beta ~ bbs.agg$CMRL)
abline(lm(bbs.agg$beta ~ bbs.agg$CMRL))
beta.matrix <- merge()