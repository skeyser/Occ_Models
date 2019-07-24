#LULc Data Processing and Cleaning 
#install.packages("pacman")
library("pacman")
pacman::p_load("here", "tidyverse", "reshape2", "stringr")


#Load in and tidy up data for link rtnames with rtenos
all.routes <- read.csv(here("Data_BBS/States_GoM/routes.csv"))
Gom.rt.id <- read.csv(here("Data_BBS/States_GoM/Gom_Routes_ID.correct.csv"))
link.to.lulc <- Gom.rt.id[, c("rteno", "rtename")]
link.to.lulc$rtename <- as.character(link.to.lulc$rtename)
link.to.lulc$rtename <- gsub(" ", "_", link.to.lulc$rtename)
link.to.lulc <- link.to.lulc[!duplicated(link.to.lulc),]

link.to.lulc[link.to.lulc$rtename == "BETAIR", 2] <- "BELAIR"
link.to.lulc[link.to.lulc$rtename == "L._ATASCOSA_NWR", 2] <- "L__ATASCOSA_NWR"
link.to.lulc[link.to.lulc$rtename == "SUGARLOAF_KEY_2", 2] <- "SUGARLOAF_KEY"
link.to.lulc[link.to.lulc$rtename == "EGLIN_A.F.B.", 2] <- "EGLIN_A_F_B_"
link.to.lulc[link.to.lulc$rtename == "SEMINOLE_HLS", 2] <- "SEMINOLE_HILLS"
link.to.lulc[link.to.lulc$rtename == "ALABAMA_PORT", 2] <- "DAUPHIN_IS_2"


#Create list for file names to read bulk files for NLCD, CCAP, and Mangrove
file_list <- list.files(path=here("Data_Envi/BBS_CCAP_50stop"), pattern="*.csv", full.names = T) # create list of all .csv files in folder
file_list_nlcd <- list.files(path=here("Data_Envi/BBS_NLCD_50stop"), pattern="*.csv", full.names = T) # create list of all .csv files in folder
file_list_man <- list.files(path=here("Data_Envi/BBS_Mangrove_50stop"), pattern="*.csv", full.names = T) # create list of all .csv files in folder

#Read all files in as list 
ccap.files <- lapply(file_list, function(x) read.csv(x, stringsAsFactors = F))
nlcd.files <- lapply(file_list_nlcd, function(x) read.csv(x, stringsAsFactors = F))
man.files <- lapply(file_list_man, function(x) read.csv(x, stringsAsFactors = F))


#####NLCD Cleaning#####
clean.nlcd <- function(df.nlcd){
year <- unique(df.nlcd$YEAR)
df.nlcd <- df.nlcd[, -1:-2]
df.nlcd <- t(df.nlcd)
df.nlcd <- as.data.frame(df.nlcd)
colnames(df.nlcd) <- as.character(unlist(df.nlcd[1,]))
df.nlcd <- df.nlcd[-1, ]
df.nlcd$site <- rownames(df.nlcd)
df.nlcd$year <- year
df.nlcd$unique_id <- paste0(df.nlcd$site, "_", df.nlcd$year)
rownames(df.nlcd) <- c()
df.nlcd
}

my.list.nlcd <- lapply(nlcd.files, clean.nlcd)

nlcd.master <- do.call(rbind, my.list.nlcd)

#Store and remove the two character variables
nlcd.append <- nlcd.master[, c("site", "unique_id")]
nlcd.master <- nlcd.master[, c(-10, -12)]

#Change all factor variables to integers
nlcd.master <- apply(nlcd.master, c(1,2), FUN = as.integer)
nlcd.master <- as.data.frame(nlcd.master)
#Bind Columns back
nlcd.master <- cbind(nlcd.master, nlcd.append)

#Remove all segments that we generating background 
nlcd.cleaned <- subset(nlcd.master, Undefined <= 5500)

#String Split Routes from segment #s 
#New variable without segment attached
nlcd.cleaned$Route <- str_sub(nlcd.cleaned$site, 1, str_length(nlcd.cleaned$site)-2)

#New Variable with only the segment number 
nlcd.cleaned$Segment <- str_sub(nlcd.cleaned$site, -1)

#Rename columns to be friendlier 
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Undefined"] <- "Background"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Low Intensity Residential/High Intensity Residential/Commercial/Industrial/Transportation/LULC Residential/NLCD/LULC Forested Re"] <- "Urban"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Orchards/Vineyards/Other/LULC Orchards/Vineyards/Other/Pasture/Hay/Row Crops/Small Grains/Fallow"] <- "Ag"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Deciduous Forest/Evergreen Forest/Mixed Forest"] <- "Forest"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Woody Wetlands"] <- "Woody_Wetlands"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Emergent Herbaceous Wetlands"] <- "Emergent_Wetlands"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Open Ocean/Open Water"] <- "Water"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Bare Rock/Sand/Clay/Quarries/Strip Mines/Gravel Pits/Transitional"] <- "Bare"
colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Shrubland/Grasslands/Herbaceous"] <- "Grassland"

#Calculate % cover for wetland in the NLCD segments
#Not in script but for creation of the 1985 averages sites were removed
#Sites don't appear in the wetland.20 DF 
#Code below calculates the 1986 means LULC for all sites to be rbinded 

# nlcd.85 <- nlcd.cleaned[nlcd.cleaned$year == 1980 | nlcd.cleaned$year == 1990,]
# nlcd.85 <- nlcd.85[!nlcd.85$site %in% bad.sites, ]
# site.list <- unique(nlcd.85$site)
# 
# nlcd.means <- data.frame(site = unique(nlcd.85$site), year = 1985, Background = NA, Water= NA, Urban = NA, Bare = NA, 
#                          Forest = NA, Grassland = NA, Ag = NA, Woody_Wetlands = NA, Emergent_Wetlands = NA)
# 
# rownames(nlcd.means) <- nlcd.means[,1]
# nlcd.means <- nlcd.means[, c(-1, -2)]
# 
# for ( i in 1:length(site.list)){
#   nlcd.temp <- nlcd.85[nlcd.85$site == site.list[i], ]
#   lulc.means <- colMeans(nlcd.temp[, 1:9])
#   nlcd.means[rownames(nlcd.means) == site.list[i],] <- lulc.means
# }
# 
# nlcd.means <- tibble::rownames_to_column(nlcd.means, var = "site")
# nlcd.means$year <- "1986"
# nlcd.means$unique_id <- paste0(nlcd.means$site, "_", nlcd.means$year)
# nlcd.site.info <- nlcd.cleaned[, c("site", "Route", "Segment")]
# nlcd.site.info <- nlcd.site.info[!duplicated(nlcd.site.info), ]
# 
# nlcd.means <- left_join(nlcd.means, nlcd.site.info, by = "site")
# nlcd.means <- nlcd.means %>% select("Background", "Water", "Urban", "Bare", "Forest",
#                                     "Grassland", "Ag", "Woody_Wetlands", "Emergent_Wetlands",
#                                     "year", "site", "unique_id", "Route", "Segment")
# 
# 
# write.csv(nlcd.means, file = here::here("Data_Envi/BBS_NLCD_50stop/nlcd1986.csv"))

nlcd.1986 <- read.csv(here::here("Data_Envi/BBS_NLCD_50stop/nlcd1986.csv"))
nlcd.1986 <- nlcd.1986[, -1]

nlcd.cleaned <- rbind(nlcd.cleaned, nlcd.1986)

nlcd.20 <- nlcd.cleaned %>% mutate(total_cover = Urban + Ag + Grassland + Forest + Woody_Wetlands + 
                                     Emergent_Wetlands + Bare + Water) %>% 
  mutate(tot_wetland = Woody_Wetlands + Emergent_Wetlands) %>% 
  mutate(pct_wetland = tot_wetland / total_cover) %>%
  select(unique_id, year, site, Background, Urban, Ag,
         Grassland, Forest, Woody_Wetlands, Emergent_Wetlands,
         Bare, Water, Route, Segment, total_cover, tot_wetland,
         pct_wetland)

#Entire Route level cover  
#Syntax for aggregate, use "." fall all variables
nlcd.agg <- within(nlcd.cleaned, rm("site", "unique_id", "Segment"))
nlcd.agg <- aggregate(data = nlcd.agg, . ~ Route + year, FUN = sum)
nlcd.agg <- within(nlcd.agg, rm("Background"))

# ######Mangrove Cleaning#######
# #Going to need to assign LULC cover classes to all the mangrove data, LA is weird#
# clean.man <- function(df.man){
#   year <- unique(df.man$YEAR)
#   df.man <- df.man[, -1:-2]
#   df.man <- t(df.man)
#   df.man <- as.data.frame(df.man)
#   colnames(df.man) <- as.character(unlist(df.man[1,]))
#   df.man <- df.man[-1, ]
#   df.man$site <- rownames(df.man)
#   df.man$year <- year
#   df.man$unique_id <- paste0(df.man$site, "_", df.man$year)
#   rownames(df.man) <- c()
#   df.man
# }
# 
# my.list.man <- lapply(man.files, clean.man)
# 
# man.master <- do.call(rbind, my.list.nlcd)
# 
# #Store and remove the two character variables
# nlcd.append <- nlcd.master[, c("site", "unique_id")]
# nlcd.master <- nlcd.master[, c(-10, -12)]
# 
# #Change all factor variables to integers
# nlcd.master <- apply(nlcd.master, c(1,2), FUN = as.integer)
# nlcd.master <- as.data.frame(nlcd.master)
# #Bind Columns back
# nlcd.master <- cbind(nlcd.master, nlcd.append)
# 
# #Remove all segments that we generating background 
# nlcd.cleaned <- subset(nlcd.master, Undefined <= 5500)
# 
# #String Split Routes from segment #s 
# #New variable without segment attached
# nlcd.cleaned$Route <- str_sub(nlcd.cleaned$site, 1, str_length(nlcd.cleaned$site)-2)
# 
# #New Variable with only the segment number 
# nlcd.cleaned$Segment <- str_sub(nlcd.cleaned$site, -1)
# 
# #Rename columns to be friendlier 
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Undefined"] <- "Background"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Low Intensity Residential/High Intensity Residential/Commercial/Industrial/Transportation/LULC Residential/NLCD/LULC Forested Re"] <- "Urban"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Orchards/Vineyards/Other/LULC Orchards/Vineyards/Other/Pasture/Hay/Row Crops/Small Grains/Fallow"] <- "Ag"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Deciduous Forest/Evergreen Forest/Mixed Forest"] <- "Forest"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Woody Wetlands"] <- "Woody_Wetlands"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Emergent Herbaceous Wetlands"] <- "Emergent_Wetlands"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Open Ocean/Open Water"] <- "Water"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Bare Rock/Sand/Clay/Quarries/Strip Mines/Gravel Pits/Transitional"] <- "Bare"
# colnames(nlcd.cleaned)[colnames(nlcd.cleaned) == "Shrubland/Grasslands/Herbaceous"] <- "Grassland"
# 
# #Syntax for aggregate, use "." fall all variables
# nlcd.agg <- within(nlcd.cleaned, rm("site", "unique_id", "Segment"))
# nlcd.agg <- aggregate(data = nlcd.agg, . ~ Route + year, FUN = sum)
# nlcd.agg <- within(nlcd.agg, rm("Background"))

#######CCAP Cleaning#########
#Create a UDF for the cleaning 
clean.ccap <- function(df.ccap){
#Pull year out for each df
year <- unique(df.ccap$YEAR)
#Remove two extraneous variables 
df.ccap <- df.ccap[, -1:-2]
#Transform data
df.ccap <- t(df.ccap)
df.ccap <- as.data.frame(df.ccap)
#Make first row colnames
colnames(df.ccap) <- as.character(unlist(df.ccap[1,]))
df.ccap <- df.ccap[-1, ]
df.ccap$site <- rownames(df.ccap)
df.ccap$year <- year
df.ccap$unique_id <- paste0(df.ccap$site, "_", df.ccap$year)
rownames(df.ccap) <- c()
df.ccap
}

#Apply the clean function to all my DFs
my.list <- lapply(ccap.files, clean.ccap)

#Rbind the DFs to each other
ccap.master <- do.call(rbind, my.list)

#Store and remove the two character variables
ccap.append <- ccap.master[, c("site", "unique_id")]
ccap.master <- ccap.master[, c(-10, -12)]

#Change all factor variables to integers
ccap.master <- apply(ccap.master, c(1,2), FUN = as.integer)
#Bind Columns back
ccap.master <- cbind(ccap.master, ccap.append)
ccap.master <- aggregate(data = ccap.master, . ~ unique_id + year + site, FUN = sum)

#Remove all segments that we generating background 
#ccap.cleaned <- subset(ccap.master, Background <= 5500)
ccap.cleaned <- ccap.master

#String Split Routes from segment #s 
#New variable without segment attached
ccap.cleaned$Route <- str_sub(ccap.cleaned$site, 1, str_length(ccap.cleaned$site)-2)

#New Variable with only the segment number 
ccap.cleaned$Segment <- str_sub(ccap.cleaned$site, -1)

#Rename columns to be friendlier 
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Developed, High Intensity/Developed, Medium Intensity/Developed, Low Intensity/Developed, Open Space"] <- "Urban"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Cultivated Crops/Pasture/Hay"] <- "Ag"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Deciduous Forest/Evergreen Forest/Mixed Forest"] <- "Forest"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Palustrine Forested Wetland/Palustrine Scrub/Shrub Wetland/Estuarine Forested Wetland/Estuarine Scrub/Shrub Wetland"] <- "Woody_Wetlands"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Palustrine Emergent Wetland/Estuarine Emergent Wetland"] <- "Emergent_Wetlands"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Open Water/Palustrine Aquatic Bed/Estuarine Aquatic Bed"] <- "Water"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Unconsolidated Shore/Bare Land"] <- "Bare"
colnames(ccap.cleaned)[colnames(ccap.cleaned) == "Grassland/Herbaceous/Scrub/Shrub"] <- "Grassland"

#Get Segments with more than 20% cover
ccap.20 <- ccap.cleaned %>% mutate(total_cover = Urban + Ag + Grassland + Forest + Woody_Wetlands + 
         Emergent_Wetlands + Bare + Water) %>% 
         mutate(tot_wetland = Woody_Wetlands + Emergent_Wetlands) %>% 
         mutate(pct_wetland = tot_wetland / total_cover)

#ccap.20 <- ccap.20[ccap.20$pct_wetland >= 0.2, ]

######Bind CCAP and NLCD together######
#Master DF with CCAP and NLCD data in it 
wetland.20 <- rbind(nlcd.20, ccap.20)

#Pull out maximum % cover segments 
wetland.max <- wetland.20 %>% group_by(site) %>%
  filter(pct_wetland == max(pct_wetland))

#Subset based on whether the max segment is equal to 
#or less than 20%
wetland.max <- wetland.max[wetland.max$pct_wetland >= 0.20,]

#Make a list of segments that fit the criteria above
max.list <- unique(wetland.max$site)

#Subset large DF based on unqiue sites to keep all sites
#that have ever has at least 20% wetland cover

wetland.20 <- wetland.20[wetland.20$site %in% max.list, ]

wetland.20 <- merge(wetland.20, link.to.lulc, by.x = "Route", by.y = "rtename")
wetland.20$rteno.x <- paste0(wetland.20$rteno, "_", wetland.20$Segment)

#Write csv for segment coverage 
#write.csv(wetland.20, file = here("Data_BBS/Generated DFs/CompleteSegments_With_20percent.csv"))



#File written to .csv for convenience
#write.csv(ccap.20, file = here("Data_BBS/Generated DFs/RouteSegs_With_20percent.csv"))

#Syntax for aggregate, use "." fall all variables
ccap.agg <- within(ccap.cleaned, rm("site", "unique_id", "Segment"))
ccap.agg <- aggregate(data = ccap.agg, . ~ Route + year, FUN = sum)
ccap.agg <- within(ccap.agg, rm("Background"))


##Cbind three different LULC datasets
lulc.agg <- rbind(ccap.agg, nlcd.agg)

#Calculate total pixel area, %cover wetland, total wetland
lulc.agg <- lulc.agg %>% mutate(total_cover = Urban + Ag + Grassland + Forest + Woody_Wetlands + 
                                  Emergent_Wetlands + Bare + Water) %>% 
  mutate(tot_wetland = Woody_Wetlands + Emergent_Wetlands) %>% 
  mutate(pct_wetland = tot_wetland / total_cover) 

#Pull out routes that have greater than 20% wetland coverage 
lulc.agg <- subset(lulc.agg, pct_wetland >= 0.2)
lulc.agg$unique_id <- paste0(lulc.agg$Route, "_", lulc.agg$year)
rownames(lulc.agg) <- lulc.agg$unique_id                                       
lulc.agg <- within(lulc.agg, rm("unique_id"))

#Calculate the % of Emergent, Woody Wetlands, and Urban 
lulc.agg <- lulc.agg %>% mutate(pct.ww = Woody_Wetlands / total_cover)
lulc.agg <- lulc.agg %>% mutate(pct.ew = Emergent_Wetlands / total_cover) %>% 
  group_by(Route) %>% mutate(lag.ww = lag(Woody_Wetlands, order_by = year))%>% 
  mutate(pct.change.ww = ((Woody_Wetlands - lag.ww) / lag.ww) * 100) %>% 
  mutate(diff.from.first.ww = (Woody_Wetlands - first(Woody_Wetlands))) %>%
  mutate(ratio.ww = (Woody_Wetlands / Emergent_Wetlands))

lulc.agg <- lulc.agg %>% mutate(pct.ur = Urban / total_cover)

lulc.agg <- lulc.agg %>% mutate(lag.ew = lag(Emergent_Wetlands, order_by = year))%>% 
  mutate(pct.change.ew = ((Emergent_Wetlands - lag.ew) / lag.ew) * 100) %>% 
  mutate(diff.from.first.ew = (Emergent_Wetlands - first(Emergent_Wetlands))) %>%
  mutate(ratio.ew = (Emergent_Wetlands / Woody_Wetlands))

lulc.agg <- lulc.agg %>% mutate(lag.ur = lag(Urban, order_by = year))%>% 
  mutate(pct.change.ur = ((Urban - lag.ur) / lag.ur) * 100) %>% 
  mutate(diff.from.first.ur = (Urban - first(Urban))) %>%
  mutate(ratio.ur = (Urban / total_cover))


##Merge Unique_id with LULC
lulc.agg$Route <- toupper(lulc.agg$Route)
lulc.agg <- merge(lulc.agg, link.to.lulc, by.x = "Route", by.y = "rtename")
lulc.agg$rteno <- as.character(lulc.agg$rteno)
lulc.agg$yr_bin <- 1
lulc.agg$yr_bin[lulc.agg$year >= 1985 & lulc.agg$year <= 1989] <- 2
lulc.agg$yr_bin[lulc.agg$year >= 1990 & lulc.agg$year <= 1994] <- 3
lulc.agg$yr_bin[lulc.agg$year >= 1995 & lulc.agg$year <= 1999] <- 4
lulc.agg$yr_bin[lulc.agg$year >= 2000 & lulc.agg$year <= 2004] <- 5
lulc.agg$yr_bin[lulc.agg$year >= 2005 & lulc.agg$year <= 2009] <- 6
lulc.agg$yr_bin[lulc.agg$year >= 2010 & lulc.agg$year <= 2014] <- 7
lulc.agg$yr_bin[lulc.agg$year >= 2015 & lulc.agg$yearlm <= 2020] <- 8
lulc.agg$unique.id <- paste0(lulc.agg$rteno, "_", lulc.agg$yr_bin)

lulc.agg$stan.ww <- scale(lulc.agg$diff.from.first.ww)
lulc.agg$stan.ew <- scale(lulc.agg$diff.from.first.ew)
lulc.agg$stan.ur <- scale(lulc.agg$diff.from.first.ur)
lulc.agg$stanr.ww <- scale(lulc.agg$ratio.ww)
lulc.agg$stanr.ww <- scale(lulc.agg$ratio.ew)
lulc.agg$stanr.ww <- scale(lulc.agg$ratio.ur)

#Stats Bullshit 
test1 <- separate(test1, unique_id, into = c("State", "Rteno", "Bin"), sep = "_")
test1$unique.id <- paste0(test1$Rteno, "_", test1$Yr_bin)

test.lulc.slice <- test.lulc %>%
  group_by(Rteno) %>%
  slice(2:n())

test.lulc <- merge(test, lulc.agg, by = "unique.id")

mod.lulc <- lmer(beta ~ anomalies + (1|Route), data = test.lulc.slice)
summary(mod.lulc)
eta_sq(mod.lulc)
tab_model(mod.lulc)
r.squaredGLMM(mod.lulc)


test1$pct.ur <- log(test1$pct.change.ur)
plot(test1$beta ~ test1$pct.ur)
lm(test1$beta ~ test1$pct.ur)


##Get residuals for partial regression plots 
ggplot (data = prplot, aes(x = pct.ur, y = y_partial)) + geom_smooth(method = "lm")
devtools::install_github("hohenstein/remef")
library(remef)
beta.t.anom <- remef(mod.lulc, fix = "pct.ur", ran = "all")
beta.t.anom <- as.data.frame(beta.t.anom)
prplot <- test1
prplot <- cbind(prplot, beta.t.anom)
urban.partial <- remef(mod.lulc, fix = "anomalies", ran = "all")
urban.partial <- as.data.frame(temp.partial)
prplot <- cbind(prplot, urban.partial)
pplot.temp <- ggplot (data = prplot, aes(x = anomalies, y = beta.t.anom)) + 
                 geom_smooth(method = "lm") + geom_point()
alpha.partial

pplot.pctur <- ggplot (data = prplot, aes(x = pct.ur, y = urban.partial)) + 
  geom_smooth(method = "lm") + geom_point()
pplot.pctur

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = Yr_Bin(fit$model)[2], y = Yr_Bin(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), x = "Residuals Temperature Anomalies", y = expression(paste("Residuals ", beta, " -Diversity")))}

temperature <- lm(y_partial ~ anomalies, data = prplot)
ggplotRegression(temperature) + labs(x.axis = "Rarefied Alpha")

precip <- lm(precip.partial ~ precip.anomalies, data = prplot)
ggplotRegression(precip)








test <- ccap.agg
hist(test$Sum)
hist(test$pct_wetland)

test$half_wetland <- 0
test$half_wetland[test$pct_wetland>=0.5]<-1

routes <- strsplit(test$site, split = "_")

route_names <- lapply(routes, function(x) unlist(x[[1]]))
route_names <- unlist(lapply(routes, function(x) unlist(x[[1]])))

test$route_names <- unlist(lapply(routes, function(x) unlist(x[[1]])))

test2 <- aggregate(half_wetland ~ route_names, data = test, FUN = "sum")
hist(test2$half_wetland)
sum(test2$half_wetland == 0)

test$half_wetland[test$pct_wetland>=0.20]<-1
test2 <- aggregate(half_wetland ~ route_names, data = test, FUN = "sum")
hist(test2$half_wetland)
sum(test2$half_wetland == 0)

max_wetland <- aggregate(pct_wetland ~ route_names, data = test, FUN = "max")
hist(max_wetland$pct_wetland)
sum(test$half_wetland)
test$half_wetland[test$pct_wetland>=0.25]<-1
sum(test$half_wetland)
shorts <- test[test$Sum <= 2000, ]

sum(test$half_wetland)
test$half_wetland <- 0
test$half_wetland[test$pct_wetland>=0.5]<-1
sum(test$half_wetland)
