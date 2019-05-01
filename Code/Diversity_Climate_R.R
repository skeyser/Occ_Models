#Script with Enviromental and Survey Together

##Clean command 
rm(list = objects(all.names = TRUE))
#Package Loading
install.packages("pacman")
library("pacman")
pacman::p_load("sjPlot", "sjstats", "nlme","reshape2", "effects","vegan", "lmerTest", "tidyverse", "lsr", "cowplot", "here", "purrr", "lme4", "car", "MASS", "MuMIn")
#When I updated R in the console, it wiped some packages like yaml and stringi from the computer
#Code below allowed me to get stringi back
#install.packages("stringi", repos="http://cran.rstudio.com/", dependencies=TRUE)


#Read in all the data for bird obs.
bbs_tx <- read.csv(here("Data_BBS/States_GoM/States/Texas.csv"))
bbs_al <- read.csv(here("Data_BBS/States_GoM/Alabama_New.csv"))
bbs_ms <- read.csv(here("Data_BBS/States_GoM/States/Mississ.csv"))
bbs_fl <- read.csv(here("Data_BBS/States_GoM/States/Florida.csv"))
bbs_la <- read.csv(here("Data_BBS/States_GoM/States/Louisia.csv"))
species <- read.csv(here("Data_BBS/States_GoM/SpeciesList.csv"))

#These were the data for the original script
#routes <- read.csv(here("Data_BBS/States_GoM/bbsrtes_gom_final_refined.csv"))
#routes66 <- read.csv(here("Data_BBS/States_GoM/bbsrts66-98.csv"))
#routes66_xy <- read.csv(here("Data_BBS/States_GoM/bbsrts66-98_xy.csv"))

#Now the script will use the new points, which is the midpoint of each route
#The 1966 and 1999 routes have only 6 routes from 1999 that are not used in 1999
routes.99 <- read.csv(here("Data_BBS/States_GoM/bbsrts_1999_GoM_xy_mid.csv")) 
routes.66 <- read.csv(here("Data_BBS/States_GoM/bbsrts_1966_GoM_xy_mid.csv"))

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

############################################################
############################################################
############***GIS Data Processing in R***##################
############################################################

#Extracting BBS Route Names to pair with recorded start routes
#Creation of the 5 segment lines
#I need the names/rteno of the routes used from GIS
# write.csv(rts_final, file = "Gom_Routes_ID.correct.csv")
# all.routes <- read.csv(here("Data_BBS/States_GoM/routes.csv"))
# Gom.rt.id <- read.csv(here("Data_BBS/States_GoM/Gom_Routes_ID.correct.csv"))
# link.to.lulc <- Gom.rt.id[, c("rteno", "rtename", "")]
# link.to.lulc$rtename <- gsub(" ", "_", link.to.lulc$rtename)
# 
# #Merge two dataframes usin the Rteno
# all.routes$RouteName <- as.character(all.routes$RouteName)
# 
# #Need to generate and merge by rteno
# #Example: Mulitple "moon lake" rts across states
# #Removing by rtename removes all "moon lakes" after the first
# all.routes$rteno <- NA
# 
# #Prepare for the rteno creation
# all.routes$stateindex <- all.routes$statenum * 1000
# all.routes$rteno <- all.routes$stateindex + all.routes$Route
# 
# Gom.rt.id$rtename <- as.character(Gom.rt.id$rtename)
# Gom.rt.id <- Gom.rt.id[, c("rteno", "rtename")]
# colnames(Gom.rt.id)[colnames(Gom.rt.id) == "rtename"] <- "RouteName"
# 
# #Now we can merge by rteno...NOT by RteName
# GoM.rts.true.starts <- merge(all.routes, Gom.rt.id, by = "rteno")
# GoM.rts.true.starts <- GoM.rts.true.starts[!duplicated(GoM.rts.true.starts),]
# GoM.rts.true.starts <- GoM.rts.true.starts[, c(-14, -15)]
# colnames(GoM.rts.true.starts)[colnames(GoM.rts.true.starts) == "RouteName.x"] <- "RteName"
# link.me.up <- Gom.rt.id[, c(3, 5, 6)]

#write.csv(GoM.rts.true.starts, file = "GoM.rt.start.csv")
######################################################
#############***End of GIS Processing***##############
######################################################

#Average Routes by rteno to generate average lat/lon
#for climate data
rts_final <- rts_final %>% group_by(rteno) %>% 
  summarise(centroid.lat = mean(latitude), centroid.lon = mean(longitude))%>%
  ungroup()


#Putting route names with RTENO and Coordinates
sites <- rts[, "rtename"]
sites <- sites[!duplicated(sites),]
rts_final <- cbind(sites, rts_final) 

#Renaming columns 
rts_final$state <- rts_final$rteno
colnames(rts_final)[colnames(rts_final)== "centroid.lon"] <- "longitude"
colnames(rts_final)[colnames(rts_final) == "centroid.lat"] <- "latitude"

#Splitting RTENO and creating state codes
rts_fl <- rts_final %>% filter(str_detect(rteno, "^25"))
rts_fl$state <- 25
bbs_fl <- bbs_fl %>% mutate(proxy = statenum * 1000) %>% mutate(rteno = proxy + Route)
bbs_fl <- merge(bbs_fl, rts_fl, "rteno")
fl <- bbs_fl[, c("rteno", "statenum", "rtename", "Year", "SpeciesTotal", "Aou", "latitude", "longitude")]
colnames(fl) <- tolower(colnames(fl))

#Alabama data slightly different 
rts_al <- rts_final %>% filter(str_detect(rteno, "^2") & !str_detect(rteno, "^25"))
rts_al$state <- 02

#Single Digit State code needs a slightly different approach for joining the rteno 
bbs_al <- mutate(bbs_al, proxy = StateNum * 1000) %>% mutate(rteno = proxy + Route)
bbs_al <- merge(bbs_al, rts_al, "rteno")
bbs_al$BCR <- 27
colnames(bbs_al) <- tolower(colnames(bbs_al))
al <- bbs_al[, c("rteno", "statenum", "rtename", "year", "speciestotal", "aou", "latitude", "longitude")]

rts_tx <- rts_final %>% filter(str_detect(rteno, "^83"))
rts_tx$state <- 83
bbs_tx <- mutate(bbs_tx, proxy = statenum * 1000) %>% mutate(rteno = proxy + Route)
bbs_tx <- merge(bbs_tx, rts_tx, "rteno")
colnames(bbs_tx) <- tolower(colnames(bbs_tx))
tx <- bbs_tx[, c("rteno", "statenum", "rtename", "year", "speciestotal", "aou", "latitude", "longitude")]

rts_la <- rts_final %>% filter(str_detect(rteno, "^42"))
rts_la$state <- 42
bbs_la <- mutate(bbs_la, proxy = statenum * 1000) %>% mutate(rteno = proxy + Route)
bbs_la <- merge(bbs_la, rts_la, "rteno")
colnames(bbs_la) <- tolower(colnames(bbs_la))
la <- bbs_la[, c("rteno", "statenum", "rtename", "year", "speciestotal", "aou", "latitude", "longitude")]

rts_ms <-rts_final %>% filter(str_detect(rteno, "^51"))
rts_ms$state <- 51
bbs_ms <- mutate(bbs_ms, proxy = statenum * 1000) %>% mutate(rteno = proxy + Route)
bbs_ms <- merge(bbs_ms, rts_ms, "rteno")
colnames(bbs_ms) <- tolower(colnames(bbs_ms)) 
ms <- bbs_ms[, c("rteno", "statenum", "rtename", "year", "speciestotal", "aou", "latitude", "longitude")]

#Join all states together 
bbs_total <- rbind(tx, al, ms, fl, la)


#Creation of species scientific name 
species$Linnean <- paste0(species$Genus, " ", species$Species)

#Merge species names with AOU
colnames(species)[2] <- "Aou"
colnames(species) <- tolower(colnames(species))
species_sub <- species[, c(-1, -4)]

#Merging species name data with observation data 
bbs_total <- merge(bbs_total, species_sub, by = "aou")

#Adding state abbreviation distinction
bbs_total$state <- 1
bbs_total$state[bbs_total$statenum == 2] <- "AL"
bbs_total$state[bbs_total$statenum == 83] <- "TX"
bbs_total$state[bbs_total$statenum == 51] <- "MS"
bbs_total$state[bbs_total$statenum == 25] <- "FL"
bbs_total$state[bbs_total$statenum == 42] <- "LA"

bbs_total$unique_id <- paste0(bbs_total$state, "_", bbs_total$rteno, "_", bbs_total$year)
bbs_total$abbrev <- paste0(bbs_total$state, "_", bbs_total$rteno)

#Create a spp by site matrix
#Keep only the columns we need
bbs_simple <- bbs_total[,c("linnean", "unique_id", "speciestotal")]
bbs_plot <- bbs_total[, c("unique_id", "statenum")]
bbs_plot <- as.data.frame(bbs_plot)

#Cast Data
bbs_cast <- dcast(data = bbs_simple, formula = unique_id ~ linnean, fun.aggregate = sum)

length(unique(bbs_total$linnean))
length(unique(bbs_total$unique_id))

rownames(bbs_cast) <- bbs_cast$unique_id
bbs_cast <- bbs_cast[,-1]


Total_Abundance <- rowSums(bbs_cast)
Total_Abundance <- cbind(rownames(bbs_cast), Total_Abundance)

#convert spp by site matric to PA matrix
bbs_pa <- bbs_cast
bbs_pa[bbs_pa > 0] <- 1

#Species richness for each route_year
Total_SR <- as.data.frame(rowSums(bbs_pa))

#Create DF same as bbs_total
bbs_master <- bbs_total

#Bind the richness and abundance
bbs_total <- cbind(Total_Abundance, Total_SR)
colnames(bbs_total) <- c("unique_id", "Total_Abundance", "Total_SR")
bbs_total$Total_Abundance <- as.numeric(bbs_total$Total_Abundance)
bbs_total$unique_id <- as.character(bbs_total$unique_id)

#plot(bbs_total$Total_Abundance, bbs_total$Total_SR)

#Create a RouteName Column for bbs_total
bbs.rt.info <- bbs_master[, c("unique_id", "abbrev", "rteno", "rtename")]
bbs_total <- merge(bbs_total, bbs.rt.info, by = "unique_id")

#Calculate rarefied species richness
min(bbs_total$Total_Abundance)

spp.accum <- specaccum(bbs_cast[1:10,])
plot(spp.accum)

#loop through 10 random sites
rand_site <- sample(1:nrow(bbs_cast), 12, replace = FALSE )

par(mfrow = c(3,4))

for (i in 1:length(rand_site)){
  site.ind <- rand_site[i]
  rand.comm <- bbs_cast[site.ind,]
  sub.sample <- seq(from= 10, to = 1500, by = 10)
  spp.accum.rand <- rarefy(rand.comm, sample = sub.sample )
  plot(x = sub.sample, y = spp.accum.rand)
}

length(unique(bbs_total$rteno)) #129 routes

#Run the rarefaction to 200N, good trade-off in terms of n sampled and not loosing surveys
rarefied.alpha <- rarefy(bbs_cast, sample = 400)
Total_N <- rowSums(bbs_cast)
Total_400 <- Total_N >= 400

sum(Total_400)
sum(!Total_400)

rarefied.alpha[!Total_400] <- NA

#Large correlation between 1000 and 400 subsample 
#not necessary to do rarefy to 1000, lose extra
# rarefied.alpha.1000 <- rarefy(bbs_cast, sample = 1000)
# Total_1000 <- Total_N >=1000
# 
# rarefied.alpha.1000[!Total_1000] <- NA
# 
# plot(rarefied.alpha, rarefied.alpha.1000)
# cor.test(rarefied.alpha, rarefied.alpha.1000)


#check names are the same so we can just cbind them
bbs_total <- bbs_total[!duplicated(bbs_total),]
sum(rownames(bbs_cast) == bbs_total$unique_id)
bbs_total <- cbind(bbs_total, rarefied.alpha)

#******Complete cases used******
bbs_total <- bbs_total[complete.cases(bbs_total), ]


#Add in year as another column
years <- strsplit(bbs_total$unique_id, split = "_")
bbs_total$Year <-as.numeric(sapply(years, function(x) x[[3]]))
bbs_total$Abbrev <- sapply(years, function(x) x[[1]])

#Correlation of ~0.83 for Total SR and rarefied alpha 
plot(bbs_total$Total_SR, bbs_total$rarefied.alpha)
cor.test(bbs_total$Total_SR, bbs_total$rarefied.alpha)

#Pull out lat longs and year
site_data <- unique(bbs_master[,c("unique_id", "year", "latitude", "longitude", "abbrev")])

site_data_merge <- merge(site_data, bbs_total, by = "unique_id")

# plot(site_data_merge$Total_SR ~ site_data_merge$year)
# plot(site_data_merge$rarefied.alpha ~ site_data_merge$year)
# 
# plot(site_data_merge$Total_SR ~ site_data_merge$latitude)
# plot(site_data_merge$rarefied.alpha ~ site_data_merge$latitude)
# 
# plot(site_means$rarefied.alpha ~ site_means$latitude)
# text(site_means$latitude , site_means$Total_Abundance, site_means$abbrev, cex = 0.8)
# abline(lm(site_means$rarefied.alpha ~ site_means$latitude))
# cor(site_means$rarefied.alpha, site_means$latitude)
# 
# plot(site_means$Total_SR ~ site_means$latitude, pch = "")
# text(site_means$latitude , site_means$Total_SR, site_means$abbrev, cex = 0.8)
# abline(lm(site_means$Total_SR ~ site_means$latitude))
# cor(site_means$Total_SR, site_means$latitude)

#Clean up some remaining names from merge 
#Could usesubset indexes, but want to make sure no hang-ups
colnames(site_data_merge)[colnames(site_data_merge) == "abbrev.y"] <- "abbrev"
colnames(site_data_merge)[colnames(site_data_merge) == "Abbrev"] <- "State"
drop.it <- "abbrev.x"
site_data_merge <- site_data_merge[, !names(site_data_merge) %in% drop.it, drop = F]

#Calculate mean SR, Abundance, and rarefied SR 
site_means <- aggregate(data = site_data_merge, cbind(Total_Abundance, Total_SR, rarefied.alpha) ~ abbrev + latitude, FUN = mean)



n.sites <- length(unique(site_data_merge$abbrev))
site.list <- as.character(unique(site_data_merge$abbrev))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp,]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$rarefied.alpha ~ site_data_merge_temp$year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites[i,2] <- slope.temp
  }
}

hist(slopes_sites$slope, breaks = 50)

slopes_sites$abbrev <- as.character(slopes_sites$site.list)

site_means2 <- merge(site_means, slopes_sites, by = "abbrev")

plot(site_means2$slope ~ site_means2$latitude)
abline(h= 0)

mean(site_means2$slope, na.rm = T)
median(site_means2$slope, na.rm = T)


#mean temporal beta diversity by site
mean.betas <- data.frame(site.list, mean.beta = NA)

#Extract sites and years from 1980 to current 
site_list_80 <- site_data$unique_id[site_data$year >= 1980] 
bbs_cast_80 <- bbs_cast[rownames(bbs_cast) %in% site_list_80,]

#Create matrix for storing Betas 
beta.matrix <- data.frame(unique_id = site_list_80, beta = NA)

#Fixing data frame to have abbrev
bbs_total$unique_id <- paste0(bbs_total$abbrev, "_", bbs_total$year)


#Loop for calculating Beta Diversity from baseline 5 years
#If I average the Beta diversity after will this generate 
#statistical biases???
#Should we loop it from year bin to year bin?
for (n in 1:n.sites){
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(bbs_total[which(bbs_total$abbrev == site.temp),"Year"])
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
site_data_merge <- merge(site_data_merge, beta.matrix, by = "unique_id")

#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, beta.slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp,]
  site_data_merge_temp <- site_data_merge_temp[complete.cases(site_data_merge_temp),]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$beta ~ site_data_merge_temp$year)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
  }
}

hist(slopes_sites_beta$beta.slope)
t.test(slopes_sites_beta$beta.slope, mu = 0)

site_data_merge <- merge(site_data_merge, bbs_plot, by = "unique_id")
site_data_merge$state <- NA
site_data_merge$state[site_data_merge$statenum == 25] <- "FL"
site_data_merge$state[site_data_merge$statenum == 83] <- "TX"
site_data_merge$state[site_data_merge$statenum == 42] <- "LA"
site_data_merge$state[site_data_merge$statenum == 51] <- "MS"
site_data_merge$state[site_data_merge$statenum == 2] <- "AL"

#Plots for Beta Diversity
test <- merge(site_data_merge, slopes_sites_beta, by = )

gghist_beta <- ggplot(data = slopes_sites_beta, aes(slopes_sites_beta$beta.slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", beta, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()

#site_data_merge$Year <- site_data_merge$count_yr + 1900


plot_beta <- (ggplot(site_data_merge, aes(x = year, y = beta, group = abbrev, 
                                          colour = factor(state, 
                                                          labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = lm, se = FALSE, aes(x = year, y = beta, group = abbrev, 
                                                         colour = factor(state, 
                                                                         labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                xlab("Year") +
                ylab(expression(paste(beta, "-diversity"))) +
                labs(colour = "States") +
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

bbs_total$unique_id <- paste0(bbs_total$unique_id, "_", bbs_total$Year)
bbs_total$site <- bbs_total$unique_id

bbs_total <- bbs_total %>% separate(site, c("StateName", "Rteno", "Yr"), sep = "_") %>% 
  unite(site, StateName, Rteno, sep = "_") 

Mean_alpha_bin <- aggregate(data = bbs_total, rarefied.alpha ~ Yr_bin + abbrev, FUN = "mean")
Mean_alpha_bin <- separate(Mean_alpha_bin, site, c("state", "rteno", "yr"), sep = "_")

plot(Mean_alpha_bin$rarefied.alpha ~ Mean_alpha_bin$Yr_bin)
lm_alpha <- lm(rarefied.alpha~ Yr_bin + abbrev, data = Mean_alpha_bin)
summary(lm_alpha)
anova(lm_alpha)
library(lme4)
lmer_alpha <- lmer(rarefied.alpha~ Yr_bin + (1|abbrev), data = Mean_alpha_bin)
summary(lmer_alpha)
anova(lmer_alpha)
dev.off()

#Slopes for mean 1980 rarefied alphas 
n.sites <- length(unique(site_data_merge$abbrev))
site.list <- as.character(unique(site_data_merge$abbrev))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp & site_data_merge$year >= 1980, ]
  site_data_merge_temp <- site_data_merge_temp[!duplicated(site_data_merge_temp$unique_id),]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$rarefied.alpha ~ site_data_merge_temp$year)
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


plot_alpha <- (ggplot(site_data_merge, aes(x = year, y = rarefied.alpha, group = abbrev, 
                                           colour = factor(State, 
                                                           labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = glm,se = FALSE, aes(x = year, y = rarefied.alpha, group = abbrev, 
                                                          colour = factor(State, 
                                                                          labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 xlab("Year") +
                 ylab(expression(paste("Rarefied ", alpha, "-diversity"))) +
                 labs(colour = "States") +
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
my_mean.alpha = aggregate(Mean_alpha_bin$rarefied.alpha , by=list(Mean_alpha_bin$Yr_bin) , mean); colnames(my_mean.alpha)=c("Yr_Bin" , "mean")
se <- function(x) sqrt(var(x)/length(x))
my_sd_alpha=aggregate(Mean_alpha_bin$rarefied.alpha , by=list(Mean_alpha_bin$Yr_bin) , se); colnames(my_sd_alpha)=c("Yr_Bin" , "sd")
my_info.rarefied=merge(my_mean.alpha , my_sd_alpha , by.x=1 , by.y=1)


box4 <- ggplot(my_info.rarefied) + 
  # geom_point(aes(x = Yr_Bin, y = CMRL) , colour=rgb(0.8,0.7,0.1,0.4) , size=5) + 
  geom_point(data = my_info.rarefied, aes(x=Yr_Bin , y = mean) , colour = rgb(0.6,0.5,0.4,0.7) , size = 8) +
  geom_errorbar(data = my_info.rarefied, aes(x = Yr_Bin, y = sd, ymin = mean - sd, ymax = mean + sd), colour = rgb(0.4,0.8,0.2,0.4) , width = 0.7 , size=1.5) +
  scale_y_continuous(name = "Rarefied Alpha-diversity", breaks = seq(35, 50, 1)) +
  scale_x_continuous(name = "Year Bin", breaks = seq(1, 8, 1), 
                     labels = c("1980 - 1984", "1985 - 1989", "1990 - 1994", "1995 - 1999", 
                                "2000 - 2004", "2005 - 2009", "2010 - 2014", "2015 - 2019")) +
  coord_cartesian(ylim = c(35, 50)) 

box4
######################################################




#Environmental Data Loading and Cleaning 
#PRISM Data Script 
prism.1 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_198001_199412.csv"))
prism.2 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_199501_200912.csv"))
prism.3 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1966_GoM_mid_201001_201712.csv"))
prism.4 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_198001_199412.csv"))
prism.5 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_199501_200912.csv"))
prism.6 <- read.csv(here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_201001_201712.csv"))
prism.1.monthly <- read.csv(here("Data_Envi/PRISM Data/PRISM_1966_GoM_mids_monthly_normals.csv")) 
prism.2.monthly <- read.csv(here("Data_Envi/PRISM Data/PRISM_1999_GoM_mid_monthly_normals.csv"))

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
betas <- beta.matrix %>% separate(unique_id, c("State", "route", "year"), sep = "_") %>% unite(site, c("State", "route"), sep = "_")
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
bbs_site_temp <- bbs_total[, c("rteno", "rtename", "abbrev")]
bbs_site_temp$rtename <- as.character(bbs_site_temp$rtename)
bbs_site_temp <- bbs_site_temp[!duplicated(bbs_site_temp$rteno), ]
colnames(bbs_site_temp)[colnames(bbs_site_temp) == "rtename"] <- "Site"

#Create a site.link dataframe 
site.link <- bbs_site_temp[, c("Site", "abbrev", "rteno")]
site.link <- site.link[!duplicated(site.link),]
site.link$site <- site.link$abbrev
min.betas <- merge(min.betas, site.link, by = "site")

colnames(min.betas)[colnames(min.betas) == "Yr_bin"] <- "min.yr.bin"

#Merge the climate averages with the site data for ref
avgs <- merge(avgs, site.link, by = "Site")
avgs$unique_id <- paste0(avgs$rteno, "_", avgs$yr_bin)
avgs <- merge(avgs, min.betas, by = "Site")
avgs <- avgs[avgs$yr_bin >= avgs$min.yr.bin,]


#Metric Conversion
avgs <- avgs %>% mutate(tmean.c = (tmean - 32) / 1.8, 
                        tmax.c = (avg.max - 32) / 1.8,
                        tmin.c = (avg.min - 32) / 1.8,
                        precipitation.c = (precipitation * 2.54))

#Calculated anomalies relative to the base year
avgs <- avgs %>% group_by(Site, Month) %>% 
  arrange(yr_bin, .by_group = T) %>% 
  mutate(anomalies = tmean.c - first(tmean.c)) %>% 
  mutate(max.anomalies = tmax.c - first(tmax.c)) %>% 
  mutate(min.anomalies = tmin.c - first(tmin.c)) %>%
  mutate(precip.anomalies = precipitation.c - first(precipitation.c))

#Same as above for anomalies, need to merge the data for betas & anomalies 
#using the join field of State_Rteno
temp.agg <- avgs[, c(16, 3, 24, 25, 26, 27)]

#Aggregate all the anomalies by year_bin for each site
temp.agg <- aggregate(. ~ site + yr_bin, temp.agg, FUN = "mean")
temp.avgs <- temp.agg
betas.mod$unique_id <- paste0(betas.mod$site, "_", betas.mod$Yr_bin)
temp.avgs$unique_id <- paste0(temp.avgs$site, "_", temp.avgs$yr_bin)

#Put betas with temperature data 
lmer.df <- merge(temp.avgs, betas.mod, by = "unique_id")
test <- lmer.df


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
t.test(test1$p.change.alpha, mu = 0)
r.squaredGLMM(mod.alpha)
plot(test1$precip.anomalies, test1$p.change.alpha)
abline(lm(test1$p.change.alpha ~ test1$precip.anomalies))

##########################################################################
###################Code Cleaned Up To Here 4/24/2019######################
##########################################################################


a.test <- test1[, c("Site", "Yr_bin", "p.change.alpha", "rarefied.alpha", "anomalies", "precip.anomalies", "min.anomalies", "max.anomalies")]
a.test <- a.test[a.test$p.change.alpha != 0, ]
s <- scale(a.test$p.change.alpha)
GHQ <- glmer(p.change.alpha ~ anomalies + precip.anomalies + (1 | site), data = a.test, family = poisson(link = "logit"), nAGQ = 25) # Set nAGQ to # of desired iterations
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
ggplot (data = prplot, aes(x = anomalies, y = precip.partial)) + geom_smooth(method = "lm")

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

##If the model include rarefied alpha need to fix for all three effects
mod.fix <- lmer(beta ~ anomalies + rarefied.alpha + precip.anomalies + (1|Site), data = test, REML = F)
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
res.beta <- residuals(lm(beta ~ precip.anomalies + rarefied.alpha, data = test))
res.temp <- residuals(lm(anomalies ~ precip.anomalies + rarefied.alpha, data = test))
res <- cbind(as.data.frame(res.beta), as.data.frame(res.temp))
res.beta.p <- residuals(lm(beta ~ anomalies + rarefied.alpha, data = test))
res.precip <- residuals(lm(precip.anomalies ~ anomalies + rarefied.alpha, data = test))
res2 <- cbind(as.data.frame(res.beta.p), as.data.frame(res.precip))
res.beta.a <- residuals(lm(beta ~ anomalies + precip.anomalies, data = test))
res.alpha <- residuals(lm(rarefied.alpha ~ anomalies + precip.anomalies, data = test))
res3 <- cbind(as.data.frame(res.beta.a), as.data.frame(res.alpha))
plot(res.temp, res.beta)
k <- lm(res$res.beta ~ res$res.temp)
p <- lm(res2$res.beta.p ~ res2$res.precip)
g <- lm(res3$res.beta.a ~ res3$res.alpha)
ggplotRegression(k)
ggplotRegression(p)
ggplotRegression(g)
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
bcr <- read.csv(here("Data_BBS/States_GoM/bbsrts_bcr.csv"))
bcr <- bcr[, c("RTENO", "BCR")]
bcr$RTENO <- as.numeric(bcr$RTENO)
bcr$BCR <- as.numeric(bcr$BCR)
test$RTENO <- as.numeric(test$RTENO)
bcr <- bcr[!duplicated(bcr$RTENO), ]
test <- merge(test, bcr, by = "RTENO")


#Table of output 
tab_model(mod.fix, p.val = "wald", dv.labels = "Beta Diversity", pred.labels = c("Intercept", "Temperature Anomalies", "Precipitation Anomalies", "Alpha Diversity"), title = "Linear Mixed-Effects Model Output", digits = 4, digits.p = 4, p.style = "asterisk", file = here("Figures and Tables"))

#Getting Dataset to make a map for route locations 
climate.map <- separate(slopes, site.list, c("State", "rteno"), sep = "_")
climate.map <- merge(climate.mao, rts, by = "rteno")

write.csv(climate.map, file = here("climate.map.csv"))



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
