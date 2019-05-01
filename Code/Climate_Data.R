#PRISM Data Script 
rm(list = objects(all.names = TRUE))
pacman::p_load("tidyverse", "reshape2", "lme4", "here")
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

#Need to link BBS data so that each survey calculates a "running anomaly" each survey starts from a different year
#Needs to be done from the survey data 
#Creation of a usable beta matrix
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

min.betas <- aggregate(Yr_bin ~ site, betas, FUN = "min")
site.link <- bbs_site_temp[, c("Site", "state", "rteno")]
site.link <- site.link[!duplicated(site.link),]
site.link$site <- paste0(site.link$state, "_", site.link$rteno)
min.betas <- merge(min.betas, site.link, by = "site")
colnames(min.betas)[colnames(min.betas) == "Yr_bin"] <- "min.yr.bin"
avgs <- merge(avgs, site.link, by = "Site")
avgs$unique_id <- paste0(avgs$site, "_", avgs$yr_bin)
avgs <- merge(avgs, min.betas, by = "site")
avgs <- avgs[avgs$yr_bin >= avgs$min.yr.bin,]
avgs <- avgs[, c(-9, -10, -11, -14)]

#Back-up file with info for joining as a link file 
#write.csv(site.link, file = "C:/Users/Spencer/Box Sync/MSI_Research/R/R_proj_Keyser_Masters/Data_BBS/States_GoM/site_link.csv")

#Calculated anomalies
avgs <- avgs %>% group_by(site, Month) %>% 
  arrange(yr_bin, .by_group = T) %>% 
  mutate(anomalies = tmean - first(tmean)) %>% 
  mutate(max.anomalies = avg.max - first(avg.max)) %>% 
  mutate(min.anomalies = avg.min - first(avg.min)) %>%
  mutate(precip.anomalies = precipitation - first(precipitation))

  
  
  
#For running regression on temps  
yr.avgs <- aggregate(anomalies ~ site + yr_bin, avgs ,FUN = "mean")
yr.avgs.max <- aggregate(max.anomalies ~ site + yr_bin, avgs ,FUN = "mean")
yr.avgs.min <- aggregate(min.anomalies ~ site + yr_bin, avgs ,FUN = "mean")
yr.avgs.pcp <- aggregate(precip.anomalies ~ site + yr_bin, avgs, FUN = "mean")

yr.avgs <- do.call("cbind", list(yr.avgs, yr.avgs.max, yr.avgs.min, yr.avgs.pcp))
yr.avgs <- yr.avgs[, c(1, 2, 3, 6, 9, 12)]
#Loop Through to run the OLS Regression 
n.sites <- length(unique(yr.avgs$site))
site.list <- as.character(unique(yr.avgs$site))

slopes <- data.frame(site.list, list(slope = "NA", slope.max = "NA", slope.min = "NA", slope.precip = "NA"))
slopes$site.list <- as.character(slopes$site.list)
slopes$slope <- as.numeric(slopes$slope)
slopes$slope.max <- as.numeric(slopes$slope.max)
slopes$slope.min <- as.numeric(slopes$slope.min)
slopes$slope.precip <- as.numeric(slopes$slope.precip)

for (i in 1:n.sites){
  site.tmp <- site.list[i]
  avgs_tmp <- yr.avgs[yr.avgs$site == site.tmp, ]
  if(nrow(avgs_tmp) > 1){
    lm.temp <- (lm(avgs_tmp$anomalies ~ avgs_tmp$yr_bin))
    lm.temp.max <- (lm(avgs_tmp$max.anomalies ~ avgs_tmp$yr_bin))
    lm.temp.min <- lm(avgs_tmp$min.anomalies ~ avgs_tmp$yr_bin)
    lm.precip <- lm(avgs_tmp$precip.anomalies ~ avgs_tmp$yr_bin)
    slope.tmp <- summary(lm.temp)$coefficients[2,1]
    slope.tmp.max <- summary(lm.temp.max)$coefficients[2,1]
    slope.tmp.min <- summary(lm.temp.min)$coefficients[2,1]
    slope.tmp.precip <- summary(lm.precip)$coefficients[2,1]
    slopes[i,2] <- slope.tmp
    slopes[i,3] <- slope.tmp.max
    slopes[i,4] <- slope.tmp.min
    slopes[i,5] <- slope.tmp.precip
  }
}
hist(avgs$anomalies)
summary(lm.temp)$coeffecients[2,1]
names(lm.temp)
slopes$slope <- as.numeric(slopes$slope)
slopes$site.list <- as.character(slopes$site.list)

temp.avgs <- separate(temp.avgs, unique_id, c("state", "Route", "year_bin"), sep = "_")
temp.avgs$year <- 1
temp.avgs$year[temp.avgs$yr_bin == 1] <- 1980
temp.avgs$year[temp.avgs$yr_bin == 2] <- 1985
temp.avgs$year[temp.avgs$yr_bin == 3] <- 1990
temp.avgs$year[temp.avgs$yr_bin == 4] <- 1995
temp.avgs$year[temp.avgs$yr_bin == 5] <- 2000
temp.avgs$year[temp.avgs$yr_bin == 6] <- 2005
temp.avgs$year[temp.avgs$yr_bin == 7] <- 2010
temp.avgs$year[temp.avgs$yr_bin == 8] <- 2015
temp.avgs <- temp.avgs[temp.avgs$yr_bin > 1,]

gghist_temp <- ggplot(data = slopes, aes(slopes$slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  labs(title = "") +
  labs(x = expression(paste("Slopes of temperature anomalies")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()


plot_temp <- (ggplot(temp.avgs, aes(x = year, y = anomalies, group = Route, 
                                          colour = factor(state, 
                                                          labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = lm, se = FALSE, aes(x = year, y = anomalies, group = Route, 
                                                         colour = factor(state, 
                                                                         labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                xlab("Year") +
                ylab(expression(paste("Temperature Anomalies " (degree*C)))) +
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
  draw_plot(plot_temp + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_temp, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 


#This is before I integrated precip in the beginning 
# #Precipitation
# avgs.pcp <- aggregate(precip ~ Site +  Month + yr_bin, prism.simple, FUN = "mean")
# 
# avgs.pcp$site.mo.yr <- paste0(avgs.pcp$Site, "_", avgs.pcp$Month, "_", avgs.pcp$yr_bin)
# 
# avgs.pcp <- avgs.pcp %>% group_by(Site, Month) %>% 
#   arrange(yr_bin, .by_group = T) %>% 
#   mutate(anomalies.pcp = precip - first(precip))
# 
# 
# yr.avgs.pcp <- aggregate(anomalies.pcp ~ Site + yr_bin, avgs.pcp ,FUN = "mean")
# 
# n.sites.pcp <- length(unique(yr.avgs.pcp$Site))
# site.list.pcp <- as.character(unique(yr.avgs.pcp$Site))
# 
# slopes.pcp <- data.frame(site.list, slope.precip = "NA")
# slopes.pcp$site.list <- as.character(slopes.pcp$site.list)
# slopes.pcp$slope.precip <- as.numeric(slopes.pcp$slope.precip)
# 
# for (i in 1:n.sites.pcp){
#   site.tmp.pcp <- site.list.pcp[i]
#   avgs_tmp_pcp <- yr.avgs.pcp[yr.avgs.pcp$Site == site.tmp.pcp, ]
#   if(nrow(avgs_tmp_pcp) > 1){
#   lm.temp.pcp <- (lm(avgs_tmp_pcp$anomalies ~ avgs_tmp_pcp$yr_bin))
#     slope.tmp.pcp <- summary(lm.temp.pcp)$coefficients[2,1]
#     slopes.pcp[i,2] <- slope.tmp.pcp
#   }
# }
# 
# #Get the two datasets together 
# total.slopes <- cbind(slopes, slopes.pcp)
# total.slopes <- total.slopes[, -5]