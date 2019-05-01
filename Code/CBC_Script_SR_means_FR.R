#Package Loading
install.packages("pacman")
library("pacman")
pacman::p_load(reshape2, vegan, tidyverse, cowplot, dplyr, MuMIn)
##Clean command 
rm(list = objects(all.names = TRUE))

# #Set working directory for work computer
# setwd("C:/Users/Spencer/Box Sync/MSI_Research/Data/CBC Data")
# 
# #Set working directory for personal computer
# setwd("C:/Users/Spencer Keyser/Box Sync/Keyser Thesis/Data/CBC Data")


cbc_df <- read.csv(file = here("Data_CBC/Surveys/CBC_Species_Sites_Total_Refined_NoCarib_Cleaned_CSV.csv"), header = TRUE)

#Creation of the Unique IDs
cbc_df$unique_id <- paste0(cbc_df$Abbrev, "_", cbc_df$Count_yr)

#make names lower case 
names(cbc_df) <- tolower(names(cbc_df))


#Create a spp by site matrix
#Keep only the columns we need
cbc_simple <- cbc_df[,c("unique_id", "sci_name", "how_many")]
cbc_plot <- cbc_df[, c("unique_id", "subnational_code")]
cbc_plot <- as.data.frame(cbc_plot)
#Cast Data
cbc_cast <- dcast(data = cbc_simple, formula = unique_id ~ sci_name, fun.aggregate = sum)
#as.matrix(cbc_cast)

length(unique(cbc_df$sci_name))
length(unique(cbc_df$unique_id))

rownames(cbc_cast) <- cbc_cast$unique_id
cbc_cast <- cbc_cast[,-1]


Total_Abundance <- rowSums(cbc_cast)
Total_Abundance <- cbind(rownames(cbc_cast), Total_Abundance)

#convert spp by site matric to PA matrix
cbc_pa <- cbc_cast
cbc_pa[cbc_pa>0] <- 1


Total_SR <- as.data.frame(rowSums(cbc_pa))
Summary_stats <- cbind(Total_Abundance, Total_SR)
colnames(Summary_stats) <- c("unique_id", "Total_Abundance", "Total_SR")
Summary_stats$Total_Abundance <- as.numeric(Summary_stats$Total_Abundance)
Summary_stats$unique_id <- as.character(Summary_stats$unique_id)


plot(Summary_stats$Total_Abundance, Summary_stats$Total_SR)

#Calculate rarefied species richness
min(Summary_stats$Total_Abundance)

spp.accum <- specaccum(cbc_cast[1:10,])
plot(spp.accum)

#loop through 10 random sites
rand_site <- sample(1:nrow(cbc_cast), 12, replace = FALSE )

par(mfrow = c(3,4))

for (i in 1:length(rand_site)){
  site.ind <- rand_site[i]
  rand.comm <- cbc_cast[site.ind,]
  sub.sample <- seq(from= 10, to = 2000, by = 10)
  spp.accum.rand <- rarefy(rand.comm, sample = sub.sample )
  plot(x = sub.sample, y = spp.accum.rand)
}
  
length(unique(cbc_df$abbrev)) #77 sites

#Run the rarefaction to 200N, good trade-off in terms of n sampled and not loosing surveys

rarefied.alpha <- rarefy(cbc_cast, sample = 200)

#check names are the same so we can just cbind them
sum(rownames(cbc_cast) == Summary_stats$unique_id)
Summary_stats <- cbind(Summary_stats, rarefied.alpha)
#Add in year as another column
years <- strsplit(Summary_stats$unique_id, split = "_")
Summary_stats$Year <-as.numeric(sapply(years, function(x) x[[2]]))+1900
Summary_stats$Abbrev <- sapply(years, function(x) x[[1]])

plot(Summary_stats$Total_SR, Summary_stats$rarefied.alpha)
cor(Summary_stats$Total_SR, Summary_stats$rarefied.alpha)


#Import effort data
cbc_effort_1 <- read.csv(file = "Keyser-CBC_Effort_Report_SQL_updated-1.csv", header = TRUE)
cbc_effort_2 <- read.csv(file = "Keyser-CBC_Effort_Report_SQL_updated-2.csv", header = TRUE)
cbc_effort_1$unique_id <- paste0(cbc_effort_1$Abbrev, "_", cbc_effort_1$Count_yr)
cbc_effort_2$unique_id <- paste0(cbc_effort_2$Abbrev, "_", cbc_effort_2$Count_yr)

cbc_effort <- merge(cbc_effort_1, cbc_effort_2, by = "unique_id")
#names(cbc_effort) <- tolower(cbc_effort)

Summary_stats2 <- merge(Summary_stats, cbc_effort, by = "unique_id", all = FALSE)
plot(y = Summary_stats2$rarefied.alpha, x = as.numeric(Summary_stats2$Distance))
cor(Summary_stats2$rarefied.alpha, as.numeric(Summary_stats2$Distance))

plot(y = Summary_stats2$rarefied.alpha, x = as.numeric(Summary_stats2$Field_counters))
cor(Summary_stats2$rarefied.alpha, as.numeric(Summary_stats2$Field_counters))

#Pull out lat longs and year
site_data <- unique(cbc_df[,c("unique_id", "count_yr", "latitude", "longitude", "abbrev")])

site_data_merge <- merge(site_data, Summary_stats, by = "unique_id")

plot(site_data_merge$Total_SR ~ site_data_merge$count_yr)
plot(site_data_merge$rarefied.alpha ~ site_data_merge$count_yr)

plot(site_data_merge$Total_SR ~ site_data_merge$latitude)
plot(site_data_merge$rarefied.alpha ~ site_data_merge$latitude)

site_means <- aggregate(data = site_data_merge, cbind(Total_Abundance, Total_SR, rarefied.alpha) ~ abbrev + latitude, FUN = mean)
plot(site_means$rarefied.alpha ~ site_means$latitude)
text(site_means$latitude , site_means$Total_Abundance, site_means$abbrev, cex = 0.8)
abline(lm(site_means$rarefied.alpha ~ site_means$latitude))
cor(site_means$rarefied.alpha, site_means$latitude)

plot(site_means$Total_SR ~ site_means$latitude, pch = "")
text(site_means$latitude , site_means$Total_SR, site_means$abbrev, cex = 0.8)
abline(lm(site_means$Total_SR ~ site_means$latitude))
cor(site_means$Total_SR, site_means$latitude)


n.sites <- length(unique(site_data_merge$abbrev))
site.list <- as.character(unique(site_data_merge$abbrev))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp,]
  if (nrow(site_data_merge_temp) > 1){
  lm.temp <- lm(site_data_merge_temp$rarefied.alpha ~ site_data_merge_temp$count_yr)
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

#Create site list for post 1975 sites 
site_list_75 <- site_data$unique_id[site_data$count_yr >= 80] 
cbc_cast_75 <- cbc_cast[rownames(cbc_cast) %in% site_list_75,]

#Create matrix for storing Betas 
beta.matrix <- data.frame(unique_id = site_list_75, beta.bc = NA)

#loop through each site B-C

for (n in 1:n.sites){
  #site.beta.temp <- vector()
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(cbc_df[which(cbc_df$abbrev == site.temp),"count_yr"])
  yrs.temp <- yrs.temp[yrs.temp >= 80]
  #only consider sites that were observed in more than 5 years 
  if (length(yrs.temp) > 5){
  #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years 
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1] + 4)]
    #This is a check but slows down the loop 
    #print(paste0("Number of years for site ", site.temp, " is ", length(first.yrs)))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years 
    cbc_cast_first <- cbc_cast[rownames(cbc_cast) %in% site.yrs.first,]
    first.comm.temp <- colMeans(cbc_cast_first)
    #get remaining years 
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1] + 4)]
    #Now loop the remaining years 
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y]
      other.site.temp <- paste0(site.temp, "_", other.yr.temp)
      cbc_cast_other <- cbc_cast[rownames(cbc_cast) == other.site.temp,]
      cbc_cast_other <- rbind(cbc_cast_other, first.comm.temp)
      beta.temp <- vegdist(sqrt(sqrt(cbc_cast_other)), method = "bray")
      beta.matrix$beta.bc[beta.matrix$unique_id == other.site.temp] <- beta.temp
    }
    }}
    

#Compute Beta-diversity with Jaccard method
beta.jaccard <- data.frame(unique_id = site_list_75, beta.jac = NA)

for (n in 1:n.sites){
  #site.beta.temp <- vector()
  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
  yrs.temp <- unique(cbc_df[which(cbc_df$abbrev == site.temp),"count_yr"])
  yrs.temp <- yrs.temp[yrs.temp >= 80]
  #only consider sites that were observed in more than 5 years 
  if (length(yrs.temp) > 5){
    #Create a baseline year from years temp
    #Sorting
    yrs.temp <- yrs.temp[order(yrs.temp)]
    #find the first 5 years 
    first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1] + 4)]
    #This is a check but slows down the loop 
    #print(paste0("Number of years for site ", site.temp, " is ", length(first.yrs)))
    site.yrs.first <- paste0(site.temp, "_", first.yrs)
    #create average community of first 5 years 
    cbc_cast_first <- cbc_cast[rownames(cbc_cast) %in% site.yrs.first,]
    first.comm.temp <- colMeans(cbc_cast_first)
    #get remaining years 
    other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1] + 4)]
    #Now loop the remaining years 
    for (y in 1:length(other.yrs.temp)){
      other.yr.temp <- other.yrs.temp[y]
      other.site.temp <- paste0(site.temp, "_", other.yr.temp)
      cbc_cast_other <- cbc_cast[rownames(cbc_cast) == other.site.temp,]
      cbc_cast_other <- rbind(cbc_cast_other, first.comm.temp)
      beta.temp <- vegdist(sqrt(sqrt(cbc_cast_other)), method = "jaccard")
      beta.jaccard$beta.jac[beta.jaccard$unique_id == other.site.temp] <- beta.temp
    }
  }}


#Put betas back in site daTA MERGE
site_data_merge <- merge(site_data_merge, beta.matrix, by = "unique_id")
site_data_merge <- merge(site_data_merge, beta.jaccard, by = "unique_id")
site_bin_agg <- aggregate(list(beta.bc = site_data_merge$beta.bc, beta.jac = site_data_merge$beta.jac),
                          by = list(Abbrev = site_data_merge$abbrev, yr_bin = site_data_merge$yr_bin), mean)
site_bin_agg$unique_id <- paste0(site_bin_agg$Abbrev, "_", site_bin_agg$yr_bin)                  
site_bin_agg <- site_bin_agg[complete.cases(site_bin_agg), ]


#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, beta.slope = NA, beta.slope.jac = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp,]
  site_data_merge_temp <- site_data_merge_temp[complete.cases(site_data_merge_temp),]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$beta.bc ~ site_data_merge_temp$count_yr)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
  }
}

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp,]
  site_data_merge_temp <- site_data_merge_temp[complete.cases(site_data_merge_temp),]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$beta.jac ~ site_data_merge_temp$count_yr)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,3] <- slope.temp
  }
}

hist(slopes_sites_beta$beta.slope)
t.test(slopes_sites_beta$beta.slope, mu = 0)

hist(slopes_sites_beta$beta.slope.jac)
t.test(slopes_sites_beta$beta.slope.jac, mu = 0)



###Create a beta diversity for each 5 year bin 
site_data_merge$yr_bin <- NA
site_data_merge$yr_bin[site_data_merge$Year <= 1984] <- 1
site_data_merge$yr_bin[site_data_merge$Year >= 1985 & site_data_merge$Year <= 1989] <- 2
site_data_merge$yr_bin[site_data_merge$Year >= 1990 & site_data_merge$Year <= 1994] <- 3
site_data_merge$yr_bin[site_data_merge$Year >= 1995 & site_data_merge$Year <= 1999] <- 4
site_data_merge$yr_bin[site_data_merge$Year >= 2000 & site_data_merge$Year <= 2004] <- 5
site_data_merge$yr_bin[site_data_merge$Year >= 2005 & site_data_merge$Year <= 2009] <- 6
site_data_merge$yr_bin[site_data_merge$Year >= 2010 & site_data_merge$Year <= 2014] <- 7
site_data_merge$yr_bin[site_data_merge$Year >= 2015 & site_data_merge$Year <= 2019] <- 8

# site_data_merge$site_bin <- paste0(site_data_merge$abbrev, "_", site_data_merge$yr_bin)
# site_data_merge_beta <- group_by(site_data_merge, site_bin) %>% summarise(beta_bin = mean(beta)) %>% cbind
# site_data_merge <- merge(site_data_merge_beta, site_data_merge, by = "site_bin")
# fl_beta_sub <- site_data_merge_beta %>% filter(str_detect(site_bin, "FL"))



####Plots for Beta Diversity##### 

site_data_merge <- merge(site_data_merge, cbc_plot, by = "unique_id")
dev.off()
gghist_beta <- ggplot(data = slopes_sites_beta, aes(slopes_sites_beta$beta.slope.jac)) + 
  geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  labs(title = "") +
  labs(x = expression(paste("Jaccard Slopes of ", beta, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2)
  # coord_flip()
gghist_beta

#site_data_merge$Year <- site_data_merge$count_yr + 1900
site_data_merge <- separate(site_data_merge, Abbrev, c("subnational_code", "site"), sep = -2)

plot_beta <- (ggplot(site_data_merge, aes(x = Year, y = beta.jac, group = abbrev, 
                                           colour = factor(subnational_code, 
                                                           labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = lm, se = FALSE, aes(x = Year, y = beta.jac, group = abbrev, 
                                                          colour = factor(subnational_code, 
                                                                          labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 xlab("Year") +
                 ylab(expression(paste("Jaccard ", beta, "-diversity"))) +
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
  draw_plot(gghist_beta, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 
  
#find the mean b-c for the site
#    mean.beta <- mean(site.beta.temp)
#    mean.betas[n,2] <- mean.beta
  

#Species Richness Trends GoM
#Add a column to Summary_Stats to indicate 5 year bin to look at long term trends
#start at 1975 as this was the year standard protocols were implemented

Summary_stats_75 <- Summary_stats[Summary_stats$Year >= 1980,]
Summary_stats_75$Yr_bin <- 1
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 1985 & Summary_stats_75$Year <= 1989] <- 2
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 1990 & Summary_stats_75$Year <= 1994] <- 3
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 1995 & Summary_stats_75$Year <= 1999] <- 4
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 2000 & Summary_stats_75$Year <= 2004] <- 5
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 2005 & Summary_stats_75$Year <= 2009] <- 6
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 2010 & Summary_stats_75$Year <= 2014] <- 7
Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 2015 & Summary_stats_75$Year <= 2019] <- 8
# Summary_stats_75$Yr_bin[Summary_stats_75$Year >= 2015 & Summary_stats_75$Year <= 2019] <- 9


Mean_alpha_bin <- aggregate(data = Summary_stats_75, rarefied.alpha ~ Yr_bin + Abbrev, FUN = "mean")

#plot(Mean_alpha_bin$rarefied.alpha ~ Mean_alpha_bin$Yr_bin)
lm_alpha <- lm(rarefied.alpha~ Yr_bin + Abbrev, data = Mean_alpha_bin)
summary(lm_alpha)
anova(lm_alpha)
library(lme4)
lmer_alpha <- lmer(rarefied.alpha~ Yr_bin + (1|Abbrev), data = Mean_alpha_bin)
summary(lmer_alpha)
anova(lmer_alpha)
dev.off()

#Slopes for mean 1975 rarefied alphas 
n.sites <- length(unique(site_data_merge$abbrev))
site.list <- as.character(unique(site_data_merge$abbrev))

slopes_sites <- data.frame(site.list, slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  site_data_merge_temp <- site_data_merge[site_data_merge$abbrev == site.temp & site_data_merge$count_yr >= 80, ]
  if (nrow(site_data_merge_temp) > 1){
    lm.temp <- lm(site_data_merge_temp$rarefied.alpha ~ site_data_merge_temp$count_yr)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites[i,2] <- slope.temp
  }
}

plot(slopes_sites$slope ~ slopes_sites$site.list)
abline(h = 0)

t.test(slopes_sites$slope, mu = 0)



##Making data frame for map in ArcGIS w/ sites, betas, and rarefied alpha 
slope.merge <- merge(slopes_sites, slopes_sites_beta, by = "site.list")
cbc_spatial <- cbc[, c("Abbrev", "Latitude", "Longitude")]
gis.df <- merge(cbc_spatial, slope.merge, by.x = "Abbrev", by.y = "site.list")
gis.df$Abbrev <- as.character(gis.df$Abbrev)
gis.df <- unique(gis.df)
gis.df <- gis.df[!is.na(gis.df$beta.slope),]
gis.df$slope_bin <- 1
gis.df$slope_bin[gis.df$slope >= -0.371463455 & gis.df$slope <= -0.0759999] <- 1
gis.df$slope_bin[gis.df$slope >= -0.076 & gis.df$slope <= 0.1299999] <- 2
gis.df$slope_bin[gis.df$slope >= 0.130 & gis.df$slope <= 0.3219999] <- 3
gis.df$slope_bin[gis.df$slope >= 0.322 & gis.df$slope <= 0.6759999] <- 4
gis.df$slope_bin[gis.df$slope >= 0.676 & gis.df$slope <= 1.412] <- 5
write.csv(gis.df, file = "gis.df")

##Time-series plots for alpha diversity
hist(slopes_sites$slope, breaks = 8, xlim = c(-1.5, 2.0), ylim = c(0, 50), main = "", 
     xlab = expression(paste("Slopes of ", alpha, "-diversity")), ylab = "# Sites")

gghist <- ggplot(data = slopes_sites, aes(slopes_sites$slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 30, binwidth = 0.25) +
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", alpha, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()

site_data_merge$Year <- site_data_merge$count_yr + 1900


plot_alpha <- (ggplot(site_data_merge, aes(x = Year, y = rarefied.alpha, group = abbrev, 
                            colour = factor(subnational_code, 
                                            labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
  #geom_point(size = 3) +
  #geom_line()+
  geom_smooth(method = lm, se = FALSE, aes(x = Year, y = rarefied.alpha, group = abbrev, 
                                           colour = factor(subnational_code, 
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

ggdraw() + 
  draw_plot(plot_alpha + theme(legend.justification = "top"), 
                               x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 
 



####Beta diversity using Baseline of first 5 years###
#Code below was placed in different area for running the script (late addition to the code here)

#site_list_75 <- site_data$unique_id[site_data$count_yr >= 75] 
#cbc_cast_75 <- cbc_cast[rownames(cbc_cast) %in% site_list_75,]

#for(i in 1:n.sites){
#  site.temp.beta <-  
#}

####Do just 1985 tp 2017
#This code is likely extraneous 


#mean temporal beta diversity by site
#long.term.betas <- data.frame(site.list, long.term.beta = NA)

#loop through each site
#for (n in 1:n.sites){
#  site.temp <- site.list[n]
  #Find all the years for which that site was surveys
#  yrs.temp <- unique(cbc_unique[which(cbc_unique$abbrev == site.temp),"count_yr"])
  #only consider sites that surveyed in 1985
#  if (sum(yrs.temp == 85, yrs.temp == 117)==2){
#      site.yrs <- c(85, 117) #pick out two years
      #add to site.temp to get unique.ids
#      sites.yrs.temp <- paste0(site.temp, "_", site.yrs)
      #find the row indices in the site x spp matrix that match the pair of surveys
#      ind.site.temp <- which(rownames(cbc_cast) %in% sites.yrs.temp)
      #subset the species x sites matrix for a given site only
#      site_spp_temp <- cbc_cast[ind.site.temp,]
      #calculate pairwise B-C
#      dissim.temp <- vegdist(site_spp_temp, method="bray", diag=FALSE, upper=TRUE)
#      long.term.betas[n,2] <- dissim.temp
#  }  
#}


#site_means3 <- merge(site_means2, long.term.betas)
#plot(site_means3$long.term.beta ~ site_means3$latitude)
#text(site_means3$latitude , site_means3$long.term.beta, site_means3$abbrev, cex = 0.8)
#abline(lm(site_means3$long.term.beta ~ site_means3$latitude))

















#Temperature vs. Latitude

cbc_envi <- read.csv("Keyser_CBC_Weather_Original_CSV.csv", header = TRUE)

names(cbc_envi) <- tolower(names(cbc_envi))

envi_unique <- cbc_envi %>% unite(unique_id, abbrev, count_yr, remove = FALSE)

#Clean up tmin
unique(envi_unique$min_temp)
envi_bad <- is.na(envi_unique$min_temp)
envi_unique$min_temp[!envi_bad]
envi_unique$min_temp <- as.character(envi_unique$min_temp)
envi_unique$min_temp <- as.numeric(envi_unique$min_temp)


#Check status of envi_unique
str(envi_unique)

#New variable of the mean temperatures
envi_unique$mean_t <- (envi_unique$max_temp + envi_unique$min_temp) / 2

#Means at sites across years
t_simple <- envi_unique[,c("abbrev", "min_temp", "max_temp", "mean_t")]
t_cast <- dcast(data = t_simple, formula = abbrev ~ mean_t, fun.aggregate = mean)

#Mean at each site for tmax and tmin
##mean_tmax <- aggregate(max_temp ~ unique_id, data = envi_unique, FUN = "mean", na.rm = TRUE)
##mean_tmin <- aggregate(min_temp ~ unique_id, data = envi_unique, FUN = "mean", na.rm = TRUE)
##envi_unique$mean_temps <- rowMeans(envi_unique[c('max_temp', 'min_temp')], na.rm=TRUE)
##envi.mean <- aggregate(cbind(max_temp, min_temp) ~ unique_id, data = envi_unique, FUN = "mean", na.rm = TRUE)
##envi.big <- merge(envi_unique, envi.mean, by = "unique_id")



