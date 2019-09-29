#################################################################################
#################################################################################
############Script for running final models for GoM biodiversity proj############
#################################################################################
#########################Created by: Spencer R Keyser############################
###################################9/25/2019#####################################
#################################################################################

#Load Packages

library("pacman")

pacman::p_load("here", "tidyverse", "MuMIn", "sjPlot", "lme4", "vegan", "viridis")

#################################################################################

#Load in data sets for analysis 

#Segments level df .4km buffer
#Need to change DF to include ratio mangrove:ww and diff.from.first.wwnm
seg.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Seg_Last.csv"), stringsAsFactors = F)

#Route level df 0.4km buffer
rt.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Rt_Last.csv"), stringsAsFactors = F)

#Occupancy Data for MDS plot & PERMANOVA 
occ <- read.csv(here::here("Data_BBS/Generated DFs/occ50.csv"), stringsAsFactors = F)

#LULC Data for MDS plot & PERMANOVA
bbs_lulc <- read.csv(here::here("Data_Envi/DF4Analysis_LULC.csv"), stringsAsFactors = F)

#Climate and Diversity only DF Segment Level
bbs_clim <- read.csv(here::here("Data_BBS/Generated DFs/DF_div_clim.csv"), stringsAsFactors = F)

#Climate and Diversity only DF Rt Level
bbs_clim_rt <- read.csv(here::here("Data_BBS/Generated DFs/DF_div_clim_rtlvl.csv"), stringsAsFactors = F)

#BBS Simple Df 
bbs_simple <- read.csv(here::here("Data_BBS/Generated DFs/BBS_Simple.csv"), stringsAsFactors = F)

###################################################################################

#MDS and PERMANOVA Analysis 

bbs.mds <- bbs_simple[bbs_simple$Year <= 1990,]
bbs.mds <- bbs.mds %>% arrange(site)

Bins <- data.frame(matrix(ncol = 2, nrow = 8))
l <- c("Year", "Yr_bin")
colnames(Bins) <- l
Bins$Year <- seq(from = 1980, to = 2015, 5)
Bins$Yr_bin <- seq(1:8)

bbs_lulc <- merge(Bins, bbs_lulc, by = "Yr_bin")

lc.mds <- bbs_lulc %>% dplyr::select(Site, Year, Emergent_Wetlands, Woody_Wetlands, mangrove, WW_NoMan)
lc.mds <- lc.mds[lc.mds$Year <= 1990,]

lc.mds <- lc.mds %>% group_by(Site) %>% rename(site = Site) %>% summarise(Emergent_Wetlands = mean(Emergent_Wetlands), Woody_Wetlands = mean(Woody_Wetlands),
                                                                          mangrove = mean(mangrove), WW_NoMan = mean(WW_NoMan)) %>% 
  mutate(ratio.ww = (Woody_Wetlands / Emergent_Wetlands), ratio.ew = (Emergent_Wetlands / Woody_Wetlands), ratio.man = (mangrove / Emergent_Wetlands)) %>% ungroup()

mds.list <- unique(bbs.mds$site)
wetland.site <- unique(lc.mds$site)

mds.full <- merge(lc.mds, bbs.mds, by = "site")

lc.mds <- dplyr::select(mds.full, c(ratio.ww, ratio.man, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$group1 <- NA
ww.mds$group1[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
#ww.mds$groups[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mixed"
ww.mds$group1[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds[, c("site", "ratio.ww", "group1")]

man.mds <- lc.mds %>% arrange(site)
man.mds$group2 <- NA
man.mds$group2[man.mds$ratio.man <= 1] <- "Non-Mangrove Dominated Wetland"
man.mds$group2[man.mds$ratio.man > 1] <- "Mangrove Dominated Wetland"
man.mds <- man.mds[, c("site", "ratio.man", "group2")]

bbs.mds <- bbs.mds %>% arrange(site)
mds_cast <- dcast(mds.full, site ~ sci_name, value.var = "Detected", fun.aggregate = sum)
rownames(mds_cast) <- mds_cast$site
mds_cast <- mds_cast[, -1]
mds_cast[mds_cast > 1] <- 1

mds <- metaMDS(mds_cast, trymax = 500)


data.scores <- as.data.frame(scores(mds))  
data.scores$site <- rownames(data.scores)  
data.scores <- merge(data.scores, ww.mds, by = "site") 
data.scores <- merge(data.scores, man.mds, by = "site")

mds_plot1 <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group1), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group1, shape = group1), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"), # "Mix",
                     values = c(0, 16))

mds_plot1

mds_plot2 <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group2), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group2, shape = group2), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Non-mangrove Dominated Wetlands", "Mangrove Dominated Wetland"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                    labels = c("Non-mangrove Dominated Wetlands", "Mangrove Dominated Wetland"), # "Mix",
                    values = c(0, 16))

mds_plot2


adon.results <- adonis(mds_cast ~ data.scores$group1, method = "jaccard", perm = 999)

adon.results2 <- adonis(mds_cast ~ data.scores$group2, method = "jaccard", perm = 999)


print(adon.results)
print(adon.results2)


########################################################################################
#########################MDS for the Occupancy Cutoffs##################################
########################################################################################

occ.mds <- occ[occ$Year <= 1990,]
occ.mds <- occ.mds %>% arrange(site)

lc.mds <- bbs_lulc %>% dplyr::select(Site, Year, Emergent_Wetlands, Woody_Wetlands, mangrove, WW_NoMan)
lc.mds <- lc.mds[lc.mds$Year <= 1990,]

lc.mds <- lc.mds %>% group_by(Site) %>% rename(site = Site) %>% summarise(Emergent_Wetlands = mean(Emergent_Wetlands), Woody_Wetlands = mean(Woody_Wetlands),
                                                                          mangrove = mean(mangrove), WW_NoMan = mean(WW_NoMan)) %>% 
          mutate(ratio.ww = (Woody_Wetlands / Emergent_Wetlands), ratio.ew = (Emergent_Wetlands / Woody_Wetlands), ratio.man = (mangrove / Emergent_Wetlands)) %>% ungroup()

occ.mds.list <- unique(occ.mds$site)
wetland.site <- unique(lc.mds$site)

occ.mds.full <- merge(lc.mds, occ.mds, by = "site")

lc.mds <- dplyr::select(occ.mds.full, c(ratio.ww, ratio.ma, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$group1 <- NA
ww.mds$group1[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
#ww.mds$groups[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mixed"
ww.mds$group1[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds[, c("site", "ratio.ww", "group1")]

man.mds <- lc.mds %>% arrange(site)
man.mds$group2 <- NA
man.mds$group2[man.mds$ratio.man <= 1] <- "Non-Mangrove Dominated Wetland"
man.mds$group2[man.mds$ratio.man > 1] <- "Mangrove Dominated Wetland"
man.mds <- man.mds[, c("site", "ratio.man", "group2")]

occ.mds <- occ.mds %>% arrange(site)
occ_mds_cast <- reshape2::dcast(occ.mds.full, site ~ Species, value.var = "Occupancy", fun.aggregate = sum)
rownames(occ_mds_cast) <- occ_mds_cast$site
occ_mds_cast <- occ_mds_cast[, -1]
occ_mds_cast[occ_mds_cast > 1] <- 1

mds <- metaMDS(occ_mds_cast, trymax = 500)


data.scores <- as.data.frame(scores(mds))  
data.scores$site <- rownames(data.scores)  
data.scores <- merge(data.scores, ww.mds, by = "site") 
data.scores <- merge(data.scores, man.mds, by = "site")


mds_plot <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group1), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group1, shape = group1), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"), # "Mix",
                     values = c(0, 16))
mds_plot #+ labs(color = "Wetland Cover Types", shape = "Wetland Cover Type") + scale_shape_manual(values = c(0, 16, 3))

mds_plot2 <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group2), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group2, shape = group2), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Mangrove Dominated Wetlands", "Non-mangrove Dominated Wetlands"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Mangrove Dominated Wetlands", "Non-mangrove Dominated Wetlands"), # "Mix",
                     values = c(0, 16))

mds_plot2



adon.results <- adonis(occ_mds_cast ~ data.scores$group1, method = "jaccard", perm = 999)

adon.results2 <- adonis(occ_mds_cast ~ data.scores$group2, method = "jaccard", perm = 999)

print(adon.results)
print(adon.results2)


#Nestedness Vs Turnover GoM
div.df <- rt.df %>% dplyr::select(c("Site", "beta50.jac", "beta50.turn", "beta50.nest", "beta.wet.jac", "beta.wet.turn", "beta.wet.nest", "BCR")) %>%
          mutate(beta50.ratio = (beta50.turn / beta50.nest), betawet.ratio = (beta.wet.turn / beta.wet.nest), Site = as.character(Site), BCR = as.character(BCR),
                 b50.ratiolog = log(beta50.ratio), bwet.ratiolog = log(betawet.ratio)) %>% arrange(BCR)

div.sum <- div.df %>% group_by(BCR) %>% summarise(turnover = mean(beta50.turn), turnover.sd = sd(beta50.turn), nest = mean(beta50.nest), nest.sd = mean(beta50.nest))

div.melt <- reshape2::melt(div.df, "BCR", c("beta50.turn", "beta50.nest"))
colnames(div.melt) <- c("BCR", "Component", "Beta")

ggplot(data = div.melt, aes(x = BCR, y = Beta, fill = Component)) + geom_boxplot() + scale_fill_viridis(discrete = T) + 

ggplot(data = div.melt, aes(y = Site, x = Beta, color = Index)) + geom_point() + scale_color_viridis("Index", discrete = T) 
ggplot(data = div.df, aes(y = Site, x = bwet.ratiolog, color = BCR)) + geom_point() + scale_color_viridis("BCR", discrete = T) + geom_vline(xintercept = 0)



#Linear Mixed Effects Models @ segment level

#Jaccard Index Wetland and Entire Community
mod1 <- lmer(data = seg.df, beta50.jac ~ p.anom.wet + min.anom.sp + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR50 + Duration + (1|rteno), REML = F)
summary(mod1)
Anova(mod1)
tab_model(mod1)

mod2 <- lmer(data = seg.df, beta.wet.jac ~ p.anom.wet + mean.anom.bird + diff.from.first.man  + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR.wet + Duration + (1|rteno), REML = F)
summary(mod2)
Anova(mod2)
tab_model(mod2)

#Turnover Component Wetland and Entire Community
mod3 <- lmer(data = seg.df, beta50.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR50 + Duration + (1|rteno), REML = F)
summary(mod3)
Anova(mod3)
tab_model(mod3)

mod4 <- lmer(data = seg.df, beta.wet.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR.wet + Duration + (1|rteno), REML = F)
summary(mod4)
Anova(mod4)
tab_model(mod4)

#Nestedness Compoentn Wetland and Entire Community
mod5 <- lmer(data = seg.df, beta50.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR50 + Duration + (1|rteno), REML = F)
summary(mod5)
Anova(mod5)
tab_model(mod5)

mod6 <- lmer(data = seg.df, beta.wet.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR.wet + Duration + (1|rteno), REML = F)
summary(mod6)
Anova(mod6)
tab_model(mod6)

#Change SR 
mod7 <- lmer(data = seg.df, alpha50.pchange ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + Duration + SR50 + (1|rteno), REML = F)
summary(mod7)
Anova(mod7)
tab_model(mod7)

mod8 <- lmer(data = seg.df, alphawet.pchange ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + Duration + SR.wet + (1|rteno), REML = F)
summary(mod8)
Anova(mod8)
tab_model(mod8)

########################################################################################################################################################################################################

#Multiple Regression Models @ Route level

#Jaccard Index Models
mod9 <- lm(data = rt.df, beta50.jac ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur  + SR50 + Duration)
summary(mod9)

mod10 <- lm(data = rt.df, beta.wet.jac ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR.wet + Duration)
summary(mod10)

#Turnover Index Models
mod11 <- lm(data = rt.df, beta50.turn ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR50 + Duration)
summary(mod11)

mod12 <- lm(data = rt.df, beta.wet.turn ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR.wet + Duration)
summary(mod12)

#Nestedness Index Models 
mod13 <- lm(data = rt.df, beta50.nest ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR50 + Duration)
summary(mod13)

mod14 <- lm(data = rt.df, beta.wet.nest ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR.wet + Duration)
summary(mod14)

#Change SR
mod15 <- lm(data = rt.df, alpha50.pchange ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + Duration)
summary(mod15)

mod16 <- lm(data = rt.df, alphawet.pchange ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + Duration)
summary(mod16)

########################################################################################################################################################################################################
bbs_clim_rt1 <- bbs_total %>% group_by(rteno.x) %>% arrange(Year) %>% filter(Year <= 1999) %>% 
  select(c("rteno.x", "Year", "Latitude", "SR50", "SR.wet")) %>% summarize_if(is.numeric, mean) %>% ungroup() 
bbs_clim_rt2 <- bbs_total %>% group_by(rteno.x) %>% arrange(Year) %>% filter(Year >= 2000) %>% 
  select(c("rteno.x", "Year", "Latitude", "SR50", "SR.wet")) %>% summarize_if(is.numeric, mean) %>% ungroup() 

bbs_clim_rt2 <- bbs_clim_rt2[bbs_clim_rt2$rteno.x %in% bbs_clim_rt1$rteno.x, ]
bbs_clim_rt1 <- bbs_clim_rt1[bbs_clim_rt1$rteno.x %in% bbs_clim_rt2$rteno.x, ]


ggplot() + 
  geom_smooth(data = bbs_clim_rt1, aes(x = Latitude, y = SR50), color = "red") +
  geom_smooth(data = bbs_clim_rt2, aes(x = Latitude, y = SR50), color = "black") +
  xlab("Latitude") +
  ylab("Detection-corrected Species Richness")
  

########################################################################################################################################################################################################

#-Wetland, -WW, +Ag  
mod22 <- lmer(data = seg.df, SR ~ pct.wetland + (1|rteno), REML = F)
summary(mod22)
Anova(mod22)

#-WW, +EW, -WWNM, +Ag
mod23 <- lm(data = rt.df, SR.wet ~ pct.wwnm)
summary(mod23)
Anova(mod23)
tab_model(mod23)

########################################################################################################################################################################################################
#CMRL - Beta div plots
mod17 <- lm(data = rt.df, beta50.jac ~ change.cmrl)
summary(mod17)

mod18 <- lm(data = rt.df, beta50.turn ~ change.cmrl)
summary(mod18)

mod19 <- lm(data = rt.df, beta50.nest ~ change.cmrl)
summary(mod19)

mod20 <- lmer(data = seg.df, change.cmrl ~ min.anom.sp + poly(min.anom.sp, 2) + (1|rteno), REML = F)
summary(mod20)
Anova(mod20)
tab_model(mod20)

mod21 <- lm(data = rt.df, change.cmrl.km ~ min.anom.sp + poly(min.anom.sp, 2))
summary(mod21)

prd <- data.frame(min.anom.sp = seq(from  = range(rt.df$min.anom.sp)[1], to = range(rt.df$min.anom.sp)[2], length.out = 100))
err <- predict(mod21, newdata = prd, se.fit = T)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

cmrl.plot <- ggplot(prd, aes(x = min.anom.sp, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#481567FF", fill = "#481567FF") +
  geom_point(data = rt.df, aes(x = min.anom.sp, y = change.cmrl.km)) +
  xlab(expression(paste("Change in Minimum Spring Temp ( ", degree ~ C, " )"))) +
  ylab('Change in CMRL (km)') +
  labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
                     "Intercept =",signif(mod21$coef[[1]],5 ),
                     " Slope =",signif(mod21$coef[[3]], 5),
                     " P =",signif(summary(mod21)$coef[3,4], 5))) +
  scale_x_continuous(breaks = seq(-0.5, 2.5, 0.5))
  
cmrl.plot
#######################################################################################################
##########################Plotting Function for partial regression#####################################
#######################################################################################################

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "black") +
    labs(title = paste("R2 = ",signif(summary(fit)$r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), 
         x = expression(paste("Residual ", alpha, "-Diversity")), 
         y = expression(paste("Residual ", beta, "-Diversity")))}

#######################################################################################################


#Mapping distribution of Beta-diversity across GoM

states <- map_data("state")
dim(states)
head(states)
gom <- subset(states, region %in% c("texas", "florida", "alabama", "arkansas", "mississippi", "louisiana", "georgia"))

map_gom_jac <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.jac), size = 3)

map_gom_jac + scale_color_viridis(name = expression(paste(beta, "-diversity"))) + theme_map() + theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  

map_gom_turn <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.turn), size = 3)

map_gom_turn + scale_color_viridis(option = 'magma', name = expression(paste(beta, "-diversity"))) + theme_map() + theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  

map_gom_nest <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.nest), size = 3)

map_gom_nest + scale_color_viridis(option = 'plasma', name = expression(paste(beta, "-diversity"))) + theme_map() + theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  

