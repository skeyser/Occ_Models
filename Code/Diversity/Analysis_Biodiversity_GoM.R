#################################################################################
#################################################################################
############Script for running final models for GoM biodiversity proj############
#################################################################################
#########################Created by: Spencer R Keyser############################
###################################9/25/2019#####################################
#################################################################################

#Load Packages

library("pacman")

pacman::p_load("here", "tidyverse", "MuMIn", "sjPlot", "lme4", "vegan")

#################################################################################

#Load in data sets for analysis 

#Segments level df .4km buffer
#Need to change DF to include ratio mangrove:ww and diff.from.first.wwnm
seg.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Seg_Last.csv"), stringsAsFactors = F)

#Route level df 0.4km buffer
rt.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Rt_Last.csv"), stringsAsFactors = F)

#Occupancy Data for MDS plot & PERMANOVA 
occ <- read.csv(here:here("Data_BBS/Generated DFs/occ50.csv"), stringsAsFactors = F)

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

bbs.mds <- bbs_simple[bbs_simple$Year == 1980,]
bbs.mds <- bbs.mds %>% arrange(site)

lc.mds <- bbs_lulc %>% dplyr::select(site, Year, ratio.ww, pct.ur)
lc.mds <- lc.mds[lc.mds$Year == 1980,]

mds.list <- unique(bbs.mds$site)
wetland.site <- unique(lc.mds$site)

mds.full <- merge(lc.mds, bbs.mds, by = "site")

lc.mds <- dplyr::select(mds.full, c(ratio.ww, ratio.man, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$groups <- NA
ww.mds$groups[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
ww.mds$groups[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds$groups

bbs.mds <- bbs.mds %>% arrange(site)
mds_cast <- dcast(mds.full, unique_id ~ sci_name, value.var = "Detected", fun.aggregate = sum)
rownames(mds_cast) <- mds_cast$unique_id
mds_cast <- mds_cast[, -1]


mds <- metaMDS(mds_cast, trymax = 50)


data.scores <- as.data.frame(scores(mds))  
data.scores$site <- rownames(data.scores)  
data.scores$grp1 <- ww.mds

mds_plot1 <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = ww.mds), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = ww.mds, shape = ww.mds), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"), # "Mix",
                     values = c(0, 16))

mds_plot1

adon.results <- adonis(mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)


########################################################################################
#########################MDS for the Occupancy Cutoffs##################################
########################################################################################

occ.mds <- occ[occ$Year == 1980,]
occ.mds <- occ.mds %>% arrange(site)

lc.mds <- bbs_lulc %>% dplyr::select(site, Year, ratio.ww)
lc.mds <- lc.mds[lc.mds$Year == 1980,]

occ.mds.list <- unique(occ.mds$site)
wetland.site <- unique(lc.mds$site)

occ.mds.full <- merge(lc.mds, occ.mds, by = "site")

lc.mds <- dplyr::select(occ.mds.full, c(ratio.ww, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$groups <- NA
ww.mds$groups[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
#ww.mds$groups[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mixed"
ww.mds$groups[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds$groups

occ.mds <- occ.mds %>% arrange(site)
occ_mds_cast <- reshape2::dcast(occ.mds.full, unique_id ~ Species, value.var = "Occupancy", fun.aggregate = sum)
rownames(occ_mds_cast) <- occ_mds_cast$unique_id
occ_mds_cast <- occ_mds_cast[, -1]


mds <- metaMDS(occ_mds_cast, trymax = 50)


data.scores <- as.data.frame(scores(mds))  
data.scores$site <- rownames(data.scores)  
data.scores$grp1 <- ww.mds
#data.scores$grp2 <- wet.mds


mds_plot <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = ww.mds), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = ww.mds, shape = ww.mds), size=2) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"),
                      values = c("#55C667FF", "#440154FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Woody Wetland Dominated"), # "Mix",
                     values = c(0, 16))
mds_plot #+ labs(color = "Wetland Cover Types", shape = "Wetland Cover Type") + scale_shape_manual(values = c(0, 16, 3))


adon.results <- adonis(occ_mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)


#Linear Mixed Effects Models @ segment level

#Jaccard Index Wetland and Entire Community
mod1 <- lmer(data = seg.df, beta50.jac ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

mod2 <- lmer(data = seg.df, beta.wet.jac ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

#Turnover Component Wetland and Entire Community
mod3 <- lmer(data = seg.df, beta50.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

mod4 <- lmer(data = seg.df, beta.wet.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

#Nestedness Compoentn Wetland and Entire Community
mod5 <- lmer(data = seg.df, beta50.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

mod6 <- lmer(data = seg.df, beta.wet.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + SR50 + Duration + (1|rteno.x.x), REML = F)

#Change SR 
mod7 <- lmer(data = seg.df, change.alpha ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration + (1|rteno.x.x), REML = F)

mod8 <- lmer(data = seg.df, change.wet.alpha ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration + (1|rteno.x.x), REML = F)

########################################################################################################################################################################################################

#Multiple Regression Models @ Route level

#Jaccard Index Models
mod9 <- lm(data = rt.df, beta50.jac ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

mod10 <- lm(data = rt.df, beta.wet.jac ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

#Turnover Index Models
mod11 <- lm(data = rt.df, beta50.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

mod12 <- lm(data = rt.df, beta.wet.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

#Nestedness Index Models 
mod13 <- lm(data = rt.df, beta50.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

mod14 <- lm(data = rt.df, beta.wet.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

#Change SR
mod15 <- lmer(data = rt.df, change.alpha ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

mod16 <- lmer(data = rt.df, change.wet.alpha ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ww + diff.from.first.ur + diff.from.first.ag + Duration)

