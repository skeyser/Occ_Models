#################################################################################
#################################################################################
############Script for running final models for GoM biodiversity proj############
#################################################################################
#########################Created by: Spencer R Keyser############################
###################################9/25/2019#####################################
#################################################################################

#Load Packages

library("pacman")

pacman::p_load("here", "tidyverse", "MuMIn", "sjPlot", "lme4", "vegan", "viridis", "ggmap", "maps", "ggfortify", "cowplot", "extrafont", "scales", "betareg")

#################################################################################

#################################################################################
#Checking Color Palettes for manual color changes 
q_colors <- 20
v_colors <- viridis(q_colors, option = "E")
show_col(v_colors)
loadfonts(device = "win")
windowsFonts()


#Load in data sets for analysis 

#Segments level df .4km buffer
#Need to change DF to include ratio mangrove:ww and diff.from.first.wwnm
#seg.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Seg_Last.csv"), stringsAsFactors = F)
seg.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Seg_Last_New.csv"), stringsAsFactors = F)

#Route level df 0.4km buffer
#rt.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Rt_Last.csv"), stringsAsFactors = F)
rt.df <- read.csv(here::here("Data_BBS/Generated DFs/DF4Analysis_Rt_Last_New.csv"), stringsAsFactors = F)

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

#BBS Total Diversity Seg Level
bbs_total_seg <- read.csv(here::here("Data_BBS/Generated DFs/bbs_div_df_seg.csv"), stringsAsFactors = F)
bbs_slopes_seg <- read.csv(here::here("Data_BBS/Generated DFs/bbs_slopes_seg.csv"), stringsAsFactors = F)

#BBS Total Diversity Rt Level
bbs_total_rt <- read.csv(here::here("Data_BBS/Generated DFs/bbs_diversity_df_rt.csv"), stringsAsFactors = F)
bbs_slopes_rt <- read.csv(here::here("Data_BBS/Generated DFs/bbs_slopes_rt.csv"))

#Rt XY
rtxy <- bbs_clim_rt[, c("Site.x", "Longitude", "Latitude")]
rtxy <- rtxy[!duplicated(rtxy),]

#Merge XY with rt.df
rt.df <- merge(rt.df, rtxy, by.x = "Site", by.y = "Site.x")

#Climate Data only
climate <- read.csv(here::here("Data_BBS/Generated DFs/climate_data.csv"),stringsAsFactors = F)
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

lc.mds <- bbs_lulc %>% dplyr::select(Site, Yr_bin, Emergent_Wetlands, Woody_Wetlands, mangrove, WW_NoMan, pct.man)
lc.mds <- lc.mds[lc.mds$Yr_bin <= 3,]

lc.mds <- lc.mds %>% group_by(Site) %>% rename(site = Site) %>% summarise(Emergent_Wetlands = mean(Emergent_Wetlands), Woody_Wetlands = mean(Woody_Wetlands),
                                                                          mangrove = mean(mangrove), WW_NoMan = mean(WW_NoMan), pct.man = mean(pct.man)) %>% 
          mutate(ratio.ww = (WW_NoMan / Emergent_Wetlands), ratio.ew = (Emergent_Wetlands / Woody_Wetlands), ratio.man = (mangrove / Emergent_Wetlands)) %>% ungroup()

occ.mds.list <- unique(occ.mds$site)
wetland.site <- unique(lc.mds$site)

occ.mds.full <- merge(lc.mds, occ.mds, by = "site")

lc.mds <- dplyr::select(occ.mds.full, c(ratio.ew, ratio.ww, ratio.man, pct.man, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$group1 <- NA
ww.mds$group1[ww.mds$ratio.ww < 0.5] <- "Emergent Wetland Dominated"
ww.mds$group1[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mixed Wetland"
ww.mds$group1[ww.mds$ratio.ww > 1.5] <- "Woody Wetland Dominated"
ww.mds <- ww.mds[, c("site", "ratio.ww", "ratio.ew", "group1")]

man.mds <- lc.mds %>% arrange(site)
man.mds$group2 <- NA
man.mds$group2[man.mds$ratio.man == 0] <- "Non-Mangrove Dominated Wetland"
man.mds$group2[man.mds$ratio.man > 0 & man.mds$ratio.man <= 1] <- "Mixed Wetland"
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
data.scores <- data.scores %>% separate(site, into = c("rteno", "seg"), sep = "_")
data.scores$site <- paste0(data.scores$rteno, "_", data.scores$seg)


mds_plot <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group1), level = 0.5, size = 1.2) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group1, shape = group1), size=3) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Emergent Wetland Dominated", "Mixed Wetland", "Woody Wetland Dominated"),
                      values = c("#FAC62DFF", "#781C6DFF", "#0C0826FF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Mixed Wetland", "Woody Wetland Dominated"), # "Mix",
                     values = c(15, 16, 17)) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black", size =1.2),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(vjust = -1, size = 12, family = "serif"),
        axis.title.y = element_text(vjust = 1.5, size = 12, family = "serif"),
        axis.ticks = element_line(size = 1.2),
        legend.text = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,2,2), "lines"),
        text = element_text(size=12, family = "serif"))
mds_plot #+ labs(color = "Wetland Cover Types", shape = "Wetland Cover Type") + scale_shape_manual(values = c(0, 16, 3))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/MDSplot1.tiff"), plot = mds_plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

mds_plot2 <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = group2), level = 0.5, size = 1.2) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = group2, shape = group2), size=3) +
  scale_colour_manual(name = "Wetland Cover Types",  
                      labels = c("Mangrove Dominated Wetland", "Mangrove-Present Wetland", "Non-mangrove Dominated Wetland"), #+ , discrete = T, option = "E") +
                      values = c("#A69D75FF", "#848279FF", "#00204DFF")) + #"#009E73"
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Mangrove Dominated Wetland", "Mangrove-Present Wetland", "Non-mangrove Dominated Wetland"), # "Mix",
                     values = c(15, 16, 17)) +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black", size =1.2),
        axis.text.x = element_text(size = 12, family = "serif"),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(vjust = -1, size = 12, family = "serif"),
        axis.title.y = element_text(vjust = 1.5, size = 12, family = "serif"),
        axis.ticks = element_line(size = 1.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,2,2), "lines"),
        text = element_text(size=12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"),
        legend.key.size = unit(0.5, "cm"))

mds_plot2

ggsave(here::here("Figures/Figures_Diversity_Manuscript/MDSplot2.tiff"), plot = mds_plot2,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

adon.results <- adonis(occ_mds_cast ~ data.scores$group1, method = "jaccard", perm = 999)

adon.results2 <- adonis(occ_mds_cast ~ data.scores$group2, method = "jaccard", perm = 999)

print(adon.results)
print(adon.results2)

mds_plot_m <- plot_grid(mds_plot2, mds_plot, nrow = 2, align = "hv", labels = c("A","B"), label_fontfamily = "serif")

ggsave(here::here("Figures/Figures_Diversity_Manuscript/MDSplot_combine.tiff"), plot = mds_plot_m,
       device = "tiff", width = 8, height = 8, units = "in", dpi = 600)

dev.off()



#Nestedness Vs Turnover GoM
div.df <- seg.df %>% dplyr::select(c("Site", "beta50.jac", "beta50.turn", "beta50.nest", "beta.wet.jac", "beta.wet.turn", "beta.wet.nest", "BCR")) %>%
          mutate(beta50.nest = (beta50.nest / beta50.jac), beta50.turn = (beta50.turn / beta50.jac), betawet.nest = (beta.wet.nest / beta.wet.jac),
                 beta.wet.turn = (beta.wet.turn / beta.wet.jac), Site = as.character(Site), BCR = as.character(BCR)) 

#div.sum <- div.df %>% group_by(BCR) %>% summarise(turnover = mean(beta50.turn), turnover.sd = sd(beta50.turn), nest = mean(beta50.nest), nest.sd = mean(beta50.nest))

div.melt <- reshape2::melt(div.df, id = c("Site", "BCR"), measure.vars = c("beta50.turn", "beta50.nest"))
colnames(div.melt) <- c("Site", "BCR", "Component", "Proportion")

div.melt <- div.melt %>% group_by(Site) %>% arrange(Component)

div.melt$BCR_name <- NA
div.melt$BCR_name[div.melt$BCR == 26] <- "Mississippi Alluvial Valley"
div.melt$BCR_name[div.melt$BCR == 27] <- "Southeastern Coastal Plain"
div.melt$BCR_name[div.melt$BCR == 31] <- "Peninsular Florida"
div.melt$BCR_name[div.melt$BCR == 37] <- "Gulf Coast Prairie"
  
comp.plot <- ggplot(data = div.melt, aes(x = Site, y = Proportion, fill = Component, group = Site)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("#00204DFF", "#777776FF"),
                    labels = c(expression(paste(" ", beta["Turn"])),
                               expression(beta["Nest"]))) +
  labs(colour = "Component") +
  ylab(expression("Proportion of Total"~~beta)) +
                    #name = expression("Compenent of"~beta*"-Diversity"),
                    #labels = c(expression(paste(beta*["Turn"])), expression(paste(beta*["Nest"])))) + 
  #scale_fill_viridis(discrete = T, option = "cividis") + 
  # theme(axis.text.x = element_text(angle = 90), 
  #       panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
  #       plot.background = element_rect(fill = "white", color = "black"),
  #       legend.title = element_text(size = 12),
  #       panel.grid.major = element_line(size = 0.5, linetype = "solid",),
  #       panel.grid.minor = element_line(size = 0.25, linetype = "solid")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size =1.2),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(vjust = -1, size = 12, family = "serif"),
        axis.title.y = element_text(vjust = 1.5, size = 12, family = "serif"),
        axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,2,2), "lines"),
        text = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 12, family = "serif"),
        legend.spacing = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0,0)) +
  #geom_text(aes(x = 1, y = 1.1, label = "Stretch it"), vjust = -1) +
  facet_wrap(~BCR_name, scales = "free")

comp.plot    
ggsave(here::here("Figures/Figures_Diversity_Manuscript/ComponentsBetaplot.tiff"), plot = comp.plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)



#ggplot(data = div.melt, aes(y = Site, x = Beta, color = Index)) + geom_point() + scale_color_viridis("Index", discrete = T) 
#ggplot(data = div.df, aes(y = Site, x = bwet.ratiolog, color = BCR)) + geom_point() + scale_color_viridis("BCR", discrete = T) + geom_vline(xintercept = 0)



#Linear Mixed Effects Models @ segment level

#Jaccard Index Wetland and Entire Community
colnames(seg.df)[colnames(seg.df) == "rteno"] <- "Route"

log.beta <- function(x){
  log(x / (1 - x))
}
beta.res <- function(x){
  (exp(x) / (1 + exp(x)))
}

seg.df$p.anom.s <- scale(seg.df$p.anom) 
seg.df$p.anom.wet.s <- scale(seg.df$p.anom.wet)
seg.df$p.anom.dry.s <- scale(seg.df$p.anom.dry)
seg.df$p.anom.f.s <- scale(seg.df$p.anom.f)
seg.df$p.anom.s.s <- scale(seg.df$p.anom.s)
seg.df$p.anom.sp.s <- scale(seg.df$p.anom.sp)
seg.df$p.anom.w.s <- scale(seg.df$p.anom.w)
seg.df$mean.anom.bird.s <- scale(seg.df$mean.anom.bird)
seg.df$diff.from.first.man.s <- scale(seg.df$diff.from.first.man)
seg.df$diff.from.first.ew.s <- scale(seg.df$diff.from.first.ew)
seg.df$diff.from.first.human.s <- scale(seg.df$diff.from.first.human)
seg.df$diff.from.first.wwnm.s <- scale(seg.df$diff.from.first.wwnm)
seg.df$diff.from.first.ww.s <- scale(seg.df$diff.from.first.ww)
seg.df$diff.from.first.nat.s <- scale(seg.df$diff.from.first.nat)
seg.df$diff.from.first.wet.s <- scale(seg.df$diff.from.first.wet)
seg.df$Duration.s <- scale(seg.df$Duration)
seg.df$SR50.s <- scale(seg.df$SR50)
seg.df$SR.wet.s <- scale(seg.df$SR.wet)
seg.df$beta50.jac.log <- log.beta(seg.df$beta50.jac)
seg.df$beta.wet.jac.log <- log.beta(seg.df$beta.wet.jac)


mod1 <- lmer(data = seg.df, beta50.jac.log ~ p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.human + Duration + SR50 + (1|Route), REML = F)
summary(mod1) 
car::Anova(mod1)
avPlots(mod1)
# mod1 <- lmer(data = seg.df, beta50.jac ~ p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.human + Duration + SR50 + (1|Route), REML = F)
# summary(mod1) 
# car::Anova(mod1)


tab_model(mod1, show.ci = 0.95, title = NULL, pred.labels = c("Intercept", "Change in Mean Wet Season Precipitation (cm)", "Change in Mean Dry Season Precipitation (cm)", "Change in Mean Breeding Season Temperature (C)",
                                                              "Change in Mangrove Cover (%)", "Change in Emergent Wetland Cover (%)", "Change in Woody Wetland Cover (%)",
                                                              "Change in Anthropogenic Cover (%)", "Species Richness", "Duration of Survey (Years)"),
          dv.labels = "Beta Diversity")

mod2 <- lmer(data = seg.df, beta.wet.jac.log ~ p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.human + Duration + SR.wet + (1|Route), REML = F)
summary(mod2)
car::Anova(mod2)
tab_model(mod2)

mod2 <- glmer(data = seg.df, beta.wet.jac ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s  + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s + SR.wet.s + (1|Route), family = (link = "log"))
summary(mod2)
car::Anova(mod2)


tab_model(mod1, mod2, show.ci = 0.95, title = NULL, pred.labels = c("Intercept", "Change in Wet Season Precipitation (cm)", "Change in Dry Season Precipitation (cm)", "Change in Mean Breeding Season Temperature (C)",
                                                                              "Change in Mangrove Cover (%)", "Change in Emergent Wetland Cover (%)", "Change in Woody Wetland Cover (%)",
                                                                              "Change in Anthropogenic Cover (%)", "Duration of Survey (Years)", "Species Richness", "Wetland Bird Species Richness"),
                          dv.labels = c("Total Bird Community Beta Diversity", "Wetland Bird Community Beta Diversity"), linebreak = T
          )

#Turnover Component Wetland and Entire Community
# mod3 <- lmer(data = seg.df, beta50.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR50 + Duration + (1|rteno), REML = F)
# summary(mod3)
# Anova(mod3)
# tab_model(mod3)
# 
# mod4 <- lmer(data = seg.df, beta.wet.turn ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR.wet + Duration + (1|rteno), REML = F)
# summary(mod4)
# Anova(mod4)
# tab_model(mod4)
# 
# #Nestedness Compoentn Wetland and Entire Community
# mod5 <- lmer(data = seg.df, beta50.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR50 + Duration + (1|rteno), REML = F)
# summary(mod5)
# Anova(mod5)
# tab_model(mod5)
# 
# mod6 <- lmer(data = seg.df, beta.wet.nest ~ p.anom.wet + mean.anom + diff.from.first.man + diff.from.first.ew + diff.from.first.ww + diff.from.first.ur + SR.wet + Duration + (1|rteno), REML = F)
# summary(mod6)
# Anova(mod6)
# tab_model(mod6)

#Change SR 
mod7 <- lmer(data = seg.df, alpha50.pchange ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s + SR50.s + (1|Route), REML = F)
summary(mod7)
car::Anova(mod7)
tab_model(mod7)

mod8 <- lmer(data = seg.df, alphawet.pchange ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s + SR.wet.s + (1|Route), REML = F)
summary(mod8)
car::Anova(mod8)
tab_model(mod8)

tab <- tab_model(mod1, mod7, mod2, mod8, show.ci = 0.95, title = NULL, pred.labels = c("Intercept", "Change in Mean Wet Season Precipitation (cm)", 
                                                                                       "Change in Mean Dry Season Precipitation (cm)", "Change in Mean Breeding Season Temperature (C)",
                                                                    "Change in Mangrove Cover (%)", "Change in Emergent Wetland Cover (%)", "Change in Woody Wetland Cover (%)",
                                                                    "Change in Anthropogenic Cover (%)", "Duration of Survey (Years)", "Species Richness", "Wetland Bird Species Richness"),
          dv.labels = c("Total Bird Community Beta Diversity", "Relative Change Alpha Diversity", "Wetland Bird Community Beta Diversity", "Relative Change Wetland Bird Alpha Diversity"), linebreak = F,
          CSS = list(css.modelcolumn1 = 'background-color: #f0f0f0;', 
                     css.modelcolumn3 = 'background-color: #f0f0f0;',
                     css.lasttablerow = 'border-bottom: 4px solid black;'), use.viewer = T)
          #file = here::here("Figures/Figures_Diversity_Manuscript/TabMod.pdf"))


tab
########################################################################################################################################################################################################

#Multiple Regression Models @ Route level
rt.df$p.anom.s <- scale(rt.df$p.anom) 
rt.df$p.anom.wet.s <- scale(rt.df$p.anom.wet)
rt.df$p.anom.dry.s <- scale(rt.df$p.anom.dry)
rt.df$p.anom.f.s <- scale(rt.df$p.anom.f)
rt.df$p.anom.s.s <- scale(rt.df$p.anom.s)
rt.df$p.anom.sp.s <- scale(rt.df$p.anom.sp)
rt.df$p.anom.w.s <- scale(rt.df$p.anom.w)
rt.df$mean.anom.bird.s <- scale(rt.df$mean.anom.bird)
rt.df$diff.from.first.man.s <- scale(rt.df$diff.from.first.man)
rt.df$diff.from.first.ew.s <- scale(rt.df$diff.from.first.ew)
rt.df$diff.from.first.human.s <- scale(rt.df$diff.from.first.human)
rt.df$diff.from.first.wwnm.s <- scale(rt.df$diff.from.first.wwnm)
rt.df$diff.from.first.ww.s <- scale(rt.df$diff.from.first.ww)
rt.df$diff.from.first.nat.s <- scale(rt.df$diff.from.first.nat)
rt.df$Duration.s <- scale(rt.df$Duration)
rt.df$SR50.s <- scale(rt.df$SR50)
rt.df$SR.wet.s <- scale(rt.df$SR.wet)
rt.df$beta50.jac.log <- log(rt.df$beta50.jac)

#Jaccard Index Models
mod9 <- lm(data = rt.df, beta50.jac.log ~ p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.human + Duration + SR50)
mod9b <- betareg::betareg(data = rt.df, beta50.jac ~ p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.human + Duration + SR50)

summary(mod9)
summary(mod9b)
tab_model(mod9)
tab_model(mod9b)
car::avPlots(mod9)
car::vif(mod9)

mod10 <- lm(data = rt.df, beta.wet.jac ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s + SR.wet)
summary(mod10)

#Turnover Index Models
# mod11 <- lm(data = rt.df, beta50.turn ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR50 + Duration)
# summary(mod11)

# mod12 <- lm(data = rt.df, beta.wet.turn ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR.wet + Duration)
# summary(mod12)

#Nestedness Index Models 
# mod13 <- lm(data = rt.df, beta50.nest ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR50 + Duration)
# summary(mod13)
# 
# mod14 <- lm(data = rt.df, beta.wet.nest ~ p.anom + mean.anom.bird + diff.from.first.man + diff.from.first.ew + diff.from.first.wwnm + diff.from.first.ur + SR.wet + Duration)
# summary(mod14)
# 
#Change SR
mod15 <- lm(data = rt.df, alpha50.pchange ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s)
summary(mod15)
car::avPlots(mod15)

mod16 <- lm(data = rt.df, alphawet.pchange ~ p.anom.wet.s + p.anom.dry.s + mean.anom.bird.s + diff.from.first.man.s + diff.from.first.ew.s + diff.from.first.wwnm.s + diff.from.first.human.s + Duration.s)
summary(mod16)
car::avPlots(mod16)

tab2 <- tab_model(mod9, mod15, mod10, mod16, show.ci = 0.95, title = NULL, pred.labels = c("Intercept", "Change in Mean Wet Season Precipitation (cm)", "Change in Mean Dry Season Preciptation (cm)", "Change in Mean Breeding Season Temperature (C)",
                                                                                       "Change in Mangrove Cover (%)", "Change in Emergent Wetland Cover (%)", "Change in Woody Wetland Cover (%)",
                                                                                       "Change in Anthropogenic Cover (%)", "Duration of Survey (Years)", "Species Richness", "Wetland Bird Species Richness"),
                 dv.labels = c("Total Bird Community Beta Diversity", "Relative Change Alpha Diversity", "Wetland Bird Community Beta Diversity", "Relative Change Wetland Bird Alpha Diversity"), linebreak = F,
                 CSS = list(css.modelcolumn1 = 'background-color: #f0f0f0;', 
                            css.modelcolumn3 = 'background-color: #f0f0f0;',
                            css.lasttablerow = 'border-bottom: 4px solid black;'), use.viewer = T)

tab2

########################################################################################################################################################################################################
bbs_clim_rt1 <- bbs_total_seg %>% group_by(rteno.x) %>% arrange(Year) %>% filter(Year <= 1999) %>% 
  select(c("rteno.x", "Year", "Latitude", "SR50", "SR.wet")) %>% mutate(Year = as.character(Year)) %>% summarize_if(is.numeric, mean) %>% ungroup() 
bbs_clim_rt2 <- bbs_total_seg %>% group_by(rteno.x) %>% arrange(Year) %>% filter(Year >= 2000) %>% 
  select(c("rteno.x", "Year", "Latitude", "SR50", "SR.wet")) %>% mutate(Year = as.character(Year)) %>% summarize_if(is.numeric, mean) %>% ungroup() 

#bbs_clim_rt2 <- bbs_clim_rt2[bbs_clim_rt2$rteno.x %in% bbs_clim_rt1$rteno.x, ]
#bbs_clim_rt1 <- bbs_clim_rt1[bbs_clim_rt1$rteno.x %in% bbs_clim_rt2$rteno.x, ]

bbs_clim_rt1$indicator <- "1980 - 1999"
bbs_clim_rt2$indicator <- "2000 - 2017"

bbs_clim_f <- rbind(bbs_clim_rt1, bbs_clim_rt2)

bbs_clim_f <- reshape2::melt(data = bbs_clim_f, id = c("indicator", "Latitude"), measure.vars = c("SR50", "SR.wet"))

bbs_clim_f$variable_f <- factor(bbs_clim_f$variable, levels = c("SR50", "SR.wet"))
levels(bbs_clim_f$variable_f) <- c(SR50 = expression(paste(alpha, "-diversity")),
                                   SR.wet = expression(paste(alpha, "-diversity"["Wetland Birds"])))


sr.peak <- ggplot() + 
  geom_smooth(data = bbs_clim_f, aes(x = Latitude, y = value, group = indicator, colour = factor(indicator), fill = factor(indicator))) +
  #geom_smooth(data = bbs_clim_f, aes(x = Latitude, y = SR.wet, colour = factor(indicator), fill = factor(indicator))) +
  #geom_smooth(data = bbs_clim_rt2, aes(x = Latitude, y = SR50), color = "black") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, family = "serif"),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  #scale_color_manual(values = c("#00204DFF", "#B9AC70FF")) +
  #scale_fill_manual(values = c("#00204DFF", "#B9AC70FF")) +
  scale_color_viridis(discrete = T, option = "E") + 
  scale_fill_viridis(discrete = T, option = "E") +
  labs(colour = "Time Series", fill = "Time Series") +
  xlab("Latitude") +
  ylab("Detection-corrected Species Richness") +
  facet_wrap(~variable_f, labeller = label_parsed)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/StartEnd_SR_SRwet.tiff"), plot = sr.peak,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Alpha and Beta Change Raw vs Occ
div.cora <- bbs_total_seg[, c("Site_div.x", "SR_MCMC", "SR50", "Year", "site", "rteno.x")]
div.cora <- reshape2::melt(data = div.cora, id = c("Year", "site", "rteno.x"), measure.vars = c("Site_div.x", "SR_MCMC", "SR50"))

div.corb <- bbs_total_seg[, c("beta.jac", "beta.mcmc", "beta50.jac", "Year", "site", "rteno.x")]
div.corb <- reshape2::melt(data = div.corb, id = c("Year", "site", "rteno.x"), measure.vars = c("beta.jac", "beta50.jac", "beta.mcmc"))


alpha <- ggplot() + 
  geom_smooth(data = div.cora, method = "loess", aes(x = Year, y = value, colour = factor(variable,
                                                                                          labels = c("Naive SR",
                                                                                                     "SR MCMC",
                                                                                                     "SR 50%")), 
                                                     fill = factor(variable, labels = c("Naive SR",
                                                                                        "SR MCMC",
                                                                                        "SR 50%")))) + #group = interaction(rteno.x, variable)
  theme_bw() +
  theme(axis.text = element_text(size = 12, family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif")) +
  xlab("Year") +
  ylab("Species Richness") +
  labs(colour = c("Species Estimation Method"), fill = c("Species Estimation Method")) +
  scale_color_viridis(discrete = T, option = "D") +
  scale_fill_viridis(discrete = T, option = "D") +
  facet_wrap(~variable, nrow = 1, labeller = )

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Occ_Plot_SR.tiff"), plot = alpha,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

  
beta.occ.plot <- ggplot() + 
  geom_smooth(data = div.corb, method = "loess", aes(x = Year, y = value, colour = factor(variable,
                                                                                          labels = c("Naive Beta",
                                                                                                     "Beta 50%",
                                                                                                     "Beta MCMC")), 
                                                     fill = factor(variable, labels = c("Naive Beta",
                                                                                        "Beta 50%",
                                                                                        "Beta MCMC")))) + #group = interaction(rteno.x, variable)
  theme_bw() +
  theme(axis.text = element_text(size = 12, family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12, family = "serif"),
        legend.title = element_text(size = 12, family = "serif")) +
  xlab("Year") +
  ylab("Beta") +
  labs(colour = c("Species Estimation Method"), fill = c("Species Estimation Method")) +
  scale_color_viridis(discrete = T, option = "E") +
  scale_fill_viridis(discrete = T, option = "E") +
  facet_wrap(~variable, nrow = 1, labeller = )

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Occ_Plot_Beta.tiff"), plot = beta.occ.plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

big.plot <- plot_grid(alpha, beta.occ.plot, nrow = 2, align = "v")
ggsave(here::here("Figures/Figures_Diversity_Manuscript/Occ_Plot_Total.tiff"), plot = big.plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


########################################################################################################################################################################################################

#-Wetland, -WW, +Ag  
mod22 <- lm(data = rt.df, SR50 ~ Latitude + I(Latitude^2))
summary(mod22)

#Plot for mod22
set.seed(100)
prd1 <- data.frame(Latitude = seq(from  = range(rt.df$Latitude)[1], to = range(rt.df$Latitude)[2], length.out = 100))
err1 <- predict(mod22, newdata = prd1, se.fit = T)

prd1$lci <- err1$fit - 1.96 * err1$se.fit
prd1$fit <- err1$fit
prd1$uci <- err1$fit + 1.96 * err1$se.fit


#model for latitude and wetland bird richness
mod23 <- lm(data = rt.df, SR.wet ~ Latitude + I(Latitude^2))
summary(mod23)

set.seed(100)
prd2 <- data.frame(Latitude = seq(from  = range(rt.df$Latitude)[1], to = range(rt.df$Latitude)[2], length.out = 100))
err2 <- predict(mod23, newdata = prd2, se.fit = T)

prd2$lci <- err2$fit - 1.96 * err2$se.fit
prd2$fit <- err2$fit
prd2$uci <- err2$fit + 1.96 * err2$se.fit

#Plots for latitude vs SR
lat.plot <- ggplot(prd1, aes(x = Latitude, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#A09877FF", fill = "#A09877FF") +
  geom_point(data = rt.df, aes(x = Latitude, y = SR50)) +
  xlab(NULL) + #expression(paste("Latitude ( ", degree, " )"
  ylab('Bird Species Richness') +
  theme(axis.title.y = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif")) +
        #panel.background = element_blank(), panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(24, 32, 1))
# labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
#                    "Intercept =",signif(mod21$coef[[1]],5 ),
#                    " Slope =",signif(mod21$coef[[3]], 5),
#                    " P =",signif(summary(mod21)$coef[3,4], 5))) +

lat.plot


lat.plot1 <- ggplot(prd2, aes(x = Latitude, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#00204DFF", fill = "#00204DFF") +
  geom_point(data = rt.df, aes(x = Latitude, y = SR.wet)) +
  xlab(expression(paste("Latitude ( ", degree, " )"))) +
  ylab('Wetland Bird Species Richness') +
  theme(axis.title.x = element_text(size = 12, family = "serif"), axis.title.y = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif")) +
  scale_x_continuous(breaks = seq(24, 32, 1))
# labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
#                    "Intercept =",signif(mod21$coef[[1]],5 ),
#                    " Slope =",signif(mod21$coef[[3]], 5),
#                    " P =",signif(summary(mod21)$coef[3,4], 5))) +

lat.plot1

lat.plot.f <- plot_grid(lat.plot, lat.plot1, labels = c('A', 'B'), align = "v", label_x = 0.08, label_y = 0.98, nrow = 2, label_size = 14, label_fontfamily = "serif")

ggsave(here::here("Figures/Figures_Diversity_Manuscript/SR_Lat.tiff"), plot = lat.plot.f,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


#-WW, +EW, -WWNM, +Ag
mod23 <- lm(data = rt.df, SR.wet ~ pct.man)
summary(mod23)
Anova(mod23)
tab_model(mod23)

########################################################################################################################################################################################################
#CMRL - Beta div plots
rt.df$change.cmrl.s <- scale(rt.df$change.cmrl, center = T, scale = T)
rt.df$change.cmrl.km.s <- scale(rt.df$change.cmrl.km, center = T, scale = T)

mod17 <- lm(data = rt.df, beta50.jac ~ change.cmrl.s + I(change.cmrl.s^2))
summary(mod17)
tab_model(mod17)
#mod17 <- lm(data = rt.df, beta50.jac ~ change.cmrl.s + I(change.cmrl.s^2))
#summary(mod17)

mod18 <- lm(data = rt.df, beta50.turn ~ change.cmrl.s + I(change.cmrl.s^2))
summary(mod18)
tab_model(mod18)

mod19 <- lm(data = rt.df, beta50.nest ~ change.cmrl.s + I(change.cmrl.s^2))
summary(mod19)
tab_model(mod19)

#Beta Jac CMRL
set.seed(100)
prd2 <- data.frame(change.cmrl.km = seq(from  = range(rt.df$change.cmrl.km)[1], to = range(rt.df$change.cmrl.km)[2], length.out = 100))
err2 <- predict(mod17, newdata = prd2, se.fit = T)

prd2$lci <- err2$fit - 1.96 * err2$se.fit
prd2$fit <- err2$fit
prd2$uci <- err2$fit + 1.96 * err2$se.fit

cmrl.plot1 <- ggplot(prd2, aes(x = change.cmrl.km, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#00204DFF", fill = "#00204DFF") +
  geom_point(data = rt.df, aes(x = change.cmrl.km, y = beta50.jac)) +
  xlab("Change in CMRL (km)") +
  ylab(expression(beta["Jac"])) +
  # labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
  #                    "Intercept =",signif(mod21$coef[[1]],5 ),
  #                    " Slope =",signif(mod21$coef[[3]], 5),
  #                    " P =",signif(summary(mod21)$coef[3,4], 5))) +
  #scale_x_continuous(breaks = seq(-400, 300, 100)) +
  theme(axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"))

cmrl.plot1

ggsave(here::here("Figures/Figures_Diversity_Manuscript/CMRL_Beta.tiff"), plot = cmrl.plot1,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Beta Turn CMRL
set.seed(100)
prd3 <- data.frame(change.cmrl.km = seq(from  = range(rt.df$change.cmrl.km)[1], to = range(rt.df$change.cmrl.km)[2], length.out = 100))
err3 <- predict(mod18, newdata = prd3, se.fit = T)

prd3$lci <- err3$fit - 1.96 * err3$se.fit
prd3$fit <- err3$fit
prd3$uci <- err3$fit + 1.96 * err3$se.fit

cmrl.plot2 <- ggplot(prd3, aes(x = change.cmrl.km, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#777776FF", fill = "#777776FF") +
  geom_point(data = rt.df, aes(x = change.cmrl.km, y = beta50.turn)) +
  xlab("Change in CMRL (km)") +
  ylab(expression(beta["Turn"])) +
  # labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
  #                    "Intercept =",signif(mod21$coef[[1]],5 ),
  #                    " Slope =",signif(mod21$coef[[3]], 5),
  #                    " P =",signif(summary(mod21)$coef[3,4], 5))) +
  #scale_x_continuous(breaks = seq(-400, 300, 100)) +
  theme(axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"))

cmrl.plot2

#Beta Nest CMRL
set.seed(100)
prd4 <- data.frame(change.cmrl.km = seq(from  = range(rt.df$change.cmrl.km)[1], to = range(rt.df$change.cmrl.km)[2], length.out = 100))
err4 <- predict(mod19, newdata = prd4, se.fit = T)

prd4$lci <- err4$fit - 1.96 * err4$se.fit
prd4$fit <- err4$fit
prd4$uci <- err4$fit + 1.96 * err4$se.fit

cmrl.plot3 <- ggplot(prd4, aes(x = change.cmrl.km, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#D7C463FF", fill = "#D7C463FF") +
  geom_point(data = rt.df, aes(x = change.cmrl.km, y = beta50.nest)) +
  xlab("Change in CMRL (km)") +
  ylab(expression(beta["Nest"])) +
  # labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
  #                    "Intercept =",signif(mod21$coef[[1]],5 ),
  #                    " Slope =",signif(mod21$coef[[3]], 5),
  #                    " P =",signif(summary(mod21)$coef[3,4], 5))) +
  #scale_x_continuous(breaks = seq(-400, 300, 100)) +
  theme(axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"))

cmrl.plot3





mod18 <- lmer(data = seg.df, beta50.jac ~ change.cmrl + (1|rteno))
summary(mod18)
car::Anova(mod18)

plot(seg.df$change.cmrl, seg.df$beta50.jac)

mod18 <- lm(data = rt.df, beta50.turn ~ change.cmrl)
summary(mod18)

mod19 <- lm(data = rt.df, beta50.nest ~ change.cmrl)
summary(mod19)

mod24 <- lm(data = rt.df, alpha50.pchange ~ change.cmrl)
summary(mod24)

mod20 <- lmer(data = seg.df, change.cmrl.km ~ mean.anom.sp + I(mean.anom.sp^2) + (1|Route), REML = F)
summary(mod20)
car::Anova(mod20)
tab_model(mod20)

# rt.df$mean.anom.sp.c <- scale(rt.df$mean.anom.sp)
# rt.df$p.anom.sp.c <- scale(rt.df$p.anom.sp)
# rt.df$diff.from.first.man.c <- scale(rt.df$diff.from.first.man)
mod21 <- lm(data = rt.df, change.cmrl.km ~ mean.anom.sp + I(mean.anom.sp^2))
summary(mod21)

set.seed(100)
prd <- data.frame(mean.anom.sp = seq(from  = range(rt.df$mean.anom.sp)[1], to = range(rt.df$mean.anom.sp)[2], length.out = 100))
err <- predict(mod21, newdata = prd, se.fit = T)

prd$lci <- err$fit - 1.96 * err$se.fit
prd$fit <- err$fit
prd$uci <- err$fit + 1.96 * err$se.fit

cmrl.plot <- ggplot(prd, aes(x = mean.anom.sp, y = fit)) +
  theme_bw() +
  geom_line() +
  geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity", color = "#A09877FF", fill = "#A09877FF") +
  geom_point(data = rt.df, aes(x = mean.anom.sp, y = change.cmrl.km)) +
  xlab(expression(paste("Change in Mean Spring Temp ( ", degree ~ C, " )"))) +
  ylab('Change in CMRL (km)') +
  # labs(title = paste("R2 = ",signif(summary(mod21)$r.squared, 5),
  #                    "Intercept =",signif(mod21$coef[[1]],5 ),
  #                    " Slope =",signif(mod21$coef[[3]], 5),
  #                    " P =",signif(summary(mod21)$coef[3,4], 5))) +
  scale_x_continuous(breaks = seq(-0.5, 2.5, 0.5)) +
  theme(axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"))
  
cmrl.plot

ggsave(here::here("Figures/Figures_Diversity_Manuscript/CMRL_SpMeanTemp.tiff"), plot = cmrl.plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

cmrl.plots <- plot_grid(cmrl.plot, cmrl.plot1, cmrl.plot2, cmrl.plot3, align = "hv", nrow = 2,
                        labels = c("A", "B", "C", "D"), label_size = 12, label_fontfamily = "serif")

ggsave(here::here("Figures/Figures_Diversity_Manuscript/CMRL_plot_tot.tiff"), plot = cmrl.plots,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


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
gom <- subset(states, region %in% c("texas", "florida", "alabama", "arkansas", "mississippi", "louisiana", "georgia", "south carolina"))

#Beta Turn
map_gom_turn <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.turn), size = 3)

map_gom_turn <- map_gom_turn + scale_color_viridis(name = expression(paste(beta, "-diversity"[italic("Turn")])), option = "plasma") + xlab("Longitude") + ylab("Latitude") + 
  theme(legend.position = c(0.12, .77), legend.title = element_text(size = 12, family = "serif"), legend.text = element_text(size = 12, family = "serif"), legend.background = element_rect(fill = "transparent"))  

# ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_Beta_Turn.tiff"), plot = map_gom_turn,
#        device = "tiff", width = 8, height = 5, units = "in", dpi = 600)
# 
# dev.off()

#Beta Nest
map_gom_nest <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.nest), size = 3)

map_gom_nest <- map_gom_nest + scale_color_viridis(option = 'magma', name = expression(paste(beta, "-diversity"[italic("Nest")])))+ xlab("Longitude") + ylab("Latitude") + theme(legend.position = c(0.12, .77), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.background = element_rect(fill = "transparent"))  

# ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_Beta_Nest.tiff"), plot = map_gom_nest,
#        device = "tiff", width = 8, height = 5, units = "in", dpi = 600)
# 
# dev.off()

#Map of Alpha Div
map_gom_alpha <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = SR50), size = 3)

map_gom_alpha <- map_gom_alpha + scale_color_viridis(name = expression(paste(alpha, "-diversity")))+ xlab("Longitude") + ylab("Latitude") + theme(legend.position = c(0.12, .77), legend.title = element_text(size = 14), legend.text = element_text(size = 12), legend.background = element_rect(fill = "transparent"))  

# ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_SR.tiff"), plot = map_gom_alpha,
#        device = "tiff", width = 8, height = 5, units = "in", dpi = 600)
# 
# dev.off()

###########################################################
########################Most important maps################
###########################################################

#Beta Jaccard
map_gom_jac <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) + #coord_map(projection = "albers", lat0 = 26, lat1 = 30)+
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta50.jac), size = 3)

map_gom_jac <- map_gom_jac + scale_color_viridis(name = expression(paste(beta, "-diversity"[italic("Jac")])), option = "E") + xlab("Longitude") + ylab("Latitude") + 
  theme_bw() + theme(legend.position = c(0,1), legend.title = element_text(size = 12, family = "serif"), legend.text = element_text(size = 12, family = "serif"), legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.title = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif"),
        legend.justification = c(0,1))   

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_Beta_Jac.tiff"), plot = map_gom_jac,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Beta Wetland 
map_gom_wet <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = beta.wet.jac), size = 3)

map_gom_wet <- map_gom_wet + scale_color_viridis(option = 'magma', name = expression(paste("Wetland Bird ", beta[italic("Jac")])))+ xlab("Longitude") + ylab("Latitude") + 
        theme_bw() + theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.title = element_text(size = 12, family = "serif"), legend.text = element_text(size = 12, family = "serif"), legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.title = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif"))   

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_Beta_Wet.tiff"), plot = map_gom_wet,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Map of change in Alpha 
map_change_alpha <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = alpha.pchange), size = 3)

map_change_alpha <- map_change_alpha + scale_color_viridis(name = expression(paste(alpha, "-diversity ", "(% ", Delta, ")")), option = "E")+ xlab("Longitude") + ylab("Latitude") + 
  theme_bw() + theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.title = element_text(size = 12, family = "serif"), legend.text = element_text(size = 12, family = "serif"), legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), axis.title = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif"))  

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_SR_change.tiff"), plot = map_change_alpha,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Map in change in alpha Wet
map_change_alpha_wet <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = rt.df, aes(x = Longitude, y = Latitude, color = alphawet.pchange), size = 3)

map_change_alpha_wet <- map_change_alpha_wet + scale_color_viridis(name = expression(paste("Wetland Bird ", alpha, " (% ", Delta, ")")), option = "magma")+ xlab("Longitude") + ylab("Latitude") + 
  theme_bw() + theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.title = element_text(size = 12, family = "serif"), legend.text = element_text(size = 12, family = "serif"), legend.background = element_rect(fill = "transparent"),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.title = element_text(size = 12, family = "serif"), axis.text = element_text(size = 12, family = "serif"))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/Map_SRwet_change.tiff"), plot = map_change_alpha_wet,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)


master_map <- plot_grid(map_gom_jac, map_gom_wet, map_change_alpha, map_change_alpha_wet, labels = c("A", "B", "C", "D"), align = "hv", label_x = 0.05, label_y = 1, axis = "rlbt", label_fontfamily = "serif")

ggsave(here::here("Figures/Figures_Diversity_Manuscript/MasterMap.tiff"), plot = master_map,
       device = "tiff", width = 15, height = 10, units = "in", dpi = 600)



#Map of Sites 
seg.df$BCR_name <- NA
seg.df$BCR_name[seg.df$BCR == 26] <- "Mississippi Alluvial Valley"
seg.df$BCR_name[seg.df$BCR == 37] <- "Gulf Coast Prairie"
seg.df$BCR_name[seg.df$BCR == 27] <- "Southeastern Coastal Plain"
seg.df$BCR_name[seg.df$BCR == 31] <- "Peninsular Florida"


map <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = seg.df, aes(x = Longitude, y = Latitude, colour = factor(BCR_name))) +#, labels = c("Gulf Coast Prairie",
                                                                                           #"Mississippi Alluvial Valley", 
                                                                                           #"Peninsular Florida",
                                                                                           #"Southeastern Coastal Plain"))), size = 3) +
  scale_color_viridis(option = "E", discrete = T) + theme_bw() + xlab("Longitude") + ylab("Latitude") + labs(color = "Bird Conservation Region") + theme(legend.title = element_text(size = 11, face = "bold", family = "serif"),
                                                                                                                                                         legend.text = element_text(size = 11, family = "serif"),
                                                                                                                                                         axis.title = element_text(size = 14, family = "serif"),
                                                                                                                                                         axis.text = element_text(size = 14, family = "serif"),
                                                                                                                                                         legend.background = element_rect("transparent"),
                                                                                                                                                         legend.key = element_rect("transparent"),
                                                                                                                                                         legend.position = c(0, 1),
                                                                                                                                                         legend.justification = c(0,1),
                                                                                                                                                         legend.key.size = unit(0.25, "cm")
                                                                                                                                                         )

ggsave(here::here("Figures/Figures_Diversity_Manuscript/BasicMap.tiff"), plot = map,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

# map_gom_cmrl <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(rt.df$Longitude) - 1.5, max(rt.df$Longitude)), ylim = c(min(rt.df$Latitude), max(rt.df$Latitude) + 2.5), ratio = 1.2) +
#   geom_point(data = seg.df, aes(x = Longitude, y = Latitude, color = change.cmrl.km), size = 3)
# 
# map_gom_cmrl + scale_color_viridis(name = expression(paste(beta, "-diversity"))) +  theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  


#Plotting Time series of Beta-diversity Seg Level Jaccard
gghist_beta <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$beta50jac.slope)) + 
  #geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  geom_density(alpha = 0.4, fill = "#444F6BFF") + #042333b2(Magma) #450A69FF(Viridis) #F98C0AFF (plasma)
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", beta[italic("Jac")])), y = "# of Sites") +
  theme_cowplot(font_size = 15, line_size = 1.2, font_family = "serif") + 
  coord_flip()

#site_data_merge$Year <- site_data_merge$count_yr + 1900

plot_beta <- (ggplot(bbs_total_seg, aes(x = Year, y = beta50.jac, group = BCR,
                                        colour = factor(BCR_name),
                                        fill = factor(BCR_name))) +
                                    # colour = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                    #                                 "Southeastern Coastal Plain", 
                                    #                                 "Peninsular Florida", 
                                    #                                 "Gulf Coastal Prairie")), 
                                    # fill = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                    #                               "Southeastern Coastal Plain", 
                                    #                               "Peninsular Florida", 
                                    #                               "Gulf Coastal Prairie")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = loess, se = T, aes(x = Year, y = beta50.jac, group = BCR,
                                                        colour = factor(BCR_name),
                                                        fill = factor(BCR_name))) +
                                                        # colour = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                        #                                 "Southeastern Coastal Plain", 
                                                        #                                 "Peninsular Florida", 
                                                        #                                 "Gulf Coastal Prairie")),
                                                        # fill = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                        #                               "Southeastern Coastal Plain", 
                                                        #                               "Peninsular Florida", 
                                                        #                               "Gulf Coastal Prairie")))) +
                #scale_color_brewer(palette = )
                xlab("Year") +
                ylab(expression(paste(beta[italic("Jac")]))) +
                labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
                #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 15, family = "serif"),
                      axis.text.y = element_text(size = 15, family = "serif"),
                      axis.title.x = element_text(vjust = -1, size = 15, family = "serif"),
                      axis.title.y = element_text(vjust = 1.5, size = 15, family = "serif"),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=15, family = "serif"))) 

plot_beta <- plot_beta + scale_color_viridis(discrete = T, option = "E") + scale_fill_viridis(discrete = T, option = "E")

betaseg_jac <- ggdraw() + 
  draw_plot(plot_beta + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.65, y = 0.05, width = .3, height = .70, scale = 1) 

ggsave(here::here("Figures/Figures_Diversity_Manuscript/BetaSegJac.tiff"), plot = betaseg_jac,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Beta Seg Jac Wet
gghist_beta <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$betawetjac.slope)) + 
  #geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
  geom_density(alpha = 0.4, fill = "#042333b2") + #042333b2(Magma) #450A69FF(Viridis) #F98C0AFF (plasma)
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", beta[italic("Jac Wetland")]), width = 10), y = "# of Sites") +
  theme_cowplot(font_size = 15, line_size = 1.2, font_family = "serif") +
  coord_flip()

#site_data_merge$Year <- site_data_merge$count_yr + 1900

plot_beta_wet <- (ggplot(bbs_total_seg, aes(x = Year, y = beta.wet.jac, group = BCR,
                                            colour = factor(BCR_name),
                                            fill = factor(BCR_name))) +
                                        # colour = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                        #                                 "Southeastern Coastal Plain", 
                                        #                                 "Peninsular Florida", 
                                        #                                 "Gulf Coastal Prairie")),
                                        # fill = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                        #                               "Southeastern Coastal Plain", 
                                        #                               "Peninsular Florida", 
                                        #                               "Gulf Coastal Prairie")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = loess, se = T, aes(x = Year, y = beta.wet.jac, group = BCR, 
                                                        colour = factor(BCR_name),
                                                        fill = factor(BCR_name))) +
                                                        # colour = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                        #                                 "Southeastern Coastal Plain", 
                                                        #                                 "Peninsular Florida", 
                                                        #                                 "Gulf Coastal Prairie")),
                                                        # fill = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                        #                               "Southeastern Coastal Plain", 
                                                        #                               "Peninsular Florida", 
                                                        #                               "Gulf Coastal Prairie")))) +
                #scale_color_brewer(palette = )
                xlab("Year") +
                ylab(expression(paste("Wetland Bird ", beta[italic("Jac")]))) +
                labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
                #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 15, family = "serif"),
                      axis.text.y = element_text(size = 15, family = "serif"),
                      axis.title.x = element_text(vjust = -1, size = 15, family = "serif"),
                      axis.title.y = element_text(vjust = 1.5, size = 15, family = "serif"),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=15, family = "serif"))) 

plot_beta_wet <- plot_beta_wet + scale_color_viridis(option = "magma", discrete = T) + scale_fill_viridis(option = "magma", discrete = T)

betasegwet_jac <- ggdraw() + 
  draw_plot(plot_beta_wet + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.65, y = 0.05, width = .3, height = .70, scale = 1) 

ggsave(here::here("Figures/Figures_Diversity_Manuscript/BetaSegWetJac.tiff"), plot = betasegwet_jac,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

# #Beta Rt Jac
# #Plotting Time series of Beta-diversity Rt & Seg Level
# gghist_beta <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$beta50jac.slope)) + 
#   #geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
#   geom_density(alpha = 0.4, fill = "#450A69FF") + #042333b2(Magma) #450A69FF(Viridis) #F98C0AFF (plasma)
#   labs(title = "") +
#   labs(x = expression(paste("Slopes of ", beta, "-diversity"[italic("Nest")])), y = "# of Sites") +
#   theme_cowplot(font_size = 15, line_size = 1.2, font_family = "serif") +
#   coord_flip()
# 
# #site_data_merge$Year <- site_data_merge$count_yr + 1900
# 
# plot_beta <- (ggplot(bbs_total_seg, aes(x = Year, y = beta50.jac, group = BCR, 
#                                         colour = factor(BCR, 
#                                                         labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
#                 #geom_point(size = 3) +
#                 #geom_line()+
#                 geom_smooth(method = loess, se = T, aes(x = Year, y = beta50.jac, group = BCR, 
#                                                         colour = factor(BCR,
#                                                                         labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
#                 #scale_color_brewer(palette = )
#                 xlab("Year") +
#                 ylab(expression(paste(beta[italic("Jaccard")]))) +
#                 labs(colour = "Bird Conservation Region") +
#                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
#                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
#                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
#                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#                 theme_bw() +
#                 theme(axis.line = element_line(colour = "black", size =1.2),
#                       axis.text.x = element_text(size = 15),
#                       axis.text.y = element_text(size = 15),
#                       axis.title.x = element_text(vjust = -1, size = 15),
#                       axis.title.y = element_text(vjust = 1.5, size = 15),
#                       axis.ticks = element_line(size = 1.2),
#                       panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(),
#                       panel.border = element_blank(),
#                       panel.background = element_blank(),
#                       plot.margin = unit(c(1,1,2,2), "lines"),
#                       text = element_text(size=15))) 
# 
# plot_beta <- plot_beta + scale_color_viridis(option = "plasma", discrete = T)
# 
# ggdraw() + 
#   draw_plot(plot_beta + theme(legend.justification = "top"), 
#             x = 0, y = 0, width = .9, height = 1) +
#   draw_plot(gghist_beta, x = 0.72, y = .030, width = .2, height = .75, scale = 1) 
# 
# 
# #Beta Rt Jac Wet
# #Plotting Time series of Beta-diversity Rt & Seg Level
# gghist_beta <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$beta50jac.slope)) + 
#   #geom_histogram(col = "black", fill = "black", bins = 10, binwidth = NULL) + 
#   geom_density(alpha = 0.4, fill = "#450A69FF") + #042333b2(Magma) #450A69FF(Viridis) #F98C0AFF (plasma)
#   labs(title = "") +
#   labs(x = expression(paste("Slopes of ", beta, "-diversity"[italic("Nest")])), y = "# of Sites") +
#   theme_cowplot(font_size = 15, line_size = 1.2) +
#   coord_flip()
# 
# #site_data_merge$Year <- site_data_merge$count_yr + 1900
# 
# plot_beta <- (ggplot(bbs_total_seg, aes(x = Year, y = beta50.jac, group = BCR, 
#                                         colour = factor(BCR, 
#                                                         labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
#                 #geom_point(size = 3) +
#                 #geom_line()+
#                 geom_smooth(method = loess, se = T, aes(x = Year, y = beta50.jac, group = BCR, 
#                                                         colour = factor(BCR,
#                                                                         labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
#                 #scale_color_brewer(palette = )
#                 xlab("Year") +
#                 ylab(expression(paste(beta[italic("Jaccard")]))) +
#                 labs(colour = "Bird Conservation Region") +
#                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
#                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
#                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
#                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
#                 theme_bw() +
#                 theme(axis.line = element_line(colour = "black", size =1.2),
#                       axis.text.x = element_text(size = 15),
#                       axis.text.y = element_text(size = 15),
#                       axis.title.x = element_text(vjust = -1, size = 15),
#                       axis.title.y = element_text(vjust = 1.5, size = 15),
#                       axis.ticks = element_line(size = 1.2),
#                       panel.grid.major = element_blank(),
#                       panel.grid.minor = element_blank(),
#                       panel.border = element_blank(),
#                       panel.background = element_blank(),
#                       plot.margin = unit(c(1,1,2,2), "lines"),
#                       text = element_text(size=15))) 
# 
# plot_beta <- plot_beta + scale_color_viridis(option = "plasma", discrete = T)
# 
# ggdraw() + 
#   draw_plot(plot_beta + theme(legend.justification = "top"), 
#             x = 0, y = 0, width = .9, height = 1) +
#   draw_plot(gghist_beta, x = 0.72, y = .030, width = .2, height = .75, scale = 1) 

#Alpha Change
gghist <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$slope50)) + 
  #geom_histogram(col = "black", fill = "black", bins = 30, binwidth = 0.25) +
  geom_density(alpha = 0.4, fill = "#450A69FF") +
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", alpha, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2, font_family = "serif") +
  coord_flip()


plot_alpha <- (ggplot(bbs_total_seg, aes(x = Year, y = SR50, group = BCR_name, 
                                     colour = factor(BCR_name), fill = factor(BCR_name))) + 
                 #labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = loess, se = T, aes(x = Year, y = SR50, group = BCR_name, 
                                                         colour = factor(BCR_name), fill = factor(BCR_name))) + 
                 #colour = factor(BCR, 
                 #                 labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 xlab("Year") +
                 ylab(expression(paste(alpha, "-diversity"))) +
                 labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black", size =1.2),
                       axis.text.x = element_text(size = 14, family = "serif"),
                       axis.text.y = element_text(size = 14, family = "serif"),
                       axis.title.x = element_text(vjust = -1, size = 14, family = "serif"),
                       axis.title.y = element_text(vjust = 1.5, size = 14, family = "serif"),
                       axis.ticks = element_line(size = 1.2),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.margin = unit(c(1,1,2,2), "lines"),
                       text = element_text(size=14, family = "serif"))) +
  scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)


#plot_alpha
plot_alpha <- ggdraw() + 
  draw_plot(plot_alpha + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist, x = 0.65, y = .030, width = .3, height = .75, scale = 1) 

ggsave(here::here("Figures/Figures_Diversity_Manuscript/AlphaBCR.tiff"), plot = plot_alpha,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


#Alpha Change
gghist <- ggplot(data = bbs_slopes_seg, aes(bbs_slopes_seg$slope.wet)) + 
  #geom_histogram(col = "black", fill = "black", bins = 30, binwidth = 0.25) +
  geom_density(alpha = 0.4, fill = "#9512A1FF") +
  labs(title = "") +
  labs(x = expression(paste("Slopes of ", alpha)), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2, font_family = "serif") +
  coord_flip()


plot_alpha <- (ggplot(bbs_total_seg, aes(x = Year, y = SR.wet, group = BCR_name, 
                                         colour = factor(BCR_name), fill = factor(BCR_name))) + 
                 #labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = loess, se = T, aes(x = Year, y = SR.wet, group = BCR_name, 
                                                         colour = factor(BCR_name), fill = factor(BCR_name))) + 
                 #colour = factor(BCR, 
                 #                 labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Gulf Coastal Prairie")))) +
                 xlab("Year") +
                 ylab(expression(paste("Wetland Bird ", alpha, "-diversity"))) +
                 labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black", size =1.2),
                       axis.text.x = element_text(size = 14, family = "serif"),
                       axis.text.y = element_text(size = 14, family = "serif"),
                       axis.title.x = element_text(vjust = -1, size = 14, family = "serif"),
                       axis.title.y = element_text(vjust = 1.5, size = 14, family = "serif"),
                       axis.ticks = element_line(size = 1.2),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.margin = unit(c(1,1,2,2), "lines"),
                       text = element_text(size=14, family = "serif"))) +
  scale_color_viridis(discrete = T, option = "plasma") + scale_fill_viridis(discrete = T, option = "plasma")


#plot_alpha
plot_alpha <- ggdraw() + 
  draw_plot(plot_alpha + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist, x = 0.65, y = .030, width = .3, height = .75, scale = 1) 

ggsave(here::here("Figures/Figures_Diversity_Manuscript/AlphaBCRwet.tiff"), plot = plot_alpha,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Alpha Change Wet

#Alpha Plot Corrected Vs Uncorrected 
plot_alpha_cor <- ggplot() + 
  geom_smooth(data = bbs_total_seg, aes(x = Year, y = Site_div.x), colour = "red") + 
  geom_smooth(data = bbs_total_seg, aes(x = Year, y = SR50), colour = "black")

plot_alpha_cor

#Beta Plot Facet Time series

plot_beta_cor <- ggplot(data = bbs_total_seg, aes(x = Year, y = beta50.jac, group = rteno.x)) + 
  geom_jitter(data = bbs_total_seg, aes(x = Year, y = beta50.jac), color = "lightgray") + 
  geom_smooth(method = lm, se = T, aes(colour = factor(BCR_name), fill = factor(BCR_name))) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta.mcmc), method = loess) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta50.jac), method = loess) 
  scale_color_viridis(discrete = T, option = "E") +
  scale_fill_viridis(discrete = T, option = "E") +
  xlab("Year") + 
  ylab(expression(paste("Detection-corrected ", beta[italic("Jac")]))) + 
  labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
  theme_bw() + 
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12, family = "serif"), legend.title = element_text(size = 12, face = "bold", family = "serif"),
        axis.title.x = element_text(size = 12, family = "serif"), axis.title.y = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif")) 
  
plot_beta_fac <- plot_beta_cor + facet_wrap(~ BCR_name, ncol = 2)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotBeta.tiff"), plot = plot_beta_fac,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Beta Plot Facet Time series Wetland Birds

plot_beta_cor_wet <- ggplot(data = bbs_total_seg, aes(x = Year, y = beta.wet.jac, group = rteno.x)) + 
  geom_jitter(data = bbs_total_seg, aes(x = Year, y = beta.wet.jac), color = "lightgray") + 
  geom_smooth(method = lm, se = T, aes(colour = factor(BCR_name), fill = factor(BCR_name))) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta.mcmc), method = loess) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta50.jac), method = loess) 
  scale_color_viridis(discrete = T, option = "magma") +
  scale_fill_viridis(discrete = T, option = "magma") +
  xlab("Year") + 
  ylab(expression(paste("Wetland Bird Detection-corrected ", " ", beta[italic("Jac")]))) + 
  labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
  theme_bw() + 
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12, family = "serif"), legend.title = element_text(size = 12, face = "bold", family = "serif"),
        axis.title.x = element_text(size = 12, family = "serif"), axis.title.y = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif")) 

plot_beta_fac_wet <- plot_beta_cor_wet + facet_wrap(~ BCR_name, ncol = 2)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotBetaWet.tiff"), plot = plot_beta_fac_wet,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Alpha Plot Facet Time series

plot_alpha_cor <- ggplot(data = bbs_total_seg, aes(x = Year, y = SR50, group = rteno.x)) + 
  geom_jitter(data = bbs_total_seg, aes(x = Year, y = SR50), color = "lightgray") + 
  geom_smooth(method = lm, se = T, aes(colour = factor(BCR_name), fill = factor(BCR_name))) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta.mcmc), method = loess) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta50.jac), method = loess) 
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  xlab("Year") + 
  ylab(expression(paste("Detection-corrected ", alpha, "-diversity"))) + 
  labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
  theme_bw() + 
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12, family = "serif"), legend.title = element_text(size = 12, face = "bold", family = "serif"),
        axis.title.x = element_text(size = 12, family = "serif"), axis.title.y = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif")) 

plot_alpha_fac <- plot_alpha_cor + facet_wrap(~ BCR_name, ncol = 2)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotAlpha.tiff"), plot = plot_alpha_fac,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

#Alpha plot wetland birds
plot_alpha_cor_wet <- ggplot(data = bbs_total_seg, aes(x = Year, y = SR.wet, group = rteno.x)) + 
  geom_jitter(data = bbs_total_seg, aes(x = Year, y = SR.wet), color = "lightgray") + 
  geom_smooth(method = lm, se = T, aes(colour = factor(BCR_name), fill = factor(BCR_name))) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta.mcmc), method = loess) + 
  #geom_smooth(data = bbs_total_seg, aes(x = Year, y = beta50.jac), method = loess) 
  scale_color_viridis(discrete = T, option = "plasma") +
  scale_fill_viridis(discrete = T, option = "plasma") +
  xlab("Year") + 
  ylab(expression(paste("Wetland Bird Detection-corrected ", " ", alpha, "-diversity"))) + 
  labs(colour = "Bird Conservation Region", fill = "Bird Conservation Region") +
  theme_bw() + 
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        legend.text = element_text(size = 12, family = "serif"), legend.title = element_text(size = 12, face = "bold", family = "serif"),
        axis.title.x = element_text(size = 12, family = "serif"), axis.title.y = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif")) 

plot_alpha_fac_wet <- plot_alpha_cor_wet + facet_wrap(~ BCR_name, ncol = 2)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotAlphaWet.tiff"), plot = plot_alpha_fac_wet,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


#Plotting Climate Data for Entire Time Series 
bcr <- seg.df[, c("Site", "BCR")]
climate <- merge(climate, bcr, by = "Site")

climate1 <- climate[, c("Year", "tmin.c", "tmax.c", "tmean.c", "tmax.bird.c", "tmean.bird.c", "tmin.bird.c")]
climate1 <- reshape2::melt(climate1, id = c("Year"))
climate1 <- climate1 %>% group_by(Year, variable) %>% summarize(mean = mean(value))

climate1$variable_f <- factor(climate1$variable, levels = c("tmax.c", "tmax.bird.c", "tmean.c", "tmean.bird.c", "tmin.c", "tmin.bird.c"))

levels(climate1$variable_f) <- c(tmax.c = expression("T"["max"]), 
                                 tmax.bird.c = expression("Breeding Season T"["max"]), 
                                 tmean.c = expression("T"["mean"]), 
                                 tmean.bird.c = expression("Breeding Season T"["mean"]), 
                                 tmin.c = expression("T"["min"]),
                                 tmin.bird.c = expression("Breeding Season T"["min"]))


climate.plot.annual <- ggplot() +
  geom_point(data = climate1, aes(x = Year, y = mean, group = variable)) +
  #geom_line(data = climate1, aes(x = Year, y = tmax.c), color = "Red") +
  #geom_line(data = climate1, aes(x = Year, y = tmean.c), color = "Black") +
  geom_smooth(data = climate1, method = lm, aes(x = Year, y = mean, group = variable_f, colour = factor(variable_f), fill = factor(variable_f))) +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Temperature ", "(", degree, "C)"))) +
  xlim(1975, 2018) +
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"),
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0, "cm")) +
  scale_color_manual(values = c("#444F6BFF", "#444F6BFF", "#838079FF", "#838079FF", "#C7B76BFF", "#C7B76BFF")) +
  scale_fill_manual(values = c("#444F6BFF", "#444F6BFF", "#838079FF", "#838079FF", "#C7B76BFF", "#C7B76BFF")) +
  facet_wrap(~variable_f, nrow = 3, scales = "free", labeller = label_parsed)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotClimate.tiff"), plot = climate.plot.annual,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()


#Precip Data
climate2 <- climate[, c("Year", "precip_Wet.c", "precip_Dry.c")]
climate2 <- reshape2::melt(climate2, id = c("Year"))
climate2 <- climate2 %>% group_by(Year, variable) %>% summarize(mean = mean(value))

climate2$variable_f <- factor(climate2$variable, levels = c("precip_Wet.c", "precip_Dry.c"))

levels(climate2$variable_f) <- c(#pmean.c = expression("Precip"["mean"]), 
                                 precip_Wet.c = expression("Wet Season Precip"["mean"]), 
                                 precip_Dry.c = expression("Dry Season Precip"["mean"]))


precip.plot.annual <- ggplot() +
  geom_point(data = climate2, aes(x = Year, y = mean, group = variable_f)) +
  #geom_line(data = climate1, aes(x = Year, y = tmax.c), color = "Red") +
  #geom_line(data = climate1, aes(x = Year, y = tmean.c), color = "Black") +
  geom_smooth(data = climate2, method = lm, aes(x = Year, y = mean, group = variable_f, colour = factor(variable_f), fill = factor(variable_f))) +
  theme_bw() +
  xlab("Year") +
  ylab("Precipitation (cm)") +
  xlim(1975, 2018) +
  theme(strip.text = element_text(size = 12, face = "bold", family = "serif"),
        strip.background = element_blank(),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"),
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0, "cm")) +
  scale_color_manual(values = c("#00204DFF", "#FFEA46FF", "#C7B76BFF")) +
  scale_fill_manual(values = c("#00204DFF", "#FFEA46FF", "#C7B76BFF")) +
  facet_wrap(~variable_f, nrow = 1, scales = "free", labeller = label_parsed)

big.clim <- plot_grid(climate.plot.annual, precip.plot.annual, nrow = 2, rel_heights = c(2.5, 1), align = "v")
 

ggsave(here::here("Figures/Figures_Diversity_Manuscript/FacetPlotBigClimate.tiff"), plot = big.clim,
       device = "tiff", width = 8, height = 8, units = "in", dpi = 600)

dev.off()

#LULC Data 
lulc.df <- seg.df[, c("Site", "diff.from.first.man", "diff.from.first.for", "diff.from.first.ww", "diff.from.first.ew", "diff.from.first.ur", "diff.from.first.wat", "diff.from.first.bare", "diff.from.first.ag", "diff.from.first.wet", "diff.from.first.grass", "diff.from.first.human", "diff.from.first.nat")]
lulc.df <- reshape2::melt(lulc.df, id = "Site")

lulc.df$variable_f <- factor(lulc.df$variable, levels = c("diff.from.first.man", "diff.from.first.ww", "diff.from.first.ew", "diff.from.first.wet", "diff.from.first.ur", "diff.from.first.ag", "diff.from.first.grass", "diff.from.first.for", "diff.from.first.wat", "diff.from.first.bare", "diff.from.first.human", "diff.from.first.nat"))

levels(lulc.df$variable_f) <- c(diff.from.first.man = expression(paste(Delta, "% Cover Mangrove")),
                                diff.from.first.wwnm = expression(paste(Delta, "% Cover Woody Wetland")),
                                diff.from.first.ew = expression(paste(Delta, "% Cover Emergent Wetland")),
                                diff.from.first.wet = expression(paste(Delta, "% Cover Total Wetland")),
                                diff.from.first.ur = expression(paste(Delta, "% Cover Urban")),
                                diff.from.first.ag = expression(paste(Delta, "% Cover Agriculture")),
                                diff.from.first.grass = expression(paste(Delta, "% Cover Grassland")),
                                diff.from.first.for = expression(paste(Delta, "% Cover Forest")),
                                diff.from.first.wat = expression(paste(Delta, "% Cover Water")),
                                diff.from.first.bare = expression(paste(Delta, "% Cover Bare")),
                                diff.from.first.human = expression(paste(Delta, "% Cover Anthropogenic")),
                                diff.from.first.nat = expression(paste(Delta, "% Cover Upland")))
                                
  


lulc.plot <- ggplot(data = lulc.df, aes(x = value, fill = variable)) +
  geom_histogram(bins = 50) + facet_wrap(~variable_f, labeller =label_parsed, nrow = 4) + scale_fill_viridis(discrete = T, option = "E") +
  xlab("Difference in Percent Cover") +
  ylab("Frequency") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text = element_text(size = 12, family = "serif"))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/lulchist.tiff"), plot = lulc.plot,
       device = "tiff", width = 8, height = 5, units = "in", dpi = 600)

dev.off()

