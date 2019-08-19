###############################################################################################
#########################Model runs for the last year of each##################################
###############################################################################################

#max.anom.s & p.anom.f explain most variation in multiple LM for climate (20% together)
#This model explains 32% variation w/o alpha
mod.last <- lm(data = bbs_last, beta.mcmc ~ alpha.mcmc + p.anom + Latitude + p.anom)
summary(mod.last)
Anova(mod.last)

ggplot(data = bbs_last, aes(x = Latitude, y = scale.cmrl)) + geom_smooth(method = loess)

summary(lm(bbs_last$beta.mcmc ~ bbs_last$mean.anom.bird))

#Again anomalies of fall precipitation interestingly explaining variation, tmax.c, tmean.c, tmax.bird.c,  
#Sig LULC terms: Pct_wetland, scale.pdew, scale.pdur (negative relationship), scale.pdwet,
#
mod.last.a <- lm(data = bbs.last.corra, alpha50.change ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ur)
summary(mod.last.a)
Anova(mod.last.a)
avPlots(mod.last.a)

cookd.a <- cooks.distance(mod.last.a)
plot(cookd.a)
abline(h = 4*mean(cookd.a))
influential.a <- as.numeric(names(cookd.a)[(cookd.a > 4*mean(cookd.a, na.rm = T))])
influential.a
bbs.last.corra <- bbs_last[-influential.a, ]


bbs_last <- bbs_last %>% mutate(alpha50.s = scale(alpha50))

betareg.last50 <- betareg::betareg(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + p.anom.dry + diff.from.first.ww, link = "logit")
summary(betareg.last50)

summary(mod.last50)

reg1 <- lm(data = bbs.last.corr, beta50 ~ p.anom.wet + diff.from.first.ww)
res1 <- resid(reg1)
reg2 <- lm(data = bbs.last.corr, alpha50 ~ p.anom.wet + diff.from.first.ww)
res2 <- resid(reg2)

reg3 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + diff.from.first.ww)
res3 <- resid(reg3)
reg4 <- lm(data = bbs.last.corr, p.anom.wet ~ alpha50 + diff.from.first.ww)
res4 <- resid(reg4)

reg5 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet)
res5 <- resid(reg5)
reg6 <- lm(data = bbs.last.corr, diff.from.first.ww ~ alpha50 + p.anom.wet)
res6 <- resid(reg6)

#reg7 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + diff.from.first.ww)
#res7 <- resid(reg7)
#reg8 <- lm(data = bbs.last.corr, change.cmrl ~ alpha50 + p.anom.wet + diff.from.first.ww)
#res8 <- resid(reg8)

mod.a <- lm(res1 ~ res2)
mod.p <- lm(res3 ~ res4)
mod.ww <- lm(res5 ~ res6)
#mod.cmrl <- lm(res7 ~ res8)

ggplotRegression(mod.a)
ggplotRegression(mod.p)
ggplotRegression(mod.ww)
ggplotRegression(mod.cmrl)

reg.a <- lm(data = bbs_last, beta50 ~ alpha50)
ggplotRegression(reg.a)

reg.p <- lm(data = bbs_last, beta50 ~ p.anom.wet)
ggplotRegression(reg.p)

reg.ww <- lm(data = bbs_last, beta50 ~ diff.from.first.ww)
ggplotRegression(reg.ww)

mod.last50 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + diff.from.first.ww)
summary(mod.last50)
avPlots(mod.last50)

cookd <- cooks.distance(mod.last50)
plot(cookd)
abline(h = 4*(mean(cookd, na.rm = T)))

influential <- as.numeric(names(cookd)[(cookd > 4*mean(cookd, na.rm = T))])
bbs.last.corr <- bbs_last[-influential, ]

library(ggplot2)
ggplot(bbs.last.corr, aes(x = diff.from.first.ww, y = beta50)) +
  geom_point(size = 4, aes(fill = NULL), shape = 21) +
  scale_fill_grey() +
  # geom_line(aes(y = predict(gy_loglog, GasolineYield),
  #               colour = "log-log", linetype = "log-log")) +
  geom_line(aes(y = predict(betareg.last50, bbs.last.corr), 
                colour = "logit", linetype = "logit")) +
  scale_colour_manual("", values = c("red", "blue")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()

pacman::p_load(ggmap, maps)

states <- map_data("state")
dim(states)
head(states)
gom <- subset(states, region %in% c("texas", "florida", "alabama", "arkansas", "mississippi", "louisiana", "georgia"))
map_gom <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(bbs_last$Longitude) - 1.5, max(bbs_last$Longitude)), ylim = c(min(bbs_last$Latitude), max(bbs_last$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = bbs_last, aes(x = Longitude, y = Latitude, color = beta50), size = 3)
box.thing <- make_bbox(lon = bbs_last$Longitude, lat = bbs_last$Latitude, f = 0.1) 

map_gom + scale_color_viridis(name = expression(paste(beta, "-diversity"))) + theme_map() + theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  
################################################################################################################################
#######################################Confirming Assumptions of Data before selecting Model####################################
################################################################################################################################
p_load(HLMdiag, DHARMa, car)

Plot.beta.linearity <- plot(resid(mod.temp.beta), bbs_full$beta.reg)

p_load(MASS)

bbs_full$mean.anom.t <- bbs_full$mean.anom + 1
qqp(bbs_full$mean.anom.t, "norm")
qqp(bbs_full$mean.anom.t, "lnorm")

################################################################################################################################
############################################***Models for Raw Species Observations***###########################################
################################################################################################################################


#mod.full.beta <- lmer(data = bbs_full, beta.reg ~ anomalies + precip.anomalies + ratio.ww + pct.ur + (1|site.x))
mod.full.beta <- lmer(data = bbs_full, beta.reg ~ alpha + mean.anom.bird + change.cmrl.km + p.anom.bird + scale.pdww + scale.pdur + (1|Site))
summary(mod.full.beta)
Anova(mod.full.beta)

sjPlot::tab_model(mod.full.beta)

#Model with full 
mod.temp.beta <- lmer(data = bbs_full, beta.reg ~ mean.anom + p.anom + (1|Site), REML = F)
summary(mod.temp.beta)
Anova(mod.temp.beta)


#Models with first and last betas 
mod.short.beta <- lmer(data = bbs_short, beta.reg ~ mean.anom.bird + p.anom.bird + scale.pdww + 
                         scale.pdew + scale.pdur + scale.ag + (1|site.x))
summary(mod.short.beta)
Anova(mod.short.beta)

temp.est <- Effect("mean.anom.bird", partial.residuals = T, mod.full.beta)
plot(temp.est)

ww.est <- Effect("scale.pdww", partial.residuals = T, mod.full.beta)
plot(ww.est)

#Models with the last beta
mod.last.beta <- lmer(data = bbs_short, beta.reg ~ anomalies + precip.anomalies + scale.ww + 
                        scale.ew + scale.ur + scale.ag + (1|site.x))

summary(mod.last.beta)
Anova(mod.last.beta)

temp.est <- Effect("anomalies", partial.residuals = T, mod.last.beta)
plot(temp.est)

ww.est <- Effect("scale.ww", partial.residuals = T, mod.last.beta)
plot(ww.est)

beta.ur <- lm(data = bbs_last, beta.mcmc ~ pct.ur)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))}#, x = "Percent Urban", y = ("Beta Diversity"))}

ggplotRegression(beta.ur)

avPlots(mod.last)

mod.last.a <- lm(data = bbs_last, alphamcmc.change ~ p.anom.wet + max.anom.wet +  p.anom.dry)
summary(mod.last.a)




##############################################################################################################################
###################################***Models for MCMC estimated Species Occurrences***########################################
##############################################################################################################################

#mod.full.betamcmc <- lmer(data = bbs_full, beta.mcmc.y ~ anomalies + precip.anomalies + pct.ww + pct.ew + pct.ur + (1|site.x)) 
mod.full.betamcmc <- lmer(data = bbs_full, beta.mcmc ~ alpha.mcmc + mean.anom.bird + change.cmrl.km + p.anom.bird + scale.pdww + scale.pdur + (1|Site))
summary(mod.full.betamcmc)
Anova(mod.full.betamcmc)

mod.h.betamcmc <- lmer(data = bbs_full, beta.mcmc ~ mean.anom.bird + scale.pdww + scale.pdur + (1|Site), REML = F)
summary(mod.h.betamcmc)
Anova(mod.h.betamcmc)

mod.betamcmc <- lmer(data = bbs_full, beta.mcmc ~ min.anom.s + (1|Site), REML = F)
summary(mod.betamcmc)
Anova(mod.betamcmc)

test <- lmer.df
test <- separate(test, site.y, into = c("rteno", "segment"), by = "_")

p_load(effects)

temp.est <- Effect("mean.anom.bird", partial.residuals = T, mod.full.betamcmc)
plot(temp.est)

ww.est <- Effect("scale.pdww", partial.residuals = T, mod.full.betamcmc)
plot(ww.est)

ur.est <- Effect("scale.pdur", partial.residuals = T, mod.full.betamcmc)
plot(ur.est)

sjPlot::tab_model(mod.full.beta)
sjPlot::tab_model(mod.full.betamcmc)

mod.temp <- lmer(data = bbs_full, tmean.c ~ Year + (1|Site), REML = F)
summary(mod.temp)
Anova(mod.temp)


mod.alpha <- lm(data = bbs_full, scale.alphamcmc ~ scale.pdww + scale.pdur)# + (1|Site), REML = F)
summary(mod.alpha)
Anova(mod.alpha)


ur.est.alpha <- Effect("scale.pdur", partial.residuals = T, mod.alpha)
plot(ur.est.alpha)

bbs.bin7 <- bbs_full[bbs_full$Yr_bin.x == 7, ]
bbs.bin7 <- bbs.bin7 %>% mutate(alpha.scale = scale(alpha.mcmc))
bbs.bin7 <- bbs.bin7 %>% mutate(ratio.ww = round(ratio.ww, digits = 2))

##################################################################################
#############################Spatial NMDS and Adonis##############################
##################################################################################

pacman::p_load(ggfortify)

#See which year had the most completed surveys
setDT(bbs_div_means)[, .(count = uniqueN(site)), by = Year] #2008 with 174

#bbs_lulc <- bbs_lulc[complete.cases(bbs_lulc), ]

bbs.mds <- bbs_simple[bbs_simple$Year == 1980,]
bbs.mds <- bbs.mds %>% arrange(site)

lc.mds <- bbs_lulc %>% dplyr::select(site, Year, ratio.ww, pct.ur)
lc.mds <- lc.mds[lc.mds$Year == 1980,]

mds.list <- unique(bbs.mds$site)
wetland.site <- unique(lc.mds$site)

mds.full <- merge(lc.mds, bbs.mds, by = "site")

lc.mds <- dplyr::select(mds.full, c(ratio.ww, pct.ur, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$groups <- NA
ww.mds$groups[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
#ww.mds$groups[ww.mds$ratio.ww >= 0.8 & ww.mds$ratio.ww <= 1.2] <- "Mix"
ww.mds$groups[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds$groups

ur.mds <- lc.mds %>% arrange(site)
ur.mds$groups <- NA
ur.mds$groups[ur.mds$pct.ur <= .33] <- "L"
ur.mds$groups[ur.mds$pct.ur > .33 & ur.mds$pct.ur < .66] <- "M"
ur.mds$groups[ur.mds$pct.ur >= .66] <- "H"
ur.mds <- ur.mds$groups

bbs.mds <- bbs.mds %>% arrange(site)
mds_cast <- dcast(mds.full, unique_id ~ sci_name, value.var = "Detected", fun.aggregate = sum)
rownames(mds_cast) <- mds_cast$unique_id
mds_cast <- mds_cast[, -1]


mds <- metaMDS(mds_cast, trymax = 50)


data.scores <- as.data.frame(scores(mds))  
data.scores$site <- rownames(data.scores)  
data.scores$grp1 <- ww.mds
#data.scores$grp2 <- ur.mds


mds_plot <- ggplot(data = data.scores) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, colour = ww.mds), level = 0.50) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = , colour = ww.mds), size=2) +
  scale_fill_manual(values = )

mds_plot + labs(color = "Wetland Cover Types")

adon.results <- adonis(mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)

#autoplot(prcomp(), data = data.scores, )

########################################################################################
#########################MDS for the Occupancy Cutoffs##################################

occ.mds <- occ50[occ50$Year == 1980,]
occ.mds <- occ.mds %>% arrange(site)

lc.mds <- bbs_lulc %>% dplyr::select(site, Year, ratio.ww, pct_wetland)
lc.mds <- lc.mds[lc.mds$Year == 1980,]

occ.mds.list <- unique(occ.mds$site)
wetland.site <- unique(lc.mds$site)

occ.mds.full <- merge(lc.mds, occ.mds, by = "site")

lc.mds <- dplyr::select(occ.mds.full, c(ratio.ww, pct_wetland, site))
lc.mds <- lc.mds[!duplicated(lc.mds), ]

#ww.mds <- ww.mds[ww.mds$site %in% mds.list,]

#ww11.site <- unique(ww.mds$site.x)
ww.mds <- lc.mds %>% arrange(site)
ww.mds$groups <- NA
ww.mds$groups[ww.mds$ratio.ww <= 1] <- "Emergent Wetland Dominated"
ww.mds$groups[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mixed"
ww.mds$groups[ww.mds$ratio.ww > 1] <- "Woody Wetland Dominated"
ww.mds <- ww.mds$groups

# wet.mds <- lc.mds %>% arrange(site)
# wet.mds$groups <- NA
# wet.mds$groups[wet.mds$pct_wetland < .5] <- "L"
#wet.mds$groups[wet.mds$pct_wetland > .33 & wet.mds$pct_wetland < .66] <- "M"
# wet.mds$groups[wet.mds$pct_wetland >= .5] <- "H"
# wet.mds <- wet.mds$groups

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
                      labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                      values = c("#0072B2", "#009E73", "#999999")) +
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                     values = c(0, 16, 3))

mds_plot #+ labs(color = "Wetland Cover Types", shape = "Wetland Cover Type") + scale_shape_manual(values = c(0, 16, 3))


adon.results <- adonis(occ_mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)








############################################################################################
##############################Plotting######################################################
############################################################################################

#CMRL 
line = "#0000CC"
fill = "#6699CC"
box1 <- ggplot(bbs_full, aes(x = Yr_bin, y = CMRL, group = site.y)) +
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





mod.1 <- lm(beta.mcmc ~ Latitude + p.anom, data = bbs_last)
res1 <- resid(mod.1)
mod.2 <- lm(alpha.mcmc ~ Latitude + p.anom, data = bbs_last)
res2 <- resid(mod.2)

mod.alpha <- lm(res1 ~ res2)

mod.3 <- lm(beta.mcmc ~ alpha.mcmc + p.anom, data = bbs_last)
res3 <- resid(mod.3)
mod.4 <- lm(Latitude ~ alpha.mcmc + p.anom, data = bbs_last)
res4 <- resid(mod.4)

mod.lat <- lm(res3 ~ res4)

mod.5 <- lm(beta.mcmc ~ alpha.mcmc + Latitude, data = bbs_last)
res5 <- resid(mod.5)
mod.6 <- lm(p.anom ~ alpha.mcmc + Latitude, data = bbs_last)
res6 <- resid(mod.6)

mod.panom <- lm(res5 ~ res6)

mod.7 <- lm(beta.mcmc ~ alpha.mcmc + p.anom + Latitude, data = bbs_last)
res7 <- resid(mod.7)
mod.8 <- lm(diff.from.first.ww ~ alpha.mcmc + p.anom + Latitude, data = bbs_last)
res8 <- resid(mod.8)

mod.ww <- lm(res7 ~ res8)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("R-squared = ",signif(summary(fit)$r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)), 
         #x = expression(paste("Detection-corrected ", alpha, "-Diversity")),
         #x = "Adjusted Precipitation Anomalies (cm)",
         x = "Adjusted Latitude (Decimal Degrees)",
         y = expression(paste("Adjusted ", 
                              beta, "-Diversity")))}


ggplotRegression(mod.alpha)
ggplotRegression(mod.panom)
ggplotRegression(mod.lat)
ggplotRegression(mod.ww)





pacman::p_load(partial_plots)




avPlots(mod.last)



cbbPalette <- c("#99FFCC","#99FFFF", "#66CCCC", "#339999", "#006666")

beta.plot <- (ggplot(bbs_total, aes(x = Year, y = beta, group = BCR, 
                                    colour = factor(BCR, 
                                                    labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Tamaulipan Brushlands", "Gulf Coastal Prairie")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(aes(x = Year, y = beta, group = BCR, 
                                colour = factor(BCR, 
                                                labels = c("Mississippi Alluvial Valley", "Southeastern Coastal Plain", "Peninsular Florida", "Tamaulipan Brushlands", "Gulf Coastal Prairie"))), linetype = 1, size = 1, method = "lm") +
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






#######################################################################
######Link Species Name with Range and looking at tropicalization######
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










########################################################################################

#Below is mostly junk

########################################################################################
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
alpha.1980 <- alpha.1980[, c("site", "yr_bin", "Site_div")]
alpha.agg <- aggregate(. ~ site + yr_bin, alpha.1980, FUN = "mean")
alpha.agg$unique_id <- paste0(alpha.agg$site, "_", alpha.agg$yr_bin)

test <- merge(test, alpha.agg, "unique_id")
test <- test %>% group_by(site) %>% arrange(Yr_bin) %>% 
  mutate(p.change.alpha = (Site_div - first(Site_div)) / first(Site_div) * 100) %>%
  ungroup(test)
colnames(test)[colnames(test) == "site.x"] <- "Site"
test <- test[, c("unique_id", "Site", "Yr_bin", "Site_div", 
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

mod.beta <- lmer(beta ~ anomalies + Site_div + precip.anomalies + (1|Rteno), data = test)
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
