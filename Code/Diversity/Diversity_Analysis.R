##############################################################################################
##################Script created for doing analyses on GoM bird diversity manuscript##########
####################################Created By: Spencer R Keyser##############################
##############################################8/19/2019#######################################
##############################################################################################

rm(list = ls())

#Load in the necessary packages
pacman::p_load("here", "tidyverse", "reshape2", "vegan", "data.table", "cowplot", "lme4", "sjPlot", 
               "sjstats", "car", "ggfortify", "ggmap", "maps", "viridis", 'MuMIn')


#load in the workspace with data 
#load(here::here("R Workspace/Diversity_Script_Clean.RData"))
load(here::here("R Workspace/Diversity_Script_Full.RData"))


#Mapping distribution of Beta-diversity across GoM

states <- map_data("state")
dim(states)
head(states)
gom <- subset(states, region %in% c("texas", "florida", "alabama", "arkansas", "mississippi", "louisiana", "georgia"))
map_gom <- ggplot(data = gom) + geom_polygon(aes(x = long, y = lat, group = group), fill = "gray", color = "black") + coord_fixed(xlim = c(min(bbs_last$Longitude) - 1.5, max(bbs_last$Longitude)), ylim = c(min(bbs_last$Latitude), max(bbs_last$Latitude) + 2.5), ratio = 1.2) +
  geom_point(data = bbs_last, aes(x = Longitude, y = Latitude, color = beta50), size = 3)
 
map_gom + scale_color_viridis(name = expression(paste(beta, "-diversity"))) + theme_map() + theme(legend.position = c(0.09, .75), legend.title = element_text(size = 15), legend.text = element_text(size = 15))  


##################################################################################
#############################Spatial NMDS and Adonis##############################
##################################################################################
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
ww.mds$groups[ww.mds$ratio.ww <= 0.5] <- "Emergent Wetland Dominated"
ww.mds$groups[ww.mds$ratio.ww >= 0.5 & ww.mds$ratio.ww <= 1.5] <- "Mix"
ww.mds$groups[ww.mds$ratio.ww > 1.5] <- "Woody Wetland Dominated"
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
  geom_point(aes(x = NMDS1, y = NMDS2, colour = ww.mds, shape = ww.mds), size=2) +
  scale_colour_manual(name = "Wetland Cover Types", 
                      labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                      values = c("#0072B2", "#009E73", "#999999")) +
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                     values = c(0, 16, 3))

mds_plot 

adon.results <- adonis(mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)


########################################################################################
#########################MDS for the Occupancy Cutoffs##################################
########################################################################################

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
ww.mds$groups[ww.mds$ratio.ww > 1.5] <- "Woody Wetland Dominated"
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
                      labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                      values = c("#0072B2", "#009E73", "#999999")) +
  scale_shape_manual(name = "Wetland Cover Types",
                     labels = c("Emergent Wetland Dominated", "Mix", "Woody Wetland Dominated"),
                     values = c(0, 16, 3))

mds_plot #+ labs(color = "Wetland Cover Types", shape = "Wetland Cover Type") + scale_shape_manual(values = c(0, 16, 3))


adon.results <- adonis(occ_mds_cast ~ ww.mds, method = "jaccard", perm = 999)


print(adon.results)

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



###############################################################################################
#########################Model runs for the last year of each##################################
###############################################################################################

#Occ 50 Results

#This model explains 35% variation
mod.last <- lm(data = bbs_last, beta50 ~ alpha50 + p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
summary(mod.last)
Anova(mod.last)

#Pull out outliers from the data set 
mod.last.outlier <- cooks.distance(mod.last)
plot(mod.last.outlier)
abline(h = 4*mean(mod.last.outlier))
influential.points <- as.numeric(names(mod.last.outlier)[(mod.last.outlier > 4*mean(mod.last.outlier, na.rm = T))])

#Remove the outliers from the data 
bbs.last.corr <- bbs_last[-influential.points, ]

#Run the model without the outliers 
mod.last <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
summary(mod.last)
#Anova(mod.last)

#Partial Regression Models 
reg1 <- lm(data = bbs.last.corr, beta50 ~ p.anom.wet + diff.from.first.ww)
res1 <- resid(reg1)
reg2 <- lm(data = bbs.last.corr, alpha50 ~ p.anom.wet + diff.from.first.ww)
res2 <- resid(reg2)

reg3 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + diff.from.first.ww)
res3 <- resid(reg3)
reg4 <- lm(data = bbs.last.corr, p.anom.wet ~ alpha50 + diff.from.first.ww)
res4 <- resid(reg4)

reg5 <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + p.anom.dry + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
res5 <- resid(reg5)
reg6 <- lm(data = bbs.last.corr, mean.anom.bird ~ alpha50 + p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.ew + diff.from.first.ur + Duration)
res6 <- resid(reg6)

mod.a <- lm(res1 ~ res2)
mod.p <- lm(res3 ~ res4)
mod.ww <- lm(res5 ~ res6)

ggplotRegression(mod.a)
ggplotRegression(mod.p)
ggplotRegression(mod.ww)

#Again anomalies of fall precipitation interestingly explaining variation, tmax.c, tmean.c, tmax.bird.c,  
#Sig LULC terms: Pct_wetland, scale.pdew, scale.pdur (negative relationship), scale.pdwet,
#
mod.last.a <- lm(data = bbs_last, alpha50.pchange ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)  #+ diff.from.first.ur + diff.from.first.ww)
summary(mod.last.a)
#Anova(mod.last.a)
avPlots(mod.last.a)
#tab_model(mod.last.a)

cookd.a <- cooks.distance(mod.last.a)
plot(cookd.a)
abline(h = 4*mean(cookd.a))
influential.a <- as.numeric(names(cookd.a)[(cookd.a > 4*mean(cookd.a, na.rm = T))])
influential.a
bbs.last.corra <- bbs_last[-influential.a, ]

mod.last.a <- lm(data = bbs.last.corra, alpha50.pchange ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration) # + diff.from.first.ur + diff.from.first.ww)
summary(mod.last.a)
Anova(mod.last.a)
#tab_model(mod.last.a)

reg.a1 <- lm(data = bbs.last.corra, alpha50.change ~ p.anom.wet + p.anom.dry) #+ diff.from.first.ur + diff.from.first.ww)
res.a1 <- resid(reg.a1)
reg.a2 <- lm(data = bbs.last.corra, max.anom.sp ~ p.anom.wet + p.anom.dry)# + diff.from.first.ur + diff.from.first.ww)
res.a2 <- resid(reg.a2)

reg.a3 <- lm(data = bbs.last.corra, alpha50.change ~ max.anom.sp + p.anom.dry)# + diff.from.first.ur + diff.from.first.ww)
res.a3 <- resid(reg.a3)
reg.a4 <- lm(data = bbs.last.corra, p.anom.wet ~ max.anom.sp + p.anom.dry)# + diff.from.first.ur + diff.from.first.ww)
res.a4 <- resid(reg.a4)

reg.a5 <- lm(data = bbs.last.corra, alpha50.change ~ max.anom.sp + p.anom.wet)# + diff.from.first.ur + diff.from.first.ww)
res.a5 <- resid(reg.a5)
reg.a6 <- lm(data = bbs.last.corra, p.anom.dry ~ max.anom.sp + p.anom.wet)# + diff.from.first.ur + diff.from.first.ww)
res.a6 <- resid(reg.a6)

# reg.a7 <- lm(data = bbs.last.corra, alpha50.change ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ww)
# res.a7 <- resid(reg.a7)
# reg.a8 <- lm(data = bbs.last.corra, diff.from.first.ur ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ww)
# res.a8 <- resid(reg.a8)
# 
# reg.a9 <- lm(data = bbs.last.corra, alpha50.change ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ur)
# res.a9 <- resid(reg.a9)
# reg.a10 <- lm(data = bbs.last.corra, diff.from.first.ww ~ max.anom.sp + p.anom.wet + p.anom.dry + diff.from.first.ur)
# res.a10 <- resid(reg.a10)

mod.t1 <- lm(res.a1 ~ res.a2)
mod.p1 <- lm(res.a3 ~ res.a4)
mod.p2 <- lm(res.a5 ~ res.a6)
# mod.ur <- lm(res.a7 ~ res.a8)
# mod.ww1 <- lm(res.a9 ~ res.a10)

ggplotRegression(mod.t1)
ggplotRegression(mod.p1)
ggplotRegression(mod.p2)
# ggplotRegression(mod.ww1)
# ggplotRegression(mod.ur)

####################Wetland Bird Analyses#################################

#Beta Div
mod.last <- lm(data = bbs_last, beta.wet ~ alpha.wet + p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
summary(mod.last)
avPlots(mod.last)

#Pull out outliers from the data set 
mod.last.outlier <- cooks.distance(mod.last)
plot(mod.last.outlier)
abline(h = 4*mean(mod.last.outlier))
influential.points <- as.numeric(names(mod.last.outlier)[(mod.last.outlier > 4*mean(mod.last.outlier, na.rm = T))])

#Remove the outliers from the data 
bbs.last.corr <- bbs_last[-influential.points, ]
weird.sites <- bbs_last[influential.points, ]


#Run the model without the outliers 
mod.last <- lm(data = bbs.last.corr, beta50 ~ alpha50 + p.anom.wet + p.anom.dry + mean.anom.bird + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
summary(mod.last)
#Anova(mod.last)


#Alpha Div
mod.last.a <- lm(data = bbs_last, alphawet.pchange ~ mean.anom.bird + p.anom.wet + p.anom.dry + diff.from.first.ww + diff.from.first.ew + diff.from.first.ur + Duration)
summary(mod.last.a)


























#####Running models with a Beta regression#####

bbs_last <- bbs_last %>% mutate(alpha50.s = scale(alpha50))

betareg.last50 <- betareg::betareg(data = bbs_last, beta50 ~ alpha50.s + p.anom.wet.s + p.anom.dry.s + scale.pdww, link = "logit")
summary(betareg.last50)

pacman::p_load(betareg)
install.packages("betareg")
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





############################################################################################
##############################Plotting######################################################
############################################################################################

#######################################################################
######Link Species Name with Range and looking at tropicalization######
#######################################################################

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




