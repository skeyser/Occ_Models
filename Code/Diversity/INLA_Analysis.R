##############################################################
#############INLA Analysis Mangrove###########################
###############Created By: Spencer R Keyser###################
####################8/18/2020#################################

#################
###Hello world###
#################

rm(list = ls())

#Install packages
if ("INLA" %in% rownames(installed.packages()) == FALSE){
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

#Load Packages 
pacman::p_load("tidyverse", "MASS", "INLA", "corrplot")
select <- dplyr::select

#Add in some predefined functions
#Logit and Expit function
logit <- function(x){
  log(x/(1-x))
}

expit <- function(x){
  exp(x)/(1 + exp(x))
}


#Bring in data 
df <- read.csv(here::here("Data_BBS/Generated DFs/INLAdfFinalRevised2.csv"))

#Look at correlations in the data 
#Pull covariates
df.pred <- df %>%
  select(pct.man, pct.ww, pct.ur, pct.ww, pct.ew, pct.wetland, pct.human, pct.ag, pct.wat, pct.bare, pct.for, pct.grass, pct.nat)

#Calculate correlations
cors <- cor(df.pred)
cor.df <- as.data.frame(round(cors, 2))

#Look at correlations
corrplot(cors, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, method = "number")

#Scale Covariates 
#DF created for running first model for species richness
model.df <- df %>%
  select(Route, Site, Year, Yr_bin, BCR, alpha, alpha.full, beta, beta.full, turn, nest, beta, tmean, tmin, tmax, ppt, 
         pct.man, pct.ur, pct.ww, pct.wwnm, pct.wetland, pct.ew, pct.ag, pct.for, pct.grass, pct.nat, ratio.ew) %>%
  mutate_at(c("tmean", "ppt", "pct.man", "pct.ur", "pct.wwnm", "ratio.ew", "pct.ww", "pct.ew", "pct.ag", "pct.for", "pct.grass", "pct.nat"), scale)

#Check out the distribution of covariates and response variables
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=F)
  return(lnew)
}


hist.df <- df %>%
  select(matches("pct")) %>%
  pivot_longer(cols = matches("pct"), names_to = "land_cov", values_to = "pct") %>%
  mutate(land_cov = gsub("pct.", "", land_cov)) %>%
  mutate(land_cov_pretty = recode(land_cov,
                           'wetland' = 'Total Wetland',
                           'nat' = 'Natural',
                           'human' = 'Anthropogenic',
                           'ww' = 'Woody Wetland',
                           'wwnm' = 'Woody Wetland No Mangrove',
                           'ew' = 'Emergent Wetland',
                           'ag' = 'Agriculture',
                           'bare' = 'Barren',
                           'wat' = 'Water',
                           'grass' = 'Grassland',
                           'man' = 'Mangrove',
                           'for' = 'Forest',
                           'ur' = 'Urban')) %>%
  mutate(pct = round(pct, digits = 2)) %>%
  filter(land_cov != "wwnm")
  

ggplot(data = hist.df, aes(pct)) +
  geom_histogram() +
  facet_wrap(~land_cov_pretty) + 
  theme_bw() + 
  xlab("% Cover") + 
  ylab("Count") + 
  theme(axis.title = element_text(family = "serif", size = 12),
        axis.text = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 12),
        strip.background = element_blank()) + 
  scale_x_continuous(labels = prettyZero)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/SupplementalInfo/LULCpct.jpg"), device = "jpeg",
       dpi = 600, width = 8, height = 8, units = "in")
  

#Climate and Land Cover Trends over time
temp.t <- inla(tmean ~ Year + f(Site, model = "iid"), data = df)
summary(temp.t)

temp.max <- inla(tmax ~ Year + f(Site, model = "iid"), data = df)
summary(temp.max)

ppt.t <- inla(ppt ~ Year + f(Site, model = "iid"), data = df)
summary(ppt.t)

man.t <- inla(pct.man ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(man.t)

ww.t <- inla(pct.ww ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ww.t)

ew.t <- inla(pct.ew ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ew.t)

ur.t <- inla(pct.ur ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ur.t)

##################################################
#########Species Richness Regression Models####### 
##################################################

sr1.ar <- inla(alpha ~ tmean +  ppt + pct.man + 
                 pct.ur + pct.ag + pct.wwnm + pct.ew + 
                 f(Site, model = "iid") + 
                 f(Route, model = "iid"), 
                 #f(BCR, model = "iid") +
                 #f(Year, model = "ar1"), 
                 data = model.df)
summary(sr1.ar)


#################################################
#############Beta Diversity Models###############
#################################################
hist.diff.df <- df %>%
  select(Route, Site, Year, Yr_bin, BCR, alpha, alpha.full, beta, beta.full, turn, nest, beta, 
         tmean.anom, tmin.anom, tmax.anom, ppt.anom, 
         diff.man, diff.ur, diff.wwnm, diff.ew, diff.nat, diff.wet, diff.grass, diff.ur, diff.human, 
         diff.for, diff.wat, diff.ag)

hist.df <- df %>%
  select(matches("diff")) %>%
  pivot_longer(cols = matches("diff."), names_to = "land_cov", values_to = "diff.pct") %>%
  mutate(land_cov = gsub("diff.", "", land_cov)) %>%
  mutate(land_cov_pretty = recode(land_cov,
                                  'wet' = 'Total Wetland',
                                  'nat' = 'Natural',
                                  'human' = 'Anthropogenic',
                                  'ww' = 'Woody Wetland',
                                  'wwnm' = 'Woody Wetland No Mangrove',
                                  'ew' = 'Emergent Wetland',
                                  'ag' = 'Agriculture',
                                  'bare' = 'Barren',
                                  'wat' = 'Water',
                                  'grass' = 'Grassland',
                                  'man' = 'Mangrove',
                                  'for' = 'Forest',
                                  'ur' = 'Urban')) %>%
  mutate(pct = round(diff.pct, digits = 2)) %>%
  filter(land_cov != "wwnm")


ggplot(data = hist.df, aes(diff.pct)) +
  geom_histogram() +
  facet_wrap(~land_cov_pretty) + 
  theme_bw() + 
  xlab(expression(Delta*" % Cover")) + 
  ylab("Count") + 
  theme(axis.title = element_text(family = "serif", size = 12),
        axis.text = element_text(family = "serif", size = 12),
        strip.text = element_text(family = "serif", size = 12),
        strip.background = element_blank()) + 
  scale_x_continuous(labels = prettyZero)

ggsave(here::here("Figures/Figures_Diversity_Manuscript/SupplementalInfo/LULCpctdiff.jpg"), device = "jpeg",
       dpi = 600, width = 8, height = 8, units = "in")



beta.df <- df %>%
  select(Route, Site, Year, Yr_bin, BCR, alpha, alpha.full, beta, beta.full, turn, nest, beta, 
         tmean.anom, tmin.anom, tmax.anom, ppt.anom, 
         diff.man, diff.ur, diff.wwnm, diff.ew, diff.nat, diff.wet, diff.grass, diff.ur, diff.human, 
         diff.for, diff.wat, diff.ag) %>%
  mutate(tmean.anom.mag = abs(tmean.anom), tmax.anom.mag = abs(tmax.anom), ppt.anom.mag = abs(ppt.anom), 
         diff.man.mag = abs(diff.man), diff.ur.mag = abs(diff.ur), diff.wwnm.mag = abs(diff.wwnm), diff.ew.mag = abs(diff.ew),
         diff.nat.mag = abs(diff.nat), diff.wet.mag = abs(diff.wet), diff.grass.mag = abs(diff.grass), diff.ur.mag = abs(diff.ur), 
         diff.human.mag = abs(diff.human), diff.ag.mag = abs(diff.ag), diff.for.mag = abs(diff.for), diff.wat.mag = abs(diff.wat)) %>%
  mutate_at(c("tmean.anom", "tmax.anom", "ppt.anom", "diff.man", "diff.ur", "diff.wwnm", "diff.ew", 
              "diff.nat", "diff.wet", "diff.grass", "diff.ur", "diff.human", "diff.ag", 
              "diff.for", "diff.wat", "tmean.anom.mag", "tmax.anom.mag", "ppt.anom.mag", "diff.man.mag", 
              "diff.ur.mag", "diff.wwnm.mag", "diff.ew.mag", "diff.nat.mag", "diff.wet.mag", "diff.grass.mag", "diff.ur.mag", 
              "diff.human.mag", "diff.ag.mag", "diff.for.mag", "diff.wat.mag"), scale)


#Transform 0s and 1s in the data 
beta.trans <- function(y, n) {
  (y * (n - 1) + 0.5)/n
}

ss <- beta.df %>% 
  filter(!is.na(beta)) %>%
  nrow()

beta.df <- beta.df %>%
  mutate(beta = if_else(beta == 0, beta.trans(beta, ss), beta),
         turn = if_else(turn == 0, beta.trans(turn, ss), turn),
         nest = if_else(nest == 0, beta.trans(nest, ss), nest),
         beta.full = if_else(beta.full == 0, beta.trans(beta.full, ss), beta.full))


#Checking out the beta diversity covariates 
hist(beta.df$diff.man, breaks = 20)  
hist(beta.df$diff.ur, breaks = 20)  
hist(beta.df$diff.wwnm, breaks = 20)  
hist(beta.df$diff.ew, breaks = 20)
hist(beta.df$diff.ag, breaks = 20)
hist(beta.df$diff.grass, breaks = 20)
hist(beta.df$diff.wat, breaks = 20)
hist(beta.df$diff.for, breaks = 20)
hist(beta.df$diff.bare, breaks = 20)

#Correlograms for the beta diversity predictors 
beta.pred <- df %>%
  select(tmean.anom, ppt.anom,
         diff.man, diff.ww, diff.ur, diff.wwnm, 
         diff.ew, diff.wet, diff.human, diff.ag, 
         diff.wat, diff.bare, diff.for, diff.grass, diff.nat)

cors <- cor(beta.pred)
cor.df <- as.data.frame(round(cors, 2))

pacman::p_load(corrplot)
corrplot(cors, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, method = "number")


#Beta Model 
beta.df$BCR1 <- beta.df$BCR2 <- beta.df$BCR3 <- beta.df$BCR4 <- beta.df$BCR5 <- beta.df$BCR6 <- beta.df$BCR7 <- beta.df$BCR

formula1 <- beta ~ 1 + tmean.anom + ppt.anom + diff.man + diff.ur + diff.ew + diff.wwnm + alpha + 
  f(BCR1, tmean.anom) + f(BCR2, ppt.anom) + f(BCR3, diff.man) + 
  f(BCR4, diff.ur) + f(BCR5, diff.ew) + f(BCR6, diff.wwnm) + 
  f(Site, model = "iid") + f(Route, model = "iid") + f(BCR, model = "iid") +
  f(Year, model = "ar1")

#####################
######Best model#####
#####################

formula2 <- beta ~ 1 + tmean.anom + ppt.anom + diff.man + diff.ag + diff.ur + diff.ew + diff.wwnm + alpha + 
  f(Site, model = "iid") + f(Route, model = "iid") +
  f(Year, model = "ar1")

#####################

formula3 <- beta ~ 1 + tmean.anom + ppt.anom + diff.man + diff.ur + diff.ew + diff.wwnm + alpha + 
  f(Site, model = "iid") + f(Route, model = "iid") +
  f(Year, model = "iid")

formula4 <- beta ~ 1 + tmean.anom + ppt.anom + diff.man + diff.ur + diff.ew + diff.wwnm + alpha + 
  f(Site, model = "iid") + f(Route, model = "iid") +
  f(Year, model = "rw1")


beta1.re <- inla(formula1,
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(beta1.re)

beta2.re <- inla(formula2,
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(beta2.re)

beta3.re <- inla(formula3,
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(beta3.re)

beta4.re <- inla(formula4,
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(beta4.re)


c(beta1.re$dic$dic, beta2.re$dic$dic, beta3.re$dic$dic, beta4.re$dic$dic)

for(i in 1:5){
m <- beta2.re$internal.marginals.hyperpar[[i]]
m.var <- inla.tmarginal(function(x) 1/exp(x), m)
print(inla.zmarginal(m.var))
}
#Full beta model

summary(beta1.re)

marginals_mode <- sapply(beta1.re$marginals.hyperpar,
                         function(x)
                           inla.mmarginal(inla.tmarginal(function(x) 1/x, x)))
names(marginals_mode) <- sapply(as.vector(as.character(names(marginals_mode2))),
                                function(y) gsub("Precision", x=y, "Mode of variance"))
marginals_mode

marginals_mode2 <- sapply(beta2.re$marginals.hyperpar,
                         function(x)
                           inla.mmarginal(inla.tmarginal(function(x) 1/x, x)))
names(marginals_mode2) <- sapply(as.vector(as.character(names(marginals_mode2))),
                                function(y) gsub("Precision", x=y, "Mode of variance"))
marginals_mode2


beta1.re.mag <- inla(beta ~ tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                   diff.ur.mag + diff.for.mag + diff.grass.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                   f(Site, model = "iid") + f(Route, model = "iid") + 
                   f(Year, model = "ar1"),
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)


summary(beta1.re.mag)


#Turnover models 
turn1.re <- inla(turn ~  1 + tmean.anom + ppt.anom + diff.man + 
                   diff.ur + diff.wwnm + diff.ew + alpha + 
                   f(Site, model = "iid") + f(Route, model = "iid") + 
                   f(Year, model = "ar1"),
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(turn1.re)

for(i in 1:5){
m <- nest1.re$internal.marginals.hyperpar[[i]]
m.var <- inla.tmarginal(function(x) 1/exp(x), m)
print(inla.zmarginal(m.var))
}


marginals_mode <- sapply(turn1.re$marginals.hyperpar,
                         function(x)
                           inla.mmarginal(inla.tmarginal(function(x) 1/x, x)))
names(marginals_mode) <- sapply(as.vector(as.character(names(marginals_mode))),
                                function(y) gsub("Precision", x=y, "Mode of variance"))
marginals_mode


turn1.re.mag <- inla(turn ~  tmax.anom.mag + ppt.anom.mag + diff.man.mag + 
                       diff.ur.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                       f(Site, model = "iid") + f(Route, model = "iid") + 
                       f(Year, model = "ar1"),
                     data = beta.df, family = "beta", 
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     verbose = F)
summary(turn1.re.mag)


#Nestedness models 
nest1.re <- inla(nest ~ 1 + tmean.anom + ppt.anom + diff.man + 
                   diff.ur + diff.wwnm + diff.ew + alpha + 
                   f(Site, model = "iid") + f(Route, model = "iid") + 
                   f(Year, model = "ar1"),
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(nest1.re)

marginals_mode_nest <- sapply(nest1.re$marginals.hyperpar,
                         function(x)
                           inla.mmarginal(inla.tmarginal(function(x) 1/x, x)))
names(marginals_mode) <- sapply(as.vector(as.character(names(marginals_mode))),
                                function(y) gsub("Precision", x=y, "Mode of variance"))
marginals_mode

nest1.re.mag <- inla(nest ~  tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                   diff.ur.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                   f(Site, model = "iid") + f(Route, model = "iid") + 
                   f(Year, model = "ar1"),
                 data = beta.df, family = "beta", 
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                 verbose = F)
summary(nest1.re)


#Data vizualization
if(!require(ggregplot)) devtools::install_github("gfalbery/ggregplot")
library('ggregplot', 'MCMCglmm')

Efxplot(list(sr1.ar))
Efxplot(list(beta2.re))
Efxplot(list(turn1.re))
Efxplot(list(nest1.re))

#Partial Prediction Plots
beta.df.pp <- df %>%
  select(Route, Site, Year, Yr_bin, BCR, alpha, alpha.full, beta, beta.full, turn, nest, beta, 
         tmean.anom, tmin.anom, tmax.anom, ppt.anom, 
         diff.man, diff.ur, diff.wwnm, diff.ew, diff.nat, diff.wet, diff.grass, diff.ur, diff.human, 
         diff.for, diff.wat, diff.ag) %>%
  mutate(tmean.anom.mag = abs(tmean.anom), tmax.anom.mag = abs(tmax.anom), ppt.anom.mag = abs(ppt.anom), 
         diff.man.mag = abs(diff.man), diff.ur.mag = abs(diff.ur), diff.wwnm.mag = abs(diff.wwnm), diff.ew.mag = abs(diff.ew),
         diff.nat.mag = abs(diff.nat), diff.wet.mag = abs(diff.wet), diff.grass.mag = abs(diff.grass), diff.ur.mag = abs(diff.ur), 
         diff.human.mag = abs(diff.human), diff.ag.mag = abs(diff.ag), diff.for.mag = abs(diff.for), diff.wat.mag = abs(diff.wat))

#Transform 0s and 1s in the data 
ss <- beta.df.pp %>% 
  filter(!is.na(beta)) %>%
  nrow()

beta.df.pp <- beta.df.pp %>%
  mutate(beta = if_else(beta == 0, beta.trans(beta, ss), beta),
         turn = if_else(turn == 0, beta.trans(turn, ss), turn),
         nest = if_else(nest == 0, beta.trans(nest, ss), nest),
         beta.full = if_else(beta.full == 0, beta.trans(beta.full, ss), beta.full))


pred.df <- beta.df.pp %>%
  dplyr::select(beta, tmean.anom.mag, ppt.anom.mag, diff.man.mag, diff.for.mag, diff.grass.mag,  
                diff.ur.mag, diff.ag.mag, diff.wwnm.mag, diff.ew.mag, alpha,
                Site, Route, Year)

  
pred.rows <- data.frame(beta = NA, tmean.anom.mag = mean(pred.df$tmean.anom.mag), ppt.anom.mag = mean(pred.df$ppt.anom.mag),
                        diff.ew.mag = seq(min(pred.df$diff.ew.mag), max(pred.df$diff.ew.mag), len = 100), 
                        diff.ag.mag = mean(pred.df$diff.ag.mag), 
                        diff.wwnm.mag = mean(pred.df$diff.wwnm.mag), diff.ur.mag = mean(pred.df$diff.ur.mag),
                        diff.man.mag = mean(pred.df$diff.man.mag), diff.for.mag = mean(pred.df$diff.for.mag), 
                        diff.grass.mag = mean(pred.df$diff.grass.mag), alpha = mean(pred.df$alpha, na.rm = T),
                        Site = NA, Route = NA, Year = NA)
pred.test <- rbind(pred.df, pred.rows)  

beta1.re.mag <- inla(beta ~  tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                       diff.ur.mag + diff.for.mag + diff.grass.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                       f(Site, model = "iid") + f(Route, model = "iid") + 
                       f(Year, model = "ar1"),
                     data = pred.test, family = "beta", 
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     verbose = F)

#summary(beta1.re.mag)

newdata <- cbind(pred.rows,
                 beta1.re.mag$summary.fitted.values[(nrow(pred.df)+1):nrow(pred.test),])
head(newdata)

newdata <- reshape:::rename(newdata, c("0.025quant"="lower", "0.975quant"="upper"))

library(grid)
ggplot(newdata, aes(y=expit(mean), x=diff.ew.mag)) +
  geom_blank()+
  geom_point(data=pred.df, aes(y = beta, x=diff.ew.mag)) +
  geom_ribbon(aes(ymin=expit(lower), ymax=expit(upper)), fill='blue', alpha=0.2) +
  geom_line() +
  scale_y_continuous(expression(Wetland~Bird~beta-diversity))+
  xlab(expression(Magnitude~Delta~of~Emergent~Wetland~Cover)) +
  theme_classic()+
  theme(legend.key.width=unit(2,'lines'), legend.position=c(0,1), legend.justification=c(0,1))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/PPplotEWnodat.jpeg"), device = "jpeg",
       width = 6, height = 6, units = "in", dpi = 400)

###Urban PP plots 
pred.rows <- data.frame(beta = NA, tmean.anom.mag = mean(pred.df$tmean.anom.mag), ppt.anom.mag = mean(pred.df$ppt.anom.mag),
                        diff.ur.mag = seq(min(pred.df$diff.ur.mag), max(pred.df$diff.ur.mag), len = 100), 
                        diff.ag.mag = mean(pred.df$diff.ag.mag), 
                        diff.wwnm.mag = mean(pred.df$diff.wwnm.mag), diff.ew.mag = mean(pred.df$diff.ew.mag),
                        diff.man.mag = mean(pred.df$diff.man.mag), diff.for.mag = mean(pred.df$diff.for.mag), 
                        diff.grass.mag = mean(pred.df$diff.grass.mag), alpha = mean(pred.df$alpha, na.rm = T),
                        Site = NA, Route = NA, Year = NA)
pred.test <- rbind(pred.df, pred.rows)  

beta1.re.mag <- inla(beta ~  tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                       diff.ur.mag + diff.for.mag + diff.grass.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                       f(Site, model = "iid") + f(Route, model = "iid") + 
                       f(Year, model = "ar1"),
                     data = pred.test, family = "beta", 
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     verbose = F)

#summary(beta1.re.mag)

newdata <- cbind(pred.rows,
                 beta1.re.mag$summary.fitted.values[(nrow(pred.df)+1):nrow(pred.test),])
head(newdata)

newdata <- reshape:::rename(newdata, c("0.025quant"="lower", "0.975quant"="upper"))

library(grid)
ggplot(newdata, aes(y=expit(mean), x=diff.ur.mag)) +
  geom_blank()+
  #geom_point(data=pred.df, aes(y=beta, x=diff.ew.mag)) +
  geom_ribbon(aes(ymin=expit(lower), ymax=expit(upper)), fill='blue', alpha=0.2) +
  geom_line() +
  scale_y_continuous(expression(Wetland~Bird~beta-diversity))+
  xlab(expression(Magnitude~Delta~of~Urban~Cover)) +
  theme_classic()+
  theme(legend.key.width=unit(2,'lines'), legend.position=c(0,1), legend.justification=c(0,1))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/PPplotURnodat.jpeg"), device = "jpeg",
       width = 6, height = 6, units = "in", dpi = 400)

#Alpha PP plots
pred.rows <- data.frame(beta = NA, tmean.anom.mag = mean(pred.df$tmean.anom.mag), ppt.anom.mag = mean(pred.df$ppt.anom.mag),
                        alpha = seq(min(pred.df$alpha, na.rm = T), max(pred.df$alpha, na.rm = T), len = 100), 
                        diff.ag.mag = mean(pred.df$diff.ag.mag), 
                        diff.wwnm.mag = mean(pred.df$diff.wwnm.mag), diff.ew.mag = mean(pred.df$diff.ew.mag),
                        diff.man.mag = mean(pred.df$diff.man.mag), diff.for.mag = mean(pred.df$diff.for.mag), 
                        diff.grass.mag = mean(pred.df$diff.grass.mag), diff.ur.mag = mean(pred.df$diff.ur.mag),
                        Site = NA, Route = NA, Year = NA)
pred.test <- rbind(pred.df, pred.rows)  

beta1.re.mag <- inla(beta ~  tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                       diff.ur.mag + diff.for.mag + diff.grass.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                       f(Site, model = "iid") + f(Route, model = "iid") + 
                       f(Year, model = "ar1"),
                     data = pred.test, family = "beta", 
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     verbose = F)

#summary(beta1.re.mag)

newdata <- cbind(pred.rows,
                 beta1.re.mag$summary.fitted.values[(nrow(pred.df)+1):nrow(pred.test),])
head(newdata)

newdata <- reshape:::rename(newdata, c("0.025quant"="lower", "0.975quant"="upper"))

library(grid)
ggplot(newdata, aes(y=expit(mean), x=alpha)) +
  geom_blank()+
  #geom_point(data=pred.df, aes(y=beta, x=diff.ew.mag)) +
  geom_ribbon(aes(ymin=expit(lower), ymax=expit(upper)), fill='blue', alpha=0.2) +
  geom_line() +
  scale_y_continuous(expression(Wetland~Bird~beta~diversity))+
  xlab(expression(alpha~diversity)) +
  theme_classic()+
  theme(legend.key.width=unit(2,'lines'), legend.position=c(0,1), legend.justification=c(0,1))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/PPplotAlphanodat.jpeg"), device = "jpeg",
       width = 6, height = 6, units = "in", dpi = 400)

#PPT PP plots
pred.rows <- data.frame(beta = NA, tmean.anom.mag = mean(pred.df$tmean.anom.mag), ppt.anom.mag = seq(min(pred.df$ppt.anom.mag), max(pred.df$ppt.anom.mag, len = 100)),
                        alpha = mean(pred.df$alpha, na.rm = T), 
                        diff.ag.mag = mean(pred.df$diff.ag.mag), 
                        diff.wwnm.mag = mean(pred.df$diff.wwnm.mag), diff.ew.mag = mean(pred.df$diff.ew.mag),
                        diff.man.mag = mean(pred.df$diff.man.mag), diff.for.mag = mean(pred.df$diff.for.mag), 
                        diff.grass.mag = mean(pred.df$diff.grass.mag), diff.ur.mag = mean(pred.df$diff.ur.mag),
                        Site = NA, Route = NA, Year = NA)
pred.test <- rbind(pred.df, pred.rows)  

beta1.re.mag <- inla(beta ~  tmean.anom.mag + ppt.anom.mag + diff.man.mag + 
                       diff.ur.mag + diff.for.mag + diff.grass.mag + diff.wwnm.mag + diff.ew.mag + alpha + 
                       f(Site, model = "iid") + f(Route, model = "iid") + 
                       f(Year, model = "ar1"),
                     data = pred.test, family = "beta", 
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                     verbose = F)

#summary(beta1.re.mag)

newdata <- cbind(pred.rows,
                 beta1.re.mag$summary.fitted.values[(nrow(pred.df)+1):nrow(pred.test),])
head(newdata)

newdata <- reshape:::rename(newdata, c("0.025quant"="lower", "0.975quant"="upper"))

library(grid)
ggplot(newdata, aes(y=expit(mean), x=ppt.anom.mag)) +
  geom_blank()+
  #geom_point(data=pred.df, aes(y=beta, x=diff.ew.mag)) +
  geom_ribbon(aes(ymin=expit(lower), ymax=expit(upper)), fill='blue', alpha=0.2) +
  geom_line() +
  scale_y_continuous(expression(Wetland~Bird~beta~diversity))+
  xlab(expression(Delta~Magnitude~Precipitation)) +
  theme_classic()+
  theme(legend.key.width=unit(2,'lines'), legend.position=c(0,1), legend.justification=c(0,1))

ggsave(here::here("Figures/Figures_Diversity_Manuscript/PPplotPPTnodat.jpeg"), device = "jpeg",
       width = 6, height = 6, units = "in", dpi = 400)


#df.pred.plots
lulc.vars <- c("pct.man", "pct.wwnm", "pct.ew", "pct.for", "pct.grass",
               "pct.ag", "pct.ur", "pct.wat")

lulc.names <- list(
  "pct.man" = "% Mangrove",
  "pct.wwnm" = "% Woody Wetland",
  "pct.ew" = "% Emergent Wetland",
  "pct.for" = "% Forest",
  "pct.grass" = "% Grassland",
  "pct.ag" = "% Agriculture",
  "pct.ur" = "% Urban",
  #"pct.bare" = "% Bare",
  "pct.wat" = "% Water"
)
lulc_labeller <- function(variable,value){
  return(lulc.names[value])
}

years <- c(1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015)
Yr_bin <- 1:8
yr_df <- data.frame(Year = years, Yr_bin = Yr_bin)

df.pp <- df %>%
  select(BCR_name, Yr_bin, pct.man, pct.wwnm, pct.ew, pct.for, pct.grass,
         pct.ag, pct.ur, pct.wat) %>%
  left_join(yr_df, by = "Yr_bin") %>%
  pivot_longer(cols = lulc.vars, names_to = "LULC", values_to = "pct")

plot_lulc <- ggplot(df.pp, aes(x = Year, y = pct, color = BCR_name)) +
  geom_smooth(method = 'loess', span = 1) +
  facet_wrap(vars(LULC), nrow = 2, labeller = lulc_labeller) + 
  theme_bw() +
  scale_color_viridis_d(option = "cividis") +
  theme(strip.background = element_blank(),
        strip.text = element_text(family = "serif", size = 12),
        legend.text = element_text(family = "serif", size = 12),
        axis.text = element_text(family = "serif", size = 12),
        legend.title = element_blank(),
        axis.title = element_text(family = "serif", size = 12)) +
  xlab("Year") + 
  ylab("Percent Cover")
  
ggsave(filename = here::here("Figures/LULC_Trends.tiff"), plot = plot_lulc, height = 8, width = 10, units = "in",
       device = "tiff", dpi = 600)
  

#New Beta Diversity - Time plots 
plot_beta <- ggplot(df, aes(x = Year, y = pct.ag, group = BCR, 
                                 colour = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                                 "Southeastern Coastal Plain", 
                                                                 "Peninsular Florida", 
                                                                 "Gulf Coastal Prairie")))) +
                #geom_jitter(size = 1, color = "grey") +
                #geom_line()+
                geom_smooth(method = loess, se = T, aes(x = Year, y = pct.ag, group = BCR,
                                                        color = factor(BCR, labels = c("Mississippi Alluvial Valley", 
                                                                                        "Southeastern Coastal Plain", 
                                                                                        "Peninsular Florida", 
                                                                                        "Gulf Coastal Prairie")))) +
                #scale_color_brewer(palette = )
                xlab("Year") +
                ylab(expression(paste(beta, "-diversity"))) +
                labs(colour = "Bird Conservation Region") +
                #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 20),
                      axis.text.y = element_text(size = 20),
                      axis.title.x = element_text(vjust = -1, size = 14),
                      axis.title.y = element_text(vjust = 1.5, size = 14),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=20)) 

plot_beta <- plot_beta + scale_color_viridis(discrete = T)

ggdraw() + 
  draw_plot(plot_beta + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.65, y = .025, width = .2, height = .75, scale = 1) 

