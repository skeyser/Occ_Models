##############################################################
#############INLA Analysis Mangrove###########################
###############Created By: Spencer R Keyser###################
####################8/18/2020#################################

#Install packages
if ("INLA" %in% rownames(installed.packages()) == FALSE){
  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

#Load Packages 
pacman::p_load("tidyverse", "MASS", "INLA")
select <- dplyr::select
#Bring in data 
df <- read.csv(here::here("Data_BBS/Generated DFs/INLAdfFinal.csv"))

df <- df %>%
  group_by(Site) %>%
  mutate(SR.change = (SR50 - first(SR50))/first(SR50)) 

#Look at correlations in the data 
df.pred <- df %>%
  select(diff.man, diff.ww, diff.wwnm, diff.ew, diff.human, diff.ag, diff.wat, diff.bare, diff.for, diff.grass, diff.nat, diff.wet)

cors <- cor(df.pred)
cor.df <- as.data.frame(round(cors, 2))

pacman::p_load(corrplot)
corrplot(cors, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

#Scale Covariates 
#DF created for running first model for species richness
model.df <- df %>%
  select(Rteno, Site, Year, Yr_bin, SR50, SR.change, SR.wet, beta.wet.jac, beta.wet.nest, beta.wet.turn, beta50.jac, tmean = tmean.bird.c, ppt.dry = precip_Dry.c, 
         ppt.wet = precip_Wet.c, pct.man, pct.ur, pct.nat) %>%
  mutate_at(c("tmean", "ppt.dry", "ppt.wet", "pct.man", "pct.ur", "pct.nat"), scale)
  
#Check out the distribution of covariates and response variables
hist(df$pct.man, breaks = 20)  

#Climate and Land Cover Trends over time
temp.t <- inla(tmean.bird.c ~ Year + f(Site, model = "iid"), data = df)
summary(temp.t)

ppt.t <- inla(precip_Wet.c ~ Year + f(Site, model = "iid"), data = df)
summary(ppt.t)

pptdry.t <- inla(precip_Dry.c ~ Year + f(Site, model = "iid"), data = df)
summary(pptdry.t)

man.t <- inla(pct.man ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(man.t)

ww.t <- inla(pct.ww ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ww.t)

ew.t <- inla(pct.ew ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ew.t)

ur.t <- inla(pct.ur ~ Yr_bin + f(Site, model = "iid"), data = df)
summary(ur.t)

#Species Richness Regression
sr1 <- inla(SR50 ~ tmean + ppt.dry + ppt.wet + pct.man + pct.nat + pct.ur, data = model.df)
summary(sr1)

sr1.me <- inla(SR50 ~ tmean + ppt.wet + pct.man + pct.nat + pct.ur + 
                 f(Site, model = "iid") + f(Rteno, model = "iid"), data = model.df)
summary(sr1.me)

sr1.ar <- inla(SR50 ~ tmean + ppt.wet + pct.man + pct.nat + pct.ur + 
                 f(Site, model = "iid") + f(Rteno, model = "iid") + 
                 f(Year, model = "ar1"), data = model.df)
summary(sr1.ar)

srwet.ar <- inla(SR.wet ~ tmean + ppt.wet + pct.man + pct.nat + pct.ur + 
                 f(Site, model = "iid") + f(Rteno, model = "iid") + 
                 f(Year, model = "ar1"), data = model.df)
summary(srwet.ar)

#Breeding season temp and precip
#Region effect
#RW1
#Model Comparison 
#Logit tranform beta and use a normal 

#Partial Prediction

