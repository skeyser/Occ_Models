#########################################################################################
#######################Script for turning MSOM Output back into##########################
#######################Community Data for Biodiversity Analyses##########################
#########################################################################################
##########################Script Created By: Spencer R Keyser############################
##################################### 7/16/2019 #########################################
#########################################################################################

#Clean the environment
rm(list = ls())

#Load in Packages

library("pacman")
pacman::p_load("tidyverse", "plyr", "reshape2", "here", "arrayhelpers")

#Finish Package Loading 

#Load in Work spaces of z-matrices 
load(file = here::here("R Workspace/Output/Output4Analysis/Muscicapida.RData"))
jdf <- read.csv(file = here::here("Data_BBS/Generated DFs/jdf.csv"))
jdf <- jdf[, -1]
yrs <- read.csv(file = here::here("Data_BBs/Generated DFs/Years.Occ.csv"))
JKmat <- readRDS(file = here::here("Data_BBS/Generated DFs/JKmat.RDS"))
JKsurv <- readRDS(file = here::here("Data_BBS/Generated DFs/JKsurv.RDS"))

#Checkout the z.prime matrix 
dim(z.prime)

#Put all of the surveys that were not conducted into the z.prime
S <- max(spp.occ$spp.id)
J <- 274
K <- 38

for(i in 1:S){
  z.prime[ i, , ] <- z.prime[ i, , ] * JKmat
}

#Create a z.prime.95 based on other literature just for keeping if need be
#Reworking all of the z.primes so the cutoffs include the value too (aka >= not >)
z.prime.5 <- means.output$z
z.prime.65 <- means.output$z
z.prime.75 <- means.output$z
z.prime.95 <- means.output$z


#Change values in matrix 
for (s in 1:S){
  for (j in 1:J){
    for (k in 1:K){
      z.prime.5[s, j, k] <- ifelse(z.prime.5[s, j, k] >= 0.5, 1, 0)
      z.prime.65[s, j, k] <- ifelse(z.prime.65[s, j, k] >= 0.65, 1, 0)
      z.prime.75[s, j, k] <- ifelse(z.prime.75[s, j, k] >= 0.75, 1, 0)
      z.prime.95[s, j, k] <- ifelse(z.prime.95[s, j, k] >= 0.95, 1, 0)
    }
  }  
}

#Quick function to calculate total observations of each species
total <- function(x, y){
  df <- x[y,,] 
  #did I just break your function by using JKmat instead?
  yr.sums <- colSums(df, 1:length(nrow), na.rm = TRUE)
  sum(yr.sums, na.rm = TRUE)
}

#Adding in the z.prime.95 into the DF
total.observed$Estimated.0.95 <- NA

for (i in 1:S){
  total.observed[i, 5] <- total(z.prime.95, i)
}

#z.prime is spp x sites x years
#Create empty matrix for handling 2-D matrices from z.prime

z2d.prime <- array2df(z.prime, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.5 <- array2df(z.prime.5, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.65 <- array2df(z.prime.65, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.75 <- array2df(z.prime.75, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.95 <- array2df(z.prime.95, levles = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")

#UDF if col[x] is a factor make it an integer, if not make it a numeric
no.factors <- function(x){
  if(is.factor(x) == T){
    as.integer(x)
  } else {
    as.numeric(x)
  }
}

#Apply the UDF (no.factors)to all the DFs
z2d.prime[] <- lapply(z2d.prime, no.factors)
z2d.prime.5[] <- lapply(z2d.prime.5, no.factors)
z2d.prime.65[] <- lapply(z2d.prime.65, no.factors)
z2d.prime.75[] <- lapply(z2d.prime.75, no.factors)
z2d.prime.95[] <- lapply(z2d.prime.95, no.factors)



#Combine z.dfs with info from spp.occ and jdf
df.prime <- z2d.prime %>% left_join( yrs, by = "year.id") %>%
            left_join(jdf, by = "site.id") %>% 
            left_join(spp.occ, by = "spp.id")

df.prime.5 <- z2d.prime.5 %>% left_join( yrs, by = "year.id") %>%
              left_join(jdf, by = "site.id") %>% 
              left_join(spp.occ, by = "spp.id")

df.prime.65 <- z2d.prime.65 %>% left_join( yrs, by = "year.id") %>%
               left_join(jdf, by = "site.id") %>% 
               left_join(spp.occ, by = "spp.id")

df.prime.75 <- z2d.prime.75 %>% left_join( yrs, by = "year.id") %>%
               left_join(jdf, by = "site.id") %>% 
               left_join(spp.occ, by = "spp.id")

df.prime.95 <- z2d.prime.95 %>% left_join( yrs, by = "year.id") %>% 
               left_join(jdf, by = "site.id") %>% 
               left_join(spp.occ, by = "spp.id")



group.name <- as.character(unique(spp.occ$Phylo.V1))
group.name1 <- paste0(group.name, "_", "CommunityDF.csv")
group.name.5 <- paste0(group.name, "_", "CommunityDF5.csv")
group.name.65 <- paste0(group.name, "_", "CommunityDF65.csv")
group.name.75 <- paste0(group.name, "_", "CommunityDF75.csv")


write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name1)))
write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.5)))
write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.65)))
write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.75)))

