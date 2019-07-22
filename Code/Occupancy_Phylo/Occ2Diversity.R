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

#library("pacman")
pacman::p_load("tidyverse", "plyr", "reshape2", "here", "arrayhelpers")

#Finish Package Loading 

#Load in Work spaces of z-matrices 
load(file = here::here("R Workspace/Output/Output4Analysis/Tits.RData"))
jdf <- read.csv(file = here::here("Data_BBS/Generated DFs/jdf.csv"))
jdf <- jdf[, -1]
yrs <- read.csv(file = here::here("Data_BBs/Generated DFs/Years.Occ.csv"))
JKmat <- readRDS(file = here::here("Data_BBS/Generated DFs/JKmat.RDS"))
JKsurv <- readRDS(file = here::here("Data_BBS/Generated DFs/JKsurv.RDS"))

#If spp.occ df doesn't exist for the given work space it needs to be read in
#saveRDS(spp.occ, file = here::here("Data_BBS/Generated DFs/Emberizidae_spp_occ.RDS"))
#spp.occ <- readRDS(here::here("Data_BBS/Generated DFs/Emberizidae_spp_occ.RDS"))

#Checkout the z.prime matrix
dim(z.prime)

#Put all of the surveys that were not conducted into the z.prime
S <- max(spp.occ$spp.id)
J <- 274
K <- 38

z.prime.95 <- means.output$z


for(i in 1:S){
  z.prime[ i, , ] <- z.prime[ i, , ] * JKmat
  z.prime.95[ i, , ] <- z.prime.95[ i, , ] * JKmat
}

#Create a z.prime.95 based on other literature just for keeping if need be
#Reworking all of the z.primes so the cutoffs include the value too (aka >= not >)

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
  #total.observed[i, 5] <- total(z.prime.95, i)
  total.observed[i, 7] <- total(z.prime.95, i)
}


#z.prime is spp x sites x years
#Create empty matrix for handling 2-D matrices from z.prime

z2d.prime <- array2df(z.prime, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.5 <- array2df(z.prime.5, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.65 <- array2df(z.prime.65, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.75 <- array2df(z.prime.75, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.95 <- array2df(z.prime.95, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")

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


#Add the z.prime.95 into the alpha.diversity calculations
group.a.div.95 <- matrix(NA, J, K)

for (j in 1:J){
  for ( k in  1:K){
    #zeros those z for unsampled segments
    group.a.div.95[ j, k ] <- sum( z.prime.95[ 1:S, j, k ], na.rm = TRUE)# * JKsurv[ j, k ] )
    #keep only average estimate
    #mean.a.div[ j, k ] <- mean(a.div[ j, k ])
  }#K
}#J

alpha.div$group.a.div.95 <- group.a.div.95 


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
group.name <- gsub("/", "_", group.name)
group.name1 <- paste0(group.name, "_", "CommunityDF.csv")
group.name.5 <- paste0(group.name, "_", "CommunityDF5.csv")
group.name.65 <- paste0(group.name, "_", "CommunityDF65.csv")
group.name.75 <- paste0(group.name, "_", "CommunityDF75.csv")
group.name.95 <- paste0(group.name, "_", "CommunityDF95.csv")


write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name1)))
write.csv(df.prime.5, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.5)))
write.csv(df.prime.65, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.65)))
write.csv(df.prime.75, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.75)))
write.csv(df.prime.95, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name.95)))

image.name <- paste0(group.name, "", ".RData")
save.image(file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", image.name)))

####Read in all of the new community matrices####


#List files based on path names 
comm_list <- list.files(path = here::here("Data_BBS/Generated DFs/OccMod_CommMats"), pattern = "*DF.csv", full.names = T)
comm_list5 <- list.files(path = here::here("Data_BBS/Generated DFs/OccMod_CommMats"), pattern = "*DF5.csv", full.names = T)
comm_list65 <- list.files(path = here::here("Data_BBS/Generated DFs/OccMod_CommMats"), pattern = "*DF65.csv", full.names = T)
comm_list75 <- list.files(path = here::here("Data_BBS/Generated DFs/OccMod_CommMats"), pattern = "*DF75.csv", full.names = T)
comm_list95 <- list.files(path = here::here("Data_BBS/Generated DFs/OccMod_CommMats"), pattern = "*DF95.csv", full.names = T)

#List RDS files 
comm_list <- list.files(path = here::here("R Workspace/Spp_Arrays"), pattern = "*.RDS", full.names =T)

#Read files in for each type
prime.files <- lapply(comm_list, function(x) readRDS(x))
prime.files <- lapply(comm_list, function(x) read.csv(x, stringsAsFactors = F))
prime.files5 <- lapply(comm_list5, function(x) read.csv(x, stringsAsFactors = F))
prime.files65 <- lapply(comm_list65, function(x) read.csv(x, stringsAsFactors = F))
prime.files75 <- lapply(comm_list75, function(x) read.csv(x, stringsAsFactors = F))
prime.files95 <- lapply(comm_list95, function(x) read.csv(x, stringsAsFactors = F))

#Turn all arrays into DF
#prime.df <- lapply(prime.files, function(x) array2df(x, levels = list(mcmc.iter = T, spp.id = T, site.id = T, year.id = T), label.x = "Occupancy"))

#Bind all of the DFs into 5 large DFs
prime.df <- do.call(rbind, prime.df)
prime.df5 <- do.call(rbind, prime.files5) 
prime.df65 <- do.call(rbind, prime.files65)
prime.df75 <- do.call(rbind, prime.files75)
prime.df95 <- do.call(rbind, prime.files95)

###############################################################################
###################################End Script##################################
###############################################################################