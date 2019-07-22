#########################################################################################
#######################Script for turning MSOM Output back into##########################
#######################Community Data for Biodiversity Analyses##########################
#############################Script includes 100 iterations############################## 
#################################for mean beta diversity#################################
#########################################################################################
##########################Script Created By: Spencer R Keyser############################
##################################### 7/22/2019 #########################################
#########################################################################################

#Clean the environment
rm(list = ls())

#Load in Packages

#library("pacman")
pacman::p_load("tidyverse", "plyr", "reshape2", "here", "arrayhelpers")

#Finish Package Loading 

#Load in data for sites 
jdf <- read.csv(file = here::here("Data_BBS/Generated DFs/jdf.csv"))
jdf <- jdf[, -1]
yrs <- read.csv(file = here::here("Data_BBs/Generated DFs/Years.Occ.csv"))
JKmat <- readRDS(file = here::here("Data_BBS/Generated DFs/JKmat.RDS"))
JKsurv <- readRDS(file = here::here("Data_BBS/Generated DFs/JKsurv.RDS"))

#Load in spp.occ DFs
spp_list <- list.files(path = here::here("Data_BBS/Generated DFs/SppOcc"), pattern = "*.csv", full.names = T)
spp_csv <- lapply(spp_list, function(x) read.csv(x, stringsAsFactors = F))

#Rename all of the DFs in the list
names(spp_csv) <- c("Kingfishers", "Cuckoos", "Falcons", "Flycatchers", "Fowl", "Swallows",
                    "Hummers", "Blackbirds", "Corvids", "Warblers", "Finches", "Rails",
                    "Wrens", "Raptors", "Tits", "Shorebirds", "Sparrows", "Owls",
                    "Mimids", "Waders")

#Load them into the global environment
list2env(spp_csv, envir = .GlobalEnv)

rm(spp_csv)

#UDF to change from factors for joins
no.factors <- function(x){
  if(is.factor(x) == T){
    as.integer(x)
   } else {
     as.numeric(x)
   }
 }

####Read in all of the new community matrices####
#List RDS files 
comm_list <- list.files(path = here::here("R Workspace/Spp_Arrays"), pattern = "*.RDS", full.names =T)


#Read files in for each type
prime.files <- lapply(comm_list, function(x) readRDS(x))

names(prime.files) <- c("BlackbirdsDF", "CorvidsDF", "CuckoosDF", "FalconsDF", "FinchesDF",
                        "FlycatchersDF", "FowlDF", "HummersDF", "KingfishersDF", "MimidsDF",
                        "OwlsDF", "RailsDF", "RaptorsDF", "ShorebirdsDF", "SparrowsDF", "SwallowsDF", 
                        "TitsDF", "WadersDF", "WarblersDF", "WrensDF")

z2df <- lapply(prime.files, function(x) array2df(x, levels = list(mcmc = T, spp.id = T,
                                                                  site.id = T, year.id = T),
                                                 label.x = "Occupancy"))

rm(prime.files)

list2env(z2df, envir = .GlobalEnv)

BlackbirdsDF <- as.data.frame(lapply(BlackbirdsDF, no.factors))
BlackbirdsDF <- BlackbirdsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Blackbirds, by = "spp.id")

CorvidsDF <- as.data.frame(lapply(CorvidsDF, no.factors))
CorvidsDF <- CorvidsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Corvids, by = "spp.id")

CuckoosDF <- as.data.frame(lapply(CuckoosDF, no.factors))
CuckoosDF <- CuckoosDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Cuckoos, by = "spp.id")

FalconsDF <- as.data.frame(lapply(FalconsDF, no.factors))
FalconsDF <- FalconsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Falcons, by = "spp.id")

FinchesDF <- as.data.frame(lapply(FinchesDF, no.factors))
FinchesDF <- FinchesDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Finches, by = "spp.id")

FlycatchersDF <- as.data.frame(lapply(FlycatchersDF, no.factors))
FlycatchersDF <- FlycatchersDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Flycatchers, by = "spp.id")

FowlDF <- as.data.frame(lapply(FowlDF, no.factors))
FowlDF <- FowlDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Fowl, by = "spp.id")

HummersDF <- as.data.frame(lapply(HummersDF, no.factors))
HummersDF <- HummersDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Hummers, by = "spp.id")

KingfishersDF <- as.data.frame(lapply(KingfishersDF, no.factors))
KingfishersDF <- KingfishersDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Kingfishers, by = "spp.id")

MimidsDF <- as.data.frame(lapply(MimidsDF, no.factors))
MimidsDF <- MimidsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Mimids, by = "spp.id")

OwlsDF <- as.data.frame(lapply(OwlsDF, no.factors))
OwlsDF <- OwlsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Owls, by = "spp.id")

RailsDF <- as.data.frame(lapply(RailsDF, no.factors))
RailsDF <- RailsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Rails, by = "spp.id")

RaptorsDF <- as.data.frame(lapply(RaptorsDF, no.factors))
RaptorsDF <- RaptorsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Raptors, by = "spp.id")

ShorebirdsDF <- as.data.frame(lapply(ShorebirdsDF, no.factors))
ShorebirdsDF <- ShorebirdsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Shorebirds, by = "spp.id")

SparrowsDF <- as.data.frame(lapply(SparrowsDF, no.factors))
SparrowsDF <- SparrowsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Sparrows, by = "spp.id")

SwallowsDF <- as.data.frame(lapply(SwallowsDF, no.factors))
SwallowsDF <- SwallowsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Swallows, by = "spp.id")

TitsDF <- as.data.frame(lapply(TitsDF, no.factors))
TitsDF <- TitsDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Tits, by = "spp.id")

WadersDF <- as.data.frame(lapply(WadersDF, no.factors))
WadersDF <- WadersDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Waders, by = "spp.id")

WarblersDF <- as.data.frame(lapply(WarblersDF, no.factors))
WarblersDF <- WarblersDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Warblers, by = "spp.id")

WrensDF <- as.data.frame(lapply(WrensDF, no.factors))
WrensDF <- WrensDF %>% left_join( yrs, by = "year.id") %>%
  left_join(jdf, by = "site.id") %>% 
  left_join(Wrens, by = "spp.id")


#Rbind dataframe back 
master <- rbind(BlackbirdsDF, CorvidsDF, CuckoosDF, FalconsDF, FinchesDF, FlycatchersDF,
                FowlDF, HummersDF, KingfishersDF, MimidsDF, OwlsDF, RailsDF, RaptorsDF,
                ShorebirdsDF, SparrowsDF, SwallowsDF, TitsDF, WadersDF, WarblersDF,
                WrensDF)

for(i in unique(df$y)) {
  nam <- paste("df", i, sep = ".")
  assign(nam, df[df$y==i,])
}

group.name <- as.character(unique(spp.occ$Phylo.V1))
group.name <- gsub("/", "_", group.name)
group.name1 <- paste0(group.name, "_", "CommunityDF.csv")


write.csv(df.prime, file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", group.name1)))


#image.name <- paste0(group.name, "", ".RData")
#save.image(file = here::here(paste0("Data_BBS/Generated DFs/OccMod_CommMats", "/", image.name)))

####Read in all of the new community matrices####
#List RDS files 
comm_list <- list.files(path = here::here("R Workspace/Spp_Arrays"), pattern = "*.RDS", full.names =T)

#Read files in for each type
prime.files <- lapply(comm_list, function(x) readRDS(x))

#Turn all arrays into DF
#prime.df <- lapply(prime.files, function(x) array2df(x, levels = list(mcmc.iter = T, spp.id = T, site.id = T, year.id = T), label.x = "Occupancy"))

#Bind all of the DFs into 5 large DFs
prime.df <- do.call(rbind, prime.df)


