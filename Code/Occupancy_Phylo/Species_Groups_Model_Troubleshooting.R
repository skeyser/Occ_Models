library("pacman")
pacman::p_load("tidyverse", "reshape2", "data.table")

final_sp_df <- read.csv(here::here("Data_BBS/Generated DFs/Final_Sp_Df_DetectionsOnly.csv"))
final_sp_df$Detected <- 1

sp.obs <- aggregate(Detected ~ spp.id, data = final_sp_df, sum)
sp.phylo <- aggregate(Detected ~ Phylo.V1.code, data = final_sp_df, sum)

sp.obs <- sp.obs %>% arrange(desc(Detected))
sp.phylo <- sp.phylo %>% arrange(desc(Detected))

hist(sp.obs$Detected, breaks = 8, col = "gray")

final_sp_df <- merge(final_sp_df, sp.obs, by = "spp.id")

#Trying to figure out the best grouping schemes
rare <- final_sp_df[final_sp_df$Detected.y < 50,]
uncommon <- final_sp_df[final_sp_df$Detected.y >= 50 & final_sp_df$Detected.y < 200,]
common <- final_sp_df[final_sp_df$Detected.y >= 200 & final_sp_df$Detected.y < 1000,]
abundant <- final_sp_df[final_sp_df$Detected.y >= 1000,] 

#Put all different DFs into a list
sp.groups <- list(rare, uncommon, common, abundant)

#Fast function to calculate the # of species in each group
sp.list <- function (x, y){
  list.temp <- x[[y]]
  species <- length(unique(list.temp$spp.id))
}

species.levels <- data.frame(Group = c("Rare", "Uncommon", "Common", "Abundant"),
                             SpeciesCount = NA)

#Loop through each DF and pull out each species count
for (i in 1:4){
  sp[i] <- sp.list(sp.groups, i)
  species.levels[i, 2] <- sp[i]
}#i

#See some functional traits 
species.info <- final_sp_df %>% select(c(sci_name, spp.id, BodyMass, Nocturnal,
                                         ForStrat_watbelowsurf, 
                                         ForStrat_wataroundsurf,
                                         ForStrat_ground, ForStrat_understory,
                                         ForStrat_midhigh, ForStrat_canopy,
                                         ForStrat_aerial, Detected.y))
species.info <- species.info[!duplicated(species.info),]

plot(species.info$BodyMass, species.info$Detected.y)

species.info$RA <- 1
species.info$RA[species.info$Detected.y < 50] <- 1
species.info$RA[species.info$Detected.y >= 50 & species.info$Detected.y < 200] <- 2
species.info$RA[species.info$Detected.y >= 200 & species.info$Detected.y < 1000] <- 3
species.info$RA[species.info$Detected.y >= 1000] <- 4

sp.bcr <- setDT(final_sp_df)[, .(count = uniqueN(sci_name)), by = BCR]
site.bcr <- setDT(final_sp_df)[, .(count = uniqueN(site)), by = BCR]
rt.bcr <- setDT(final_sp_df)[, .(count = uniqueN(rteno.x)), by = BCR]
yr.bcr <- setDT(final_sp_df)[, .(count = uniqueN(Year)), by = BCR]


#save.image(here::here("R Workspace/SpeciesGroupsModTroubleShooting.RData"))

