#Hi Git

#Cleaning command 
rm(list = objects(all.names = TRUE))
install.packages("pacman")
pacman::p_load("reshape2", "lefse", "genpathmox", "tidyverse", "dplyr", "compare")

setwd("C:/Users/Spencer/Box Sync/Keyser Thesis/Data/Functional Traits")


#Trait files 
elton <- read.csv(file = "Functional_Traits_ESA_Jetz_updated_csv.csv", header = TRUE)
#lcta <- read.csv(file = "Army_Corps_Bird_Data.csv", header = TRUE)
#lislevand <- read.csv(file = "Bird_Data_Mating.csv", header = TRUE)

##Standardize the names for merging 
colnames(elton)[colnames(elton) == "Scientific"] <- "sci_name"
colnames(elton)[colnames(elton) == "English"] <- "COM_NAME"
#colnames(lcta)[colnames(lcta) == "SCIENTIFIC_NAME"] <- "sci_name"
#colnames(lislevand)[colnames(lislevand) == "Species_name"] <- "sci_name"

##Merge trait data
#trait_data <- merge(elton, lcta, by = "sci_name")
#trait_data <- merge(trait_data, lislevand, by = "sci_name")

##Upload community data
cbc <- read.csv(file = "CBC_Species_Sites_Total_Refined_NoCarib_Cleaned_CSV.csv", header =TRUE)
length(unique(cbc$SCI_NAME))
cbc <- cbc[which(cbc$SCI_NAME != "Porphyrio porphyrio"),]
cbc$unique_id <- paste0(cbc$Abbrev, "_", cbc$Count_yr)

##Removing subspecies issues 
cbc_sep <- tidyr::separate(cbc, SCI_NAME, c("genus", "species", "subspecies"))
cbc_sep$subspecies <- NA
cbc_sep$species[sub(" +$", "", cbc_sep)]
cbc_sep <- cbc_sep[ ,-12]
cbc_sep <- tidyr::unite(cbc_sep, sp, genus, species, sep = " ", remove = FALSE)
colnames(cbc_sep)[colnames(cbc_sep) == "sp"] <- "sci_name"

#Sum species observations in 5 year bins
#set year bins
years <- data.frame(yr = 80:117, yr_bin = NA)
years$yr_bin[years$yr <= 84] <- 1
years$yr_bin[years$yr >= 85 & years$yr <= 89] <- 2
years$yr_bin[years$yr >= 90 & years$yr <= 94] <- 3
years$yr_bin[years$yr >= 95 & years$yr <= 99] <- 4
years$yr_bin[years$yr >= 100 & years$yr <= 104] <- 5
years$yr_bin[years$yr >= 105 & years$yr <= 109] <- 6
years$yr_bin[years$yr >= 110 & years$yr <= 114] <- 7
years$yr_bin[years$yr >= 115 & years$yr <= 120] <- 8

cbc_binned <- merge(cbc_sep, years, by.x = "Count_yr", by.y = "yr")
cbc_binned <- cbc_binned[,c("Abbrev", "sci_name", "how_many", "yr_bin", "unique_id")]
cbc_binned_agg <- aggregate(how_many ~ Abbrev +yr_bin +sci_name, cbc_binned, mean)
cbc_binned_agg$unique_id <- paste0(cbc_binned_agg$Abbrev, "_", cbc_binned_agg$yr_bin)

#Casting dataset 
cbc_simple <- cbc_binned[, c("unique_id", "sci_name", "how_many")]
cbc_cast <- dcast(data = cbc_simple, formula = unique_id ~ sci_name, fun.aggregate = sum)
rownames(cbc_cast) <- cbc_cast[, 1]
cbc_cast <- cbc_cast[, -1]


##Extract birds from the CBC for species pool 
##Mismatch in the length of the vectors for species recorded and the traits database
##likely from subspecies and "group" reporting, should be fine to leave out for dendrogram 
cbc_trait <- cbc_sep[, c("sci_name", "COM_NAME", "Abbrev")]
colnames(cbc_trait)[colnames(cbc_trait) == "SCI_NAME"] <- "sci_name"

cbc_unique <- unique(cbc_trait[,c("sci_name", "COM_NAME")])
cbc_unique <- as.data.frame(cbc_unique)
colnames(cbc_unique)[colnames(cbc_unique) == "cbc_unique"] <- "sci_name"

#test <- gsub(" ", "", cbc_unique$sci_name, fixed = TRUE)
#test <- as.data.frame(test)
#colnames(test)[colnames(test) == "test"] <- "sci_name"
#length(unique(test$sci_name))
#test.e <- gsub(" ","", elton$sci_name, fixed = TRUE)
#test.e <- as.data.frame(test.e)
#colnames(test.e)[colnames(test.e) == "test.e"] <- "sci_name"
#length(unique)


cbc_unique <- cbc_unique[which(cbc_unique$sci_name != "Porphyrio porphyrio"),]
elton_merge <- merge(cbc_unique, elton, by = "sci_name")
testing <- anti_join(cbc_unique, elton_merge, by = "sci_name")
#Removr Western swamphen, it does not have trait data and is a rare species in the US 
#At this point the elton merge df is set to go 


#Need t isolate the observations that are different between both of these data frames
length(unique(cbc$SCI_NAME))
length(unique(elton$sci_name))
length(unique(elton_merge$sci_name))
length(unique(cbc_unique$sci_name))

elton_diet <- elton_merge[, c(1, 10:20)]
elton_diet <- elton_diet[, -2]
elton_diet <- elton_diet %>% distinct(sci_name, .keep_all = TRUE)
elton_diet$Inverts <- elton_diet$Inverts * .01
elton_diet$Mammal.Bird <- elton_diet$Mammal.Bird * .01
elton_diet$Herps <- elton_diet$Herps * .01
elton_diet$Fish <- elton_diet$Fish * .01
elton_diet$Unknown_Verts <- elton_diet$Unknown_Verts * .01
elton_diet$Scavengers <- elton_diet$Scavengers * .01
elton_diet$Fruit <- elton_diet$Fruit * .01
elton_diet$Nect <- elton_diet$Nect * .01
elton_diet$Seed <- elton_diet$Seed * .01
elton_diet$Other_plant <- elton_diet$Other_plant * .01

#********#Come back to this##******************
#Need to fix the duplication of species 
#dimnames(elton_diet) <- list("sci_name"=as.character(elton_diet[,1]),"traits"=dimnames(elton_diet)[[3:11]])
rownames(elton_diet) <- elton_diet[,1]
elton_diet <- elton_diet[, -1]

##Attempting to calculate a distance matrix 
# distance <- dist(elton_diet)
# mat <- as.matrix(distance)
# na <- is.na(mat)
# sum(na)
#elton_matrix <- as.matrix(elton_diet)
#diet_dist <- dist(elton_matrix, method = "euclidean", diag = FALSE)
#diet_mat <- as.matrix(diet_dist)
#rownames(diet_mat) <- diet_mat$species

#Create a spp by site matrix
#Keep only the columns we need
#cbc_simple <- cbc_sep[,c("unique", "sci_name", "how_many")]
#cbc_plot <- cbc_sep[, c("unique_id", "subnational_code")]
#cbc_plot <- as.data.frame(cbc_plot)
#Cast Data
#cbc_casted <- dcast(data = cbc_simple, formula = unique ~ sci_name, fun.aggregate = sum)
#rownames(cbc_casted) <- cbc_casted$unique
#cbc_casted <- cbc_casted[, -1]
#as.matrix(cbc_cast)

##Exclude nocturnal species
#remove caprimulgiformes and strigiformes?
#It may be best not to remove these species unless they are going to be removed all together 
# test <- elton_merge$Order[which(elton_merge$Order != c("Strigiformes", "Caprimulgiformes"))]
# test <- as.data.frame(test)

### Nate's Code

##Creating a distance matrix
#Fmntd <- function(dist.mat, my.sample){


#   Fmntd.sub = function(x){ 

## Get the names of the species present in a 
## community.
#    com.names = names(x[x > 0])

## Make the community phylogenetic distance 
## matrix by extracting those rows and columns 
## that have species present in our community.
#   my.com.dist = dist.mat[com.names, com.names]

## Set all diagonal values to NA so that the 
## zeros for conspecific comparisons do not 
## interfere with our calculation of nearest 
## neighbors.
#  diag(my.com.dist) = NA

## Use apply() to calculate the minimum value in 
## each row of the community phylogenetic 
## distance matrix and take a mean of those 
## values.
# mean(apply(my.com.dist, MARGIN = 1, min, na.rm=T), na.rm=T)

#}

#apply(my.sample, MARGIN = 1, Fmntd.sub)

#}

##Traits
##Code is generating issues for subscript out of bounds issue
##This might be because of a mismatch in the community data with the trai data? 
##Or this code is set for the data it was designed under 
my.sample <- cbc_cast
#my.dist.mat <- distance

mat = as.matrix(dist(elton_diet))

mat <- mat[rownames(mat) %in% colnames(cbc_cast), colnames(mat)%in%colnames(cbc_cast)]



#Loop throught each site
site.list <- unique(cbc_binned[,c("Abbrev", "unique_id")])
site.list$Abbrev <- as.character(site.list$Abbrev)
site.list$unique_id[sub(" +$", "", site.list)]
sites <- as.character(unique(site.list$Abbrev))

#create data.frame to hold beta output
betas.out <- data.frame(second_bin = site.list$unique_id, beta.from.baseline = NA)

# Loop through sites and extract the Functional Beta diversity from each year bin site 
# to compare to the baseline year bin

for (u in 1:length(sites)){
  site.temp <- sites[u]
  surveys.temp <- site.list$unique_id[site.list$Abbrev == site.temp]
  
  
  cbc_temp <- cbc_cast[row.names(cbc_cast) %in% surveys.temp,]
  
  
  betas.temp <- new.comdistnn(my.sample = cbc_temp, my.dist.mat = mat)
  
  for (y in 1:nrow(betas.temp)){
    temp.unique.id <- row.names(betas.temp)[y]
    betas.out[betas.out$second_bin == temp.unique.id,"beta.from.baseline"] <- betas.temp[y,1]}
  
}

betas.out$second_bin <- as.character(betas.out$second_bin)  
betas.out <- separate(betas.out, second_bin, c("Abbrev", "Count_yr"))
betas.out <- unite(betas.out, unique_id, Abbrev, Count_yr, sep = "_")
betas.out$unique2 <- betas.out$unique_id
betas.out <- separate(betas.out, unique2, c("Abbrev", "Count_yr"))
betas.reserve <- betas.out
#lefse example
#

new.comdistnn = function(my.sample, my.dist.mat){
  
  get.presents = function(x){
    names(x)[x>0]
  }
  
  list.of.names = apply(my.sample, 1, get.presents)
  
  Dnn.apply.function = function(x){
    tmp.function = function(z){
      
      mean(c(apply(as.matrix(my.dist.mat[x,z]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[x,z]), MARGIN=2, min, na.rm=T)))
      
    }
    
    lapply(list.of.names, FUN=tmp.function)
    
    # for (t in 1:length(list.of.names)){
    #   tmp.row <- list.of.names[[t]]
    #   tmp.function(tmp.row)
    # }
  }
  
  dnn.output = lapply(list.of.names,Dnn.apply.function)
  
  
  outt = do.call(rbind,lapply(dnn.output,unlist))
  
  outt
  
}


n.sites <- length(unique(betas.out$Abbrev))
site.list <- as.character(unique(betas.out$Abbrev))

# #Put betas back in site daTA MERGE
# site_data_merge <- merge(site_data_merge, beta.matrix, by = "unique_id")

#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, beta.slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  betas.out <- betas.out[betas.out$Abbrev == site.temp,]
  betas.out <- betas.out[complete.cases(betas.out),]
  if (nrow(betas.out) > 1){
    lm.temp <- lm(betas.out$beta.from.baseline ~ betas.out$Count_yr)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
  }
}

