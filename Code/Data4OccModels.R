################################################
########Cleaning Data for Occupancy Models######
##########Author: Spencer R Keyser##############
################################################

#####Package Loading#####
library(pacman)
pacman::p_load(tidyverse, reshape2, ggplot2)

#####Packages Loaded######

#####Adding in UDFs for Occ######

expit <- function(p){
  return(exp(p) / (exp(p) + 1))
}

logit <- function(p, e = 0.0005){
  p = ifelse(p == 1, p- e, p)
  p = ifelse(p == 0, p + e, p)
  return(log(p / (1 - p)))
}

######Standardizing functions######
######Generic Function Gelman's Suggestion######

standardise <- function(xmat, stdev = 1, marg = c(1, 2, 3)){
  mean.xmat = mean(as.vector(xmat), na.rm = T)
  sd.xmat = sd(as.vector(xmat), na.rm = T)
  std.xmat = apply(xmat, marg, function(x){
    (x - mean.xmat) / (stdev * sd.xmat)})
  retun(std.xmat)
}

scale2 <- function(x, na.rm = F) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


#####Standardizing Vectors#####
standardise.vector <- function(x, stdev = 1){(x - mean(x, na.rm = T))/
    (stdev * sd(x, na.rm = T))}

#####Ending Standardization Functions#####

###Load in DFs###

final_sp_df <- read.csv(file = here::here("Data_BBS/Generated DFs/Final_Sp_Df_DetectionsOnly.csv"), header = T)
site.occ.df <- read.csv(file = here::here("Data_BBS/Generated DFs/Site.Occ.csv"), header = T)
spp.occ <- read.csv(file = here::here("Data_BBS/Generated DFs/Spp.Occ.csv"), header = T)
Years.Occ <- read.csv(file = here::here("Data_BBS/Generated DFs/Years.Occ.csv"), header = T)
bcr.occ <- read.csv(file = here::here("Data_BBS/Generated DFs/bcr.detection.raw.csv"), header = T)

#Standardize variables for analysis

##Long way to scale variables
#site.occ.df <- site.occ.df %>% mutate(TOD.scaled = scale(TOD, center = T, scale = T))
#site.occ.df$TOD.scaled <- as.numeric(site.occ.df$TOD.scaled)
#site.occ.df <- site.occ.df %>% mutate(Ordinal.scaled = scale(OrdinalDate, center = T, scale = T))
#site.occ.df$Ordinal.scaled <- as.numeric(site.occ.df$Ordinal.scaled)
#site.occ.df$ChangeD.scaled <- site.occ.df$ChangeD

site.occ.scaled <- site.occ.df %>% mutate_at(c("TOD", "OrdinalDate"), scale2)

site.occ.scaled$ChangeD.scaled[site.occ.scaled$ChangeD.scaled == 0] <- -1

spp.occ <- spp.occ %>% mutate(Mass.scaled = scale(BodyMass, center = T, scale = T))

#Finished Standardizing and Centering

####Create Matrices
site.occ.ma <- site.occ.scaled %>% select(site, site.id, Year,
                                     year.id, TOD, OrdinalDate,
                                     ChangeD.scaled) %>%
               arrange(site.id)

##Site x Year Matrix with scaled.tod values and NAs for missing

#Time of day matrix
TOD.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum, 
                value.var = "TOD", fill = NA_real_, drop = F)
TOD.ma <- TOD.ma[, -1]
TOD.ma <- as.matrix(TOD.ma)
colnames(TOD.ma) <- NULL

#Ordinal Date Matrix
Ord.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum,
                value.var = "OrdinalDate", fill = NA_real_, drop = F)
Ord.ma <- Ord.ma[, -1]
Ord.ma <- as.matrix(Ord.ma)
colnames(Ord.ma) <- NULL

#1st year observer matrix (1 new obs, -1 same observer)
Obs.ma <- dcast(site.occ.ma, site.id ~ year.id, fun.aggregate = sum,
                value.var = "ChangeD.scaled", fill = NA_real_, drop = F)
Obs.ma <- Obs.ma[, -1]
Obs.ma <- as.matrix(Obs.ma)
colnames(Obs.ma) <- NULL


#####Import species specific data#####
#Creation of the complete ydf matrix
#alternative approach ## 
#head( spp.df )
head(  final_sp_df )
head( site.occ.df )
glimpse( site.occ.df )
head( spp.occ)

#check for duplicates
site.occ.df[ duplicated( site.occ.df ), ] #none!

#create reduced dataframe from final_sp_df
spp.df <- final_sp_df %>% dplyr::select( sci_name, site, Year, Detected )
#view
head( spp.df ); dim( spp.df )
#remove zero detections:
spp.df <- spp.df %>% dplyr::filter( Detected == 1 )
#view
head( spp.df ); dim( spp.df )
glimpse( spp.df )
#append species data
#first relabel first column
colnames( spp.df )[ 1 ] <- "Species"
colnames( spp.occ )[ 1 ] <- "Species"
spp.df <- left_join( spp.df, spp.occ, by = "Species" )
#check
#view
head( spp.df ); dim( spp.df ) #it didn't add columns..hooray!
#check for duplicates:
spp.df[ which( (spp.df$spp.id == 1) & ( spp.df$site == '2141_4') ),  ]
#append siteXyear info #ensure we turn all.y=F so that sitesXyear unsurveyed are not
# added
spp.df <- left_join( spp.df, site.occ.df, by = c( "site", "Year" ), all.x=T, all.y=F )
#view
tail( spp.df ); dim( spp.df ) #row numbers stayed the same
#did it add rows
sum( spp.df$Detected ) 
#check for duplicates
spp.df[ which( (spp.df$spp.id == 1) & ( spp.df$site.id == 12) ),  ]


#remove duplicates for species for now #you need to not do this with the new
#final df. #####
spp.df <- spp.df %>% group_by( spp.id, site.id, year.id ) %>%
  slice( 1 ) #only keeps one record for each speciesXsiteXyear
tail( spp.df ); dim( spp.df )
#removed ~200,000 records

# setdimensions ###
#total number of species
S <- max(spp.df$spp.id)
#total number of segments
J <- max(spp.df$site.id)
#total number of routes
M <- max(spp.df$rteno.id)
#total number of sampling years
K <- max(spp.df$year.id)


### working out missing sampling years for a given segment:
#we will work from site.occ.df which has all details of which segments were surveyed 
# which year:
#view
tail( site.occ.df ); dim( site.occ.df )
site.occ.df[ which( site.occ.df$site.id == 1), ]
spp.df[ which( spp.df$site.id == 1), ]
glimpse( site.occ.df )
#use it to create complete dataframe
JKdf <- site.occ.df %>% dplyr::select( site.id, year.id )
#check for duplicates
JKdf[ duplicated( JKdf ), ] #some duplicates present!!!!
#remove
#JKdf <- JKdf[ !duplicated( JKdf ), ]
#convert year code to factor
JKdf$year.id <- as.factor( JKdf$year.id )
#check
levels( JKdf$year.id )
#convert site code to factor
JKdf$site.id <- as.factor( JKdf$site.id )
#check
levels( JKdf$site.id )
#add surveyed column
JKdf$surveyed <- 1 
#view
head( JKdf ); dim( JKdf )
#add missing combinations
JKdf <- JKdf %>% tidyr::complete( site.id, year.id )
#view
head( JKdf ); dim( JKdf )
#dimensions should equal: 
J*K
#turn to wide format
JKmat <- tidyr::spread( JKdf, key = year.id, value = surveyed )
#it has to have J rows and K columns
head( JKmat); dim( JKmat )
#convert to matrix:
JKmat <- as.matrix( JKmat[ ,2:dim( JKmat )[2] ] )
# #remove column names from matrix
colnames( JKmat ) <- NULL
#view
head( JKmat); dim( JKmat )
#define how many segments were surveyed each year:
surveyedJ <- colSums( JKmat, na.rm = TRUE )
####### now we create the ydf ##### 
#create observations dataframe
ydf <- array( 0, dim = c(S, J, K) )
ydf[ 1, , ] * JKmat
# assigned 1 when species was detected on given year and route
for( i in 1:dim( spp.df )[1] ){
  ydf[ as.numeric( spp.df[i,'spp.id'] ), as.numeric(spp.df[i,'site.id']),
       as.numeric( spp.df[i,'year.id'] ) ] <- 1 #as.numeric( spp.df[i, 'Detected'] )
}

# #check that it worked for species 1:
table( ydf[4,,] )
table( spp.df$Detected[ which( spp.df$spp.id == 4) ] )  #77 
spp.df[ which( spp.df$spp.id == 4), ]
ydf[1,8,32]

#add missing values 
for( i in 1:S ){
  ydf[ i, , ] <- ydf[ i, , ] * JKmat
}

ydf[1, 17, ]

#Code below for augmented dataset
#Taken from Zipkin et al. 2010 
#(github.com/zipkinlab/Community_model_examples-covariate_model/blob/master/covariate%20model%20code.r)
#Create all zero encounter histories to add to the detection array X 
#as part of the data augmentation to account for additional 
#species (beyond the n observed species). 

#nzeroes is the number of all zero encounter histories to be added
nzeroes = 50
#X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled  
X.zero = matrix(0, nrow=70, ncol=4)
X.zero[1:36,4] = NA;  X.zero[38:56,4] = NA;  
X.zero[59:61,4] = NA;  X.zero[66:70,4] = NA; 

ydf.aug <- array(0, dim=c(dim(ydf.aug)[1],dim(ydf.aug)[2],dim(ydf.aug)[3]+nzeroes))
ydf.aug[,,(dim(ydf.aug)[3]+1):dim(ydf.aug)[3]] = rep(ydf.zero, nzeroes)
dimnames(X)=NULL
Xaug[,,1:dim(X)[3]] <-  X

# #total number of species
# S <- length(unique(spp.df$spp.id))
# #total number of segments
# J <- length(unique(spp.df$site.id))
# #total number of routes
# M <- length(unique(spp.df$rt.id))
# #total number of sampling years
# K <- length(unique(spp.df$year.id))

# #observations dataframe
# ydf <- array(0, dim = c(S, J, K))
# 
# # assigned 1 when species was detected on given year and route
# for( i in 1:dim( spp.df )[1] ){
#   ydf[ as.numeric( spp.df[i,'spp.id'] ), as.numeric(spp.df[i,'site.id']), 
#        as.numeric( spp.df[i,'year.id'] ) ] <- 1
# }
# #view
# head(spp.df)
# #check that it worked for species 1:
# sum( colSums( ydf[1,,] ) ) #72
# nrow( spp.df[ which( spp.df$spp.id == 1), ] ) #77 
# #what is missing?
# check <- spp.df[ which( spp.df$spp.id == 1), ]
# #check
# table( check$year.id )
# colSums( ydf[1,,] )
# #year.id =32 species 1 has 11 records in spp.df & only 6 in ydf. How?
# #tried reruning and rechecking
# #check year.id = 32
# ydf[1,,32]
# check[ which( check$year.id == 32),]

############################################

#save.image(file = here::here("R Workspace/Data4OccModels_5_24.RData"))