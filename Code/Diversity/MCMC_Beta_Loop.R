#Bring in all of the scripts (work this up with 5)
df.list <- list.files(here::here("Data_BBS/Generated DFs/MCMC DFs"), pattern = "*.csv", full.names = T)

#Load in the necessary DFs to subset with
JKsurv <- readRDS(here::here("Data_BBS/Generated DFs/JKmat.RDS"))
jdf <- read.csv(here::here("Data_BBS/Generated DFs/jdf.csv"))
yrs <- read.csv(here::here("Data_BBs/Generated DFs/Years.Occ.csv"))

#Yr_bins
Bins <- data.frame(matrix(ncol = 2, nrow = 40))
l <- c("Year", "Yr_bin")
colnames(Bins) <- l
Bins$Year <- 1980:2019
Bins$Yr_bin <- rep(1:8, times = 1, each = 5)

#Melt the JKsurv DFs
#Keeps the dfs in their format without needing to convert back to matrices
JKdf <- arrayhelpers::array2df(JKsurv, levels = list(site.id = T, year.id = T), label.x = "Surveyed")
JKdf <- lapply(JKdf, as.numeric)
JKdf <- as.data.frame(JKdf)

#Make JKdf usable
JKdf <- JKdf %>% dplyr::left_join(jdf, by = "site.id") %>% dplyr::left_join(yrs, by = "year.id")

#Remove all non-surveyed years
JKclean <- JKdf[complete.cases(JKdf),]
JKclean$site_yr <- paste0(JKclean$site, "_", JKclean$Year)

#Pull out the unique survey_yr combos
surveys <- unique(JKclean$site_yr)

#Create a similar variable in the DFs
df.cleaner <- function(x){
x$site_yr <- paste0(x$site, "_", x$Year)

#Subset DF with the surveys 
x <- x[x$site_yr %in% surveys, ]
x <- x %>% left_join(Bins, by = "Year") 
#df1.new$site_yr <- paste0(df1.new$site, "_", df1.new$Yr_bin)
}


#mean temporal beta diversity by site
mcmc <- paste0("Beta", "_", seq(100))
beta.matrix <- data.frame(matrix(ncol = 100, nrow = 4925))
colnames(beta.matrix) <- mcmc
rownames(beta.matrix) <- surveys


#Create matrix for storing Betas 
#beta.matrix <- data.frame(unique_id = bbs_simple$unique_id, beta = NA)
#beta.matrix$unique_id <- as.character(beta.matrix$unique_id)

#Reove duplites
#beta.matrix <- beta.matrix[!duplicated(beta.matrix),]

n.sites <- length(unique(JKclean$site))
site.list <- as.character(unique(JKclean$site))

#Loop through the sites so that we calculate beta based on year bins 
for (q in 1:length(df.list)){
  df.temp <- read.csv(file = df.list[q]) 
  df.temp <- df.cleaner(df.temp)
  bbs_cast <- dcast(df.temp, site_yr ~ Species, value.var = "Occupancy", fun.aggregate = sum)
  rownames(bbs_cast) <- bbs_cast[, 1]
  bbs_cast <- bbs_cast[, -1]
  bbs_cast[bbs_cast >= 1] <- 1
  
  for (n in 1:n.sites){
    site.temp <- site.list[n]
    #Find all the years for which that site was surveys
    yrs.temp <- unique(JKclean[which(JKclean$site == site.temp),"Year"])
    yrs.temp <- yrs.temp[yrs.temp >= 1]
    #only consider sites that were observed in more than 5 years 
    if (length(yrs.temp) > 1){ #select all sites that have mre than 1 year bin 
      #Create a baseline year from years temp
      #Sorting
      yrs.temp <- yrs.temp[order(yrs.temp)]
      #find the first 5 years 
      first.yrs <- yrs.temp[yrs.temp <= (yrs.temp[1])] #pulls first year bin
      print(paste0("Baseline Year for ", site.temp, " is ", first.yrs))
      site.yrs.first <- paste0(site.temp, "_", first.yrs)
      #create average community of first 5 years 
      bbs_cast_first <- bbs_cast[rownames(bbs_cast) %in% site.yrs.first,]
      first.comm.temp <- bbs_cast_first 
      #get remaining years 
      other.yrs.temp <- yrs.temp[yrs.temp > (yrs.temp[1])]
      #Now loop the remaining years 
      for (y in 1:length(other.yrs.temp)){
        other.yr.temp <- other.yrs.temp[y] #pull out another yr
        other.site.temp <- paste0(site.temp, "_", other.yr.temp) #make the unique id
        bbs_cast_other <- bbs_cast[rownames(bbs_cast) == other.site.temp,] #pull the comm data
        bbs_cast_other <- rbind(bbs_cast_other, first.comm.temp) #put the current yr and baseline
        beta.temp <- vegdist(sqrt(sqrt(bbs_cast_other)), method = "jaccard") #calc jaccard disim
        beta.matrix[rownames(beta.matrix) == other.site.temp, q] <- beta.temp #store it in df
      } #y
    } #if statement
  } #n
} #q

#Calculate a Mean Beta for all of the iterations
beta.matrix2 <- beta.matrix[complete.cases(beta.matrix),]

#Create a DF to store summary stats
beta.means <- data.frame(Sites = beta.matrix2$Sites, mean.beta = NA, beta.sd = NA)

#Calculate means and sds for beta diversity
rownames(beta.matrix2) <- beta.matrix2$Sites #makes the sites the rows
beta.matrix2 <- beta.matrix2[, -1] #Remove sites foe calculations

beta.means$mean.beta <- rowMeans(beta.matrix2) #calculate means across rows
beta.means$beta.sd <- apply(beta.matrix2, 1, sd) #calculate SD
#write.csv(beta.means, file = here::here("Data_BBS/Generated DFs/MCMC DFs/beta_means.csv"))





