#Cleaning command 
rm(list = objects(all.names = TRUE))
install.packages("pacman")
pacman::p_load("reshape2", "cluster","lefse", "tidyverse", "cowplot", "here")

#Trait files 
elton <- read.csv(file = here("Functional_Traits/Functional_Traits_ESA_Jetz_updated_csv.csv"), header = TRUE) 

#lcta <- read.csv(file = "Army_Corps_Bird_Data.csv", header = TRUE)
#lislevand <- read.csv(file = "Bird_Data_Mating.csv", header = TRUE)

##Standardize the names for merging 
colnames(elton)[colnames(elton) == "Scientific"] <- "sci_name"
colnames(elton)[colnames(elton) == "English"] <- "COM_NAME"
#colnames(lcta)[colnames(lcta) == "SCIENTIFIC_NAME"] <- "sci_name"
#colnames(lislevand)[colnames(lislevand) == "Species_name"] <- "sci_name"


##Upload community data
temp = list.files(path = here("Data_BBS/States_GoM/States"), pattern="*.csv", full.names = T)
myfiles = lapply(temp, read.csv)
bbs.gom <- do.call(rbind, myfiles)

Species <- read.csv(here("Data_BBS/States_GoM/SpeciesList.csv"))


#Bin the Years 
bbs.gom <- bbs.gom[bbs.gom$Year >= 1980, ]

bbs.gom$yr.bin[bbs.gom$Year >= 1980 & bbs.gom$Year < 1985] <- 1 
bbs.gom$yr.bin[bbs.gom$Year >= 1985 & bbs.gom$Year < 1990] <- 2 
bbs.gom$yr.bin[bbs.gom$Year >= 1990 & bbs.gom$Year < 1995] <- 3 
bbs.gom$yr.bin[bbs.gom$Year >= 1995 & bbs.gom$Year < 2000] <- 4 
bbs.gom$yr.bin[bbs.gom$Year >= 2000 & bbs.gom$Year < 2005] <- 5 
bbs.gom$yr.bin[bbs.gom$Year >= 2005 & bbs.gom$Year < 2010] <- 6
bbs.gom$yr.bin[bbs.gom$Year >= 2010 & bbs.gom$Year < 2015] <- 7 
bbs.gom$yr.bin[bbs.gom$Year >= 2015 & bbs.gom$Year < 2017] <- 8 

colnames(bbs.gom)[colnames(bbs.gom) == "Aou"] <- "AOU"

#Shape and organize species data for merging
Species <- Species[, c(2, 3, 5, 6, 7, 8)]

#Change variables from factors to characters 
Species.character <- lapply(Species[2:6], as.character)

#Extract AOU
aou <- Species[,1]
Species.character <- as.data.frame(Species.character, stringsAsFactors = F)

#Bind the AOU with Names
Species <- cbind(aou, Species.character)
colnames(Species)[colnames(Species) == "aou"] <- "AOU"

#Merge Species identifiers with observation data
bbs.gom <- merge(bbs.gom, Species, by = "AOU")
bbs.gom$sci_name <- paste0(bbs.gom$Genus, " ", bbs.gom$Species)

# Creat Unique ID for aggregation, state + route + yr
bbs.gom$route.id <- paste0(bbs.gom$statenum, "_", bbs.gom$Route)
bbs.gom$uid <- paste0(bbs.gom$route.id, "_", bbs.gom$Year)

#Clean Scienific Names so that all unid. birds are dropped
#Drop species that are unidentified and that are genus level
bbs.clean <- bbs.gom
bbs.clean <- bbs.clean[!grepl("/" , bbs.clean$sci_name),]
bbs.clean <- bbs.clean[!grepl("sp.", bbs.clean$sci_name),]
bbs.clean <- bbs.clean[!grepl(" x ", bbs.clean$sci_name),]
bbs.clean <- bbs.clean[!bbs.clean$sci_name == "Porphyrio porphyrio",]
bbs.clean$sci_name[bbs.clean$sci_name == "Setophaga coronata audoboni"] <- "Setophaga coronata"
bbs.clean$sci_name[bbs.clean$sci_name == "Colaptes auratus cafer"] <- "Colaptes auratus"
bbs.clean$sci_name[bbs.clean$sci_name == "Colaptes auratus auratus"] <- "Colaptes auratus"
bbs.clean$sci_name[bbs.clean$sci_name == "Aphelocoma woodhouseii"] <- "Aphelocoma californica"
bbs.clean$sci_name[bbs.clean$sci_name == "Aphelocoma wollweberi"] <- "Aphelocoma ultramarina"
bbs.clean$sci_name[bbs.clean$sci_name == "Antrostomus arizonae"] <- "Antrostomus vociferus"
bbs.clean$sci_name[bbs.clean$sci_name == "Anas platyrhynchos diazi"] <- "Anas platyrhynchos"
bbs.clean$sci_name[bbs.clean$sci_name == "Ardea herodias occidentalis"] <- "Ardea herodias"

#Change some species names in Elton
elton.clean <- elton
elton.clean$sci_name <- as.character(elton.clean$sci_name)
elton.clean$sci_name[elton.clean$sci_name == "Pipilo fuscus"] <- "Melozone fusca"
elton.clean$sci_name[elton.clean$sci_name == "Grus canadensis"] <- "Antigone canadensis"
elton.clean$sci_name[elton.clean$sci_name == "Carduelis tristis"] <- "Spinus tristis"
elton.clean$sci_name[elton.clean$sci_name == "Carduelis psaltria"] <- "Spinus psaltria"

#Merge the trait data with the observation data 
elton_merge <- merge(bbs.clean, elton.clean, by = "sci_name")

#This will pull the species that were in obs data, but not found in trait data
#BBS Clean has 15 species that are not matched 
testing <- anti_join(bbs.clean, elton_merge, by = "sci_name")

#Average Total Amount for 5 year bins
bbs.agg <- aggregate(data = bbs.clean, SpeciesTotal ~ sci_name + route.id + yr.bin, FUN = "mean")

bbs.agg$rt.bin <- paste0(bbs.agg$route.id, "_", bbs.agg$yr.bin)

#Cast data to create a site x species matrix 
bbs.cast <- bbs.agg[, c("sci_name", "rt.bin", "SpeciesTotal")]
bbs.cast <- dcast(data = bbs.cast, formula = rt.bin ~ sci_name, fun.aggregate = sum)

#Make the site the rownames 
rownames(bbs.cast) <- bbs.cast[, 1]
#Remove rownames
bbs.cast <- bbs.cast[, -1]

#USe lapply function to perform operation on specific columns 
traits <- elton_merge[, c(1, 41, 56, 57, c(45:52))]

traits <- distinct(traits, sci_name, .keep_all = TRUE)

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

##Army Corps Traits




##Code is generating issues for subscript out of bounds issue
##This might be because of a mismatch in the community data with the trai data? 
##Or this code is set for the data it was designed under 
bbs.cast <- as.matrix(bbs.cast)
my.sample <- bbs.cast
rownames(traits) <- traits[, 1]
traits <- traits[, -1]

mat = as.matrix(daisy(traits), metric = "gower")
mat <- mat[rownames(mat) %in% colnames(bbs.cast), colnames(mat)%in%colnames(bbs.cast)]


my.dist.mat <- mat
#Loop throught each site
site.list.fun <- unique(bbs.agg[,c("route.id", "rt.bin")])

sites.fun <- unique(site.list.fun$route.id)

#create data.frame to hold beta output
betas.out <- data.frame(second_bin = site.list.fun$rt.bin, beta.from.baseline = NA)
betas.out.prime <- data.frame(second_bin = site.list.fun$rt.bin, beta.from.baseline = NA)
betas.out.prime.null <- data.frame(second_bin = site.list.fun$rt.bin, mean.from.null = NA, sd.from.null = NA)
# Loop through sites and extract the Functional Beta diversity from each year bin site 
# to compare to the baseline year bin

for (u in 1:length(sites.fun)){
  site.temp <- sites.fun[u]
  surveys.temp <- site.list.fun$rt.bin[site.list.fun$route.id == site.temp]
  
  bbs.temp <- bbs.cast[row.names(bbs.cast) %in% surveys.temp,]
  
  
  betas.temp <- new.comdistnn(my.sample = bbs.temp, my.dist.mat = mat)
  
  for (y in 1:nrow(betas.temp)){
    temp.unique.id <- row.names(betas.temp)[y]
    betas.out[betas.out$second_bin == temp.unique.id,"beta.from.baseline"] <- betas.temp[y,1]}
  
}

##With prime code 
for (k in 1:length(sites.fun)){
  site.temp <- sites.fun[k]
  surveys.temp <- site.list.fun$unique_id[site.list.fun$Abbrev == site.temp]
  
  cbc_temp <- (cbc_cast[row.names(cbc_cast) %in% surveys.temp,])
  
  if (is.matrix(cbc_temp)){
    betas.temp.prime <- new.comdistnn.prime(my.sample = cbc_temp, my.dist.mat = mat2)
    
    for (q in 1:nrow(betas.temp.prime)){
      temp.unique.id <- row.names(betas.temp.prime)[q]
      betas.out.prime[betas.out.prime$second_bin == temp.unique.id,"beta.from.baseline"] <- betas.temp.prime[q,1]}
    
  }
  
}  


fl.merge <- merge(fl.merge, betas.out, by.x = "site_bin", by.y = "second_bin")
mod2 <- lmer(data = fl.merge, beta.from.baseline ~ base.pct.diff + yr_bin + (1|Abbrev))
summary(mod2)
coeffs <- coef(summary(mod2))
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2
cbind(coeffs, "p value" = round(p,3))
coeff_abbrev <- coef(mod2)$Abbrev

#lefse code
#NN community distance without abundance weighting 

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

#Prime methods utilize community abundances in the calculations 
new.comdistnn.prime = function(my.sample, my.dist.mat){
  
  my.ra.sample = my.sample/rowSums(my.sample)
  
  get.presents = function(x){
    names(x[x>0])
  }
  
  get.weights = function(x){
    x[x>0]
  }
  
  
  list.of.weights = apply(my.ra.sample, 1, get.weights)
  
  Dnn.apply.function = function(x){
    tmp.function = function(z){
      
      weighted.mean(c(apply(as.matrix(my.dist.mat[names(x), names(z)]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[names(x), names(z)]), MARGIN=2, min, na.rm=T)),
                    c(x,z)
      )		
      
    }
    
    lapply(list.of.weights, FUN=tmp.function)
    
  }
  
  dnn.output = lapply(list.of.weights, Dnn.apply.function)
  
  outt = do.call(rbind,lapply(dnn.output,unlist))
  
  
  outt
  
}

#NN Prime Null Model 
my.sample = cbc_tmp2
my.dist.mat = mat2
new.comdistnn.prime.null = function(my.sample, my.dist.mat){
  
  row.names(my.dist.mat) = sample(row.names(my.dist.mat), length(row.names(my.dist.mat)), replace=F)
  
  my.ra.sample = my.sample/rowSums(my.sample)
  
  get.presents = function(x){
    names(x[x>0])
  }
  
  get.weights = function(x){
    as.list(x[x>0])
  }
  
  
  list.of.weights = apply(my.ra.sample, 1, get.weights)
  
  Dnn.apply.function = function(x){
    tmp.function = function(z){
      
      weighted.mean(c(apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=2, min, na.rm=T)),
                    c(x,z)
      )		
      
    }
    
    lapply(list.of.weights, FUN=tmp.function)
    
  }
  
  dnn.output = lapply(list.of.weights, Dnn.apply.function)
  
  
  dnn.output = do.call(rbind,lapply(dnn.output,unlist))
  
}
t <- new.comdistnn.prime.null(my.sample, my.dist.mat)

##With null prime code 
#Need to figure out a way to generate means of replicates and 
#store the mean values for each replicate 
for (r in 1:length(sites.fun)){
  site.temp <- sites.fun[r]
  surveys.temp <- site.list.fun$unique_id[site.list.fun$Abbrev == site.temp]
  
  cbc_temp <- cbc_cast[row.names(cbc_cast) %in% surveys.temp,]
  
  if (is.matrix(cbc_temp)){
    #betas.temp.prime.null <-  cbc_temp[c(1, r+1),] #How to tell R to take the next recurring
    #new.comdistnn.prime.null(my.sample = cbc_temp, my.dist.mat = mat2)
    
    for (j in 1:(nrow(cbc_temp) - 1)){
      temp.unique.id <- row.names(cbc_temp)[j+1]
      cbc_tmp2 <- cbc_temp[c(1,j + 1),]
      null.tmp <- vector(length = 1000)
      for (x in 1:1000){
        tmp <- new.comdistnn.prime.null(cbc_tmp2, mat2)
        null.tmp[x] <- tmp
      }
      mean.null <- mean(null.tmp)
      sd.null <- sd(null.tmp)
      #matrix[j,] <- data.frame(row.names = temp.unique.id, mean.null, sd.null)
      #clean this up how save mean and sd column
      betas.out.prime.null[betas.out.prime.null$second_bin == temp.unique.id, "mean.from.null"] <- mean.null
      betas.out.prime.null[betas.out.prime.null$second_bin == temp.unique.id, "sd.from.null"] <- sd.null
      }
  
  }
  
}

beta.merge <- merge(betas.out.prime.null, betas.out.prime)
beta.merge[is.na(beta.merge)] <- 0
beta.merge <- separate(beta.merge, second_bin, c("Abbrev", "yr_bin"), sep = "_")
beta.merge$unique_id <- paste0(beta.merge$Abbrev, "_", beta.merge$yr_bin)
beta.merge.agg <- aggregate(list(null.mean = beta.merge$mean.from.null, sd.null = beta.merge$sd.from.null, beta = beta.merge$beta.from.baseline), by = list(Abbrev = beta.merge$Abbrev), mean) 
beta.merge.agg <- mutate(beta.merge.agg, SES = (beta - null.mean) / sd.null)
beta.merge.agg <- beta.merge.agg[complete.cases(beta.merge.agg), ]
#betas.temp.prime.null[j+1, 1:2] this was in place of "matrix[j, 1:2]"

betas.temp.prime.null
new.comdistnn.prime.null(cbc_tmp2, mat2)
new.comdistnn.prime.null(cbc_tmp2, mat2)
for (q in 1:length(list.of.weights)){
  print(paste("first site", q))
  p <- list.of.weights[[q]]
  for (l in 1:length(list.of.weights)){
    print(l)
    m <- list.of.weights[[l]]

    temp <- weighted.mean(c(apply(as.matrix(my.dist.mat[names(p), names(m)]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[names(p), names(m)]), MARGIN=2, min, na.rm=T)),
                  c(p,m))
    print(temp)
  }
}
#   





####
cbc_red <- select(cbc_binned_agg, unique_id, Abbrev, yr_bin)

fun_beta_merge <- merge(betas.out.prime, cbc_red, by.x = "second_bin", by.y = "unique_id")

n.sites <- length(unique(fun_beta_merge$Abbrev))
site.list <- as.character(unique(fun_beta_merge$Abbrev))

#Calculate slopes for B diversity
slopes_sites_beta <- data.frame(site.list, beta.slope = NA)

for (i in 1:n.sites){
  site.temp <- site.list[i]
  fun_beta_merge_temp <- fun_beta_merge[fun_beta_merge$Abbrev == site.temp,]
  fun_beta_merge_temp <- fun_beta_merge_temp[complete.cases(fun_beta_merge_temp),]
  if (nrow(fun_beta_merge_temp) > 1){
    lm.temp <- lm(fun_beta_merge_temp$beta.from.baseline ~ fun_beta_merge_temp$yr_bin)
    slope.temp <- summary(lm.temp)$coefficients[2,1]
    slopes_sites_beta[i,2] <- slope.temp
  }
}

t.test(slopes_sites_beta$beta.slope, mu = 0)

cbc_abbrev <- unique(cbc[, c("Abbrev", "Subnational_code")])

cbc_abbrev$Abbrev <- as.character(cbc_abbrev$Abbrev)
fun_beta_merge$Abbrev <- as.character(fun_beta_merge$Abbrev)
fun_beta_merge <- merge(fun_beta_merge, cbc_abbrev, by = "Abbrev")
fun_beta_merge$yrs <- fun_beta_merge$yr_bin
fun_beta_merge$yrs[fun_beta_merge$yrs == 1] <- "1980-1984"
fun_beta_merge$yrs[fun_beta_merge$yrs == 2] <- "1985-1989"
fun_beta_merge$yrs[fun_beta_merge$yrs == 3] <- "1990-1994"
fun_beta_merge$yrs[fun_beta_merge$yrs == 4] <- "1995-1999"
fun_beta_merge$yrs[fun_beta_merge$yrs == 5] <- "2000-2004"
fun_beta_merge$yrs[fun_beta_merge$yrs == 6] <- "2005-2009"
fun_beta_merge$yrs[fun_beta_merge$yrs == 7] <- "2010-2014"
fun_beta_merge$yrs[fun_beta_merge$yrs == 8] <- "2015-2017"

gghist_beta <- ggplot(data = slopes_sites_beta, aes(slopes_sites_beta$beta.slope)) + 
  geom_histogram(col = "black", fill = "black", bins = 20, binwidth = NULL) + 
  labs(title = "") +
  labs(x = expression(paste("Slopes of Functional ", beta, "-diversity")), y = "# of Sites") +
  theme_cowplot(font_size = 14, line_size = 1.2) +
  coord_flip()

plot_beta <- (ggplot(fun_beta_merge, aes(x = yrs, y = beta.from.baseline, group = Abbrev, 
                                          colour = factor(Subnational_code, 
                                                          labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                #geom_point(size = 3) +
                #geom_line()+
                geom_smooth(method = lm, se = FALSE, aes(x = yrs, y = beta.from.baseline, group = Abbrev, 
                                                         colour = factor(Subnational_code, 
                                                                         labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                xlab("Five-year Bins") +
                ylab(expression(paste("Functional ", beta, "-diversity"))) +
                labs(colour = "States") +
                #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                theme_bw() +
                theme(axis.line = element_line(colour = "black", size =1.2),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12),
                      axis.title.x = element_text(vjust = -1, size = 14),
                      axis.title.y = element_text(vjust = 1.5, size = 14),
                      axis.ticks = element_line(size = 1.2),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,2,2), "lines"),
                      text = element_text(size=14)))

ggdraw() + 
  draw_plot(plot_beta + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  draw_plot(gghist_beta, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 

