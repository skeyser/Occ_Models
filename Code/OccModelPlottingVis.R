###Summary Stats for models###
means.output <- upm1$mean
Rhats <- upm1$Rhat
model.sum <- upm1$summary

#rm("upm1")

J <- 236
K <- 38
S <- max(spp.occ$spp.id)
#Mean detection across species 
mean.p.sp <- data.frame(Species = 1:length(unique(spp.occ$Species)), "Mean p No RE" = NA,
                        "Mean p RE" = NA)

mean.p.sp$Species <- spp.occ$Species
#Loop throug all species
for( s in 1:length(unique(spp.occ$spp.id)) ){ 
  #Pull out one species matrix * by non-surveyed sites and years
  p.sp.mat <- (means1$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
  p.sp.mat2 <- (means2$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
  #Sum across sites for each year
  p.sp.mat <- apply(p.sp.mat, 2, sum, na.rm = T)
  p.sp.mat2 <- apply(p.sp.mat2, 2, sum, na.rm = T)
  #divide by total number of sampled segments each year
  p.sp <-  p.sp.mat[ 1:K ] / surveyedJ[ 1:K ]
  p.sp2 <-  p.sp.mat2[ 1:K ] / surveyedJ[ 1:K ]
  #Put values into empty DF
  mean.p.sp [s, 2] <- mean(p.sp)
  mean.p.sp [s, 3] <- mean(p.sp2)
}#S

#Pull out just the z matrix
z.prime1 <- means1$z
z.prime2 <-means2$z

#Loop through each species, site, and year to set 0s and 1s 
#manually using cut-off values

for (s in 1:S){
  for (j in 1:J){
    for (k in 1:K){
      #z.prime <- ifelse( mean$z[s, j, k] >= 0.50, 1, 0) 
      #round function rounds to nearest interger so it would do the same as above line
      z.prime1[s, j, k] <- ifelse(means1$z[s, j, k] > .5, 1, 0)
      z.prime2[s, j, k] <- ifelse(means2$z[s, j, k] > .5, 1, 0) 
    }#K
  }#J
}#S

#Multiply each species matrix by the NA df for species totals
for( i in 1:S ){
  z.prime1[ i, , ] <- z.prime1[ i, , ] * JKsurv
  z.prime2[ i, , ] <- z.prime2[ i, , ] * JKsurv
}

z.prime[1,,]
z.prime[2,,]
sum(is.na(z.prime[1,,]))

#Find yearly alpha diversity per subgroup
group.a.div <- matrix(NA, 274, 38)
group.a.div2 <- matrix(NA, 274, 38)

for (j in 1:J){
  for ( k in  1:K){
    #zeros those z for unsampled segments
    group.a.div[ j, k ] <- sum( z.prime1[ 1:S, j, k ] * JKsurv[ j, k ] )
    group.a.div2[ j, k ] <- sum( z.prime2[ 1:S, j, k ] * JKsurv[ j, k ] )
    #keep only average estimate
    #mean.a.div[ j, k ] <- mean(a.div[ j, k ])
  }#K
}#J


#Quick function to calculate total observations of each species
total <- function(x, y){
  df <- x[y,,]
  yr.sums <- colSums(df, 1:length(nrow))
  sum(yr.sums)
}

total.observed <- data.frame("Observed" = 1:max(spp.occ$spp.id), "Estimated RE" = NA, 
                             "Estimated No Sp RE" = NA)

for (i in 1:max(spp.occ$spp.id)){
  total.observed[i, 1] <- total(ydf, i)
  total.observed[i, 2] <- total(z.prime1, i)
  total.observed[i, 3] <- total(z.prime2, i)
}

#View All non-converged parameters 
bad.params <- function(x) {
  rhats <- unlist(Rhats)
  buggers <- rhats[rhats > 1.15]
  buggers[complete.cases(buggers)]
}

noconverge <- bad.params(Rhats)

###########################################################################
####################################Plotting###############################
###########################################################################

##Traceplots for relevant parameters

MCMCtrace(upm1, c("alpha", "delta", "int.psi", "int.p", "epsID.psi", "sigma.delta", "sigma.psi", "sigmaID.psi"), Rhat = T, pdf = T, file = here::here("Figues and Tables/TraceSwallows"))


#Parameter Plots 
MCMCplot(upm1, params = c("alpha"), rank = T, main = "Swallows Phylo Group Alphas", labels = c("Observer Effect", "Time of Day",
                                                                                               "Day of Year", "Body Mass"))

MCMCplot(upm1, params = c("int.psi", "int.p"), rank = T, main = "Swallows Phylo Group Intercepts", labels = c("Intercept Occupancy", 
                                                                                                              "Intercept Detection")) 
MCMCplot(upm1, params = c("epsID.psi", "sigmaID.psi"), main = "Swallows Phylo Group Site Effects")
                                                                                                                
MCMCplot(upm1, params = c("eps.psi", "sigma.psi"), main = "Swallows Phylo Group Year Effects")

MCMCplot(upm1, params = c("delta", "sigma.delta"), main = "Swallows Species RE", labels = c("BASW", "CASW", "CLSW", 
                                                                                            "PUMA", "BASW", "NRWS", 
                                                                                            "Sigma Delta"))


#Save workspace
#save.image(file = here::here("R Workspace/SwallowModelOut.RData"))

         