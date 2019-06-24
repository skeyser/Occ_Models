###Summary Stats for models###
means.output <- upm1$mean
Rhats <- upm1$Rhat
model.sum <- upm1$summary

#rm("upm1")

J <- 274
K <- 38
S <- max(spp.occ$spp.id)
#Mean detection across species 
mean.p.sp <- data.frame(matrix(NA, nrow = 5, ncol = 1))
names(mean.p.sp) <- c("mean.p")

#Loop throug all species
for( s in 1:length(unique(spp.occ$spp.id)) ){ 
  #Pull out one species matrix * by non-surveyed sites and years
  p.sp.mat <- (means.output$p[ s, 1:J, 1:K ]  * JKsurv[ 1:J, 1:K ])
  #Sum across sites for each year
  p.sp.mat <- apply(p.sp.mat, 2, sum)
  #divide by total number of sampled segments each year
  p.sp <-  p.sp.mat[ 1:K ] / surveyedJ[ 1:K ]
  #Put values into empty DF
  mean.p.sp [s,] <- mean(p.sp)
}#S

#Pull out just the z matrix
z.prime <- means.output$z

#Loop through each species, site, and year to set 0s and 1s 
#manually using cut-off values

for (s in 1:S){
  for (j in 1:J){
    for (k in 1:K){
      #z.prime <- ifelse( mean$z[s, j, k] >= 0.50, 1, 0) 
      #round function rounds to nearest interger so it would do the same as above line
      z.prime[s, j, k] <- ifelse(means.output$z[s, j, k] > .5, 1, 0) 
    }#K
  }#J
}#S

#Multiply each species matrix by the NA df for species totals
for( i in 1:S ){
  z.prime[ i, , ] <- z.prime[ i, , ] * JKmat
}

z.prime[1,,]
z.prime[2,,]
sum(is.na(z.prime[1,,]))

#Find yearly alpha diversity per subgroup
group.a.div <- matrix(NA, 274, 38)

for (j in 1:J){
  for ( k in  1:K){
    #zeros those z for unsampled segments
    group.a.div[ j, k ] <- sum( z.prime[ 1:S, j, k ] * JKsurv[ j, k ] )
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

         