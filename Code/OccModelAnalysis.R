##########################################################
#########Single Visit Multispecies Occupancy Model########
######Script has code for JAGS specified model############
###############Created By: Spencer R Keyser###############
##########################################################

######Load Packages######
#rm(list = ls())

#library(pacman)
#pacman::p_load(tidyverse, ggplot2, ggfortify, JagsUI)

#install.packages("jagsUI")
library("jagsUI")

#########Package Loading Complete########
#########################################

###Load in the Data4OccModels Workspace###
##########################################

#load(here::here('R Workspace/Data4OccModels_6_3.RData'))


#########################################################################################
###########################################################################
####################### define MCMC settings ##############################

ni <- 2; nt <- 1; nb <- 0; nc <- 1 #iterations, thinning, burnin, chains

##### end of MCMC parameters definition ##############
############################################################################
#################    run alternative models ################################

####################################################################
#

##########################################
##### Multi-species occ model single######
#####             visit            #######
#####site-level characterisitics   #######
#####        of time of day        #######
#####and first-year observer effect#######
#####Random intercepts for species #######
#####         and for sites        #######
##########################################

###########################################
###Load in Data###
#bcr.occ <- read.csv(bcr_occ)

# jdf <- read.csv(jdf)
# 
# spp.occ <- read.csv(spp_occ)
# 
# TOD.ma <- read.csv(TOD_ma)
# TOD.ma <- as.matrix(TOD.ma)
# 
# Ord.ma <- read.csv(Ord_ma)
# Ord.ma <- as.matrix(Ord.ma)
# # 
# Obs.ma <- read.csv(Ord_ma)
# Obs.ma <- as.matrix(Obs.ma)
# # 
# load(ydf)
#### end data load ####
########Model specification###############
sink( "om1.txt" )
cat("
    model{
    
    #ecological model
    #for occupancy
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #relate occupancy probability to random intercepts for route and year:
    logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
    eps.psi[ k ]
    for( s in 1:S ){ #loop species
    #modeling true occupancy
    muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
    z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
    }#S
    }#K
    }#J
    
    #Detection Model 
    for (s in 1:S){ #loop species
    for (j in 1:J){ #loop sites (rt_segment)
    for (k in 1:K){ #loop years 
    
    logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
    alpha[2] * Ord.ma[j,k] + 
    alpha[3] * TOD.ma[j, k] + 
    alpha[ 4 ] * Mass.scaled[ s ] +
    delta[ s ]
    
    mup[ s, j, k ] <- z[s, j, k] * p[s, j, k]
    
    ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
    
    }#K
    }#J
    }#S
    
    #defining data augmentation indicator prior
    for( s in 1:S ){ #loop species
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #restrict data augmentation indicator to the corresponding BCR species      
    w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
    
    #indicator of whether species should be added to a particular segment
    w[ s, j, k ] ~ dbern(  omega )
    } #k
    } #j
    } #s
    
    
    #priors
    #data augmentation prior
    omega ~ dbeta( 4, 4 )
    
    
    #priors for detection model 
    
    #Defining the intercept prior for alpha.0
    int.p <- log(mean.p / (1 - mean.p))
    mean.p ~ dbeta(4, 4)
    
    for(q in 1:Q){ #loop through predictors
    alpha[q] ~ dnorm(0, 0.1)
    }#Q
    
    # Priors for random species intercepts
    for(g in 1:G ){
    #phylo group mean response and error
    mu.spp[ g ] ~ dnorm(0, prec.spp[ g ] )
    prec.spp[ g ] <- 1 / ( sigma.spp[ g ] * sigma.spp[ g ] )
    sigma.spp[ g ] ~ dt(0, 1, 4) T(0, )
    } #g
    
    #loop through species
    for( s in 1:S ){
    #estimate species-specific response based on phylo group mean
    delta[ s ] ~ dnorm( mu.spp[ Phylo.V1.code[s] ], 
    prec.spp[ Phylo.V1.code[s] ] )
    } #s
    
    
    #define intercepts for ecological model as mean probs:
    #logit main intercept
    int.psi <- log(mean.psi / (1 - mean.psi))
    #mean occupancy probability
    mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
    
    #random year intercepts
    for( k in 1:K ){
    eps.psi[k] ~ dnorm( 0, prec.psi ) T( -6, 6 )
    } #K loop
    
    #associated precision of radom year intercepts:
    prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
    sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
    
    #random route effect
    for( m in 1:M ){#loop routes
    epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
    }#M
    
    #associated precision for random route effect
    precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
    sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
    
    
    # 
    # #derived parameters 
    #   for( k in 1:K ){
    #     for( j in 1:J ){
    #       #estimate yearly alpha diversity
    #       a.div[ j, k ] <- sum( z[ 1:S, j, k ] * JKmat[ j, k ] )
    #     } #J
    #   }#K
    # 
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "om1.txt"

#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             #, 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             #             , 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             , 'delta' #random intercept for phylogenetic group
             #             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
fm1 <- jags( win.data, inits = inits, params, modelname, 
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE) 
fm1.time <- proc.time() - ptm

# #auto update the model
# upm1 <- autojags( win.data, inits = inits, params, modelname, 
#           n.chains = 3, n.thin = 5, n.burnin = 0,
#           iter.increment = 10, max.iter = 100,
#           Rhat.limit = 1.1, save.all.iter=FALSE, parallel = TRUE )




##############################################################################
######### different data augmentation set up ##############################
########Model specification###############
sink( "om2.txt" )
cat("
    model{
    
    #ecological model
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #relate occupancy probability to random intercepts for route and year:
    logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
    eps.psi[ k ]
    for( s in 1:S ){ #loop species
    #modeling true occupancy
    muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
    z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
    }#S
    }#K
    }#J
    
    #Detection Model 
    for (s in 1:S){ #loop species
    for (j in 1:J){ #loop sites (rt_segment)
    for (k in 1:K){ #loop years 
    #relate p to predictors
    logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
    alpha[2] * Ord.ma[j,k] + 
    alpha[3] * TOD.ma[j, k] + 
    alpha[ 4 ] * Mass.scaled[ s ] +
    delta[ s ]
    mup[ s, j, k ] <- z[s, j, k] * p[s, j, k]
    ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
    }#K
    }#J
    }#S
    
    #defining data augmentation indicator prior
    for( s in 1:S ){ #loop species
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #restrict data augmentation only if undetected species is found in BCR      
    w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
    #indicator of whether undetected species should be added to j segment and k year
    w[ s, j, k ] ~ dbern(  omega.y )
    #only use omega for undetected species:
    omega.y <- ifelse( det[ s, j, k ] == 1, 1, omega )
     } #k
    } #j
    } #s
    
    #priors
    #data augmentation prior
    omega ~ dbeta(4,4)

    #priors for detection model 
    #Defining the intercept prior for alpha.0
    int.p <- log(mean.p / (1 - mean.p))
    mean.p ~ dbeta(4, 4)
    
    for(q in 1:Q){ #loop through predictors
    alpha[q] ~ dnorm(0, 0.1)
    }
    
    # Priors for random species intercepts
    for(g in 1:G ){
    #phylo group mean response and error
    mu.spp[ g ] ~ dnorm(0, prec.spp[ g ] )
    prec.spp[ g ] <- 1 / ( sigma.spp[ g ] * sigma.spp[ g ] )
    sigma.spp[ g ] ~ dt(0, 1, 4) T(0, )
    } #g
    
    #loop through species
    for( s in 1:S ){
    #estimate species-specific response based on phylo group mean
    delta[ s ] ~ dnorm( mu.spp[ Phylo.V1.code[s] ], 
    prec.spp[ Phylo.V1.code[s] ] )
    } #s
    
    
    #define intercepts for ecological model as mean probs:
    #logit main intercept
    int.psi <- log(mean.psi / (1 - mean.psi))
    #mean occupancy probability
    mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
    
    #random year intercepts
    for( k in 1:K ){
    eps.psi[k] ~ dnorm( 0, prec.psi ) T( -6, 6 )
    } #K loop
    
    #associated precision of radom year intercepts:
    prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
    sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
    
    #random route effect
    for( m in 1:M ){#loop routes
    epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
    }#M
    
    #associated precision for random route effect
    precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
    sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
    
    
    # 
    # #derived parameters 
    #   for( k in 1:K ){
    #     for( j in 1:J ){
    #       #estimate yearly alpha diversity
    #       a.div[ j, k ] <- sum( z[ 1:S, j, k ] * JKmat[ j, k ] )
    #     } #J
    #   }#K
    # 
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "om2.txt"

#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             , 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             #, 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             , 'delta' #random intercept for phylogenetic group
             #             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       det = ydf, #to avoid cyclic error
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
fm2 <- jags( win.data, inits = inits, params, modelname, 
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE) 
fm2.time <- proc.time() - ptm

# #auto update the model
# upm2 <- autojags( win.data, inits = inits, params, modelname, 
#                   n.chains = 3, n.thin = 5, n.burnin = 0,
#                   iter.increment = 10, max.iter = 100,
#                   Rhat.limit = 1.1, save.all.iter=FALSE, parallel = TRUE )







##############################################################################
######### from om2: commenting out different aspects ##############################
########Model specification###############
sink( "om3.txt" )
cat("
    model{
    
    #ecological model
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #relate occupancy probability to random intercepts for route and year:
    logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
                            eps.psi[ k ]
    for( s in 1:S ){ #loop species
    #modeling true occupancy
    muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
    #muz[ s, j, k ] <- psi[ j, k ] * bcr.occ[ s, bcr.id[j] ]
    z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
    }#S
    }#K
    }#J
    
    #Detection Model 
    for (s in 1:S){ #loop species
    for (j in 1:J){ #loop sites (rt_segment)
    for (k in 1:K){ #loop years 
    #relate p to predictors
    logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
    alpha[2] * Ord.ma[j,k] + 
    alpha[3] * TOD.ma[j, k] + 
    alpha[ 4 ] * Mass.scaled[ s ] #+
    # delta[ s ]
    mup[ s, j, k ] <- z[s, j, k] * p[s, j, k]
    ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
    }#K
    }#J
    }#S
    
    #defining data augmentation indicator prior
    for( s in 1:S ){ #loop species
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #restrict data augmentation only if undetected species is found in BCR
    w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
    #indicator of whether undetected species should be added to j segment and k year
    w[ s, j, k ] ~ dbern(  omega.y )
    #only use omega for undetected species:
    omega.y <- ifelse( det[ s, j, k ] == 1, 1, omega )
    } #k
    } #j
    } #s
    
    #priors
    #data augmentation prior
    omega ~ dbeta(4,4)

    #priors for detection model 
    #Defining the intercept prior for alpha.0
    int.p <- log(mean.p / (1 - mean.p))
    mean.p ~ dbeta(4, 4)
    
    for(q in 1:Q){ #loop through p predictors
    alpha[q] ~ dnorm(0, 0.1)
    }
    
    # # Priors for random species intercepts
    # for(g in 1:G ){
    #   #phylo group mean response and error
    #   mu.spp[ g ] ~ dnorm(0, prec.spp[ g ] )
    #   prec.spp[ g ] <- 1 / ( sigma.spp[ g ] * sigma.spp[ g ] )
    #   sigma.spp[ g ] ~ dt(0, 1, 4) T(0, )
    # } #g
    # 
    #   #loop through species
    #   for( s in 1:S ){
    #       #estimate species-specific response based on phylo group mean
    #       delta[ s ] ~ dnorm( mu.spp[ Phylo.V1.code[s] ],
    #                                 prec.spp[ Phylo.V1.code[s] ] )
    #     } #s
    
    #define intercepts for ecological model as mean probs:
    #logit main intercept
    int.psi <- log(mean.psi / (1 - mean.psi))
    #mean occupancy probability
    mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
    
    #random year intercepts
    for( k in 1:K ){
    eps.psi[k] ~ dnorm( 0, prec.psi ) T( -6, 6 )
    } #K loop
    
    #associated precision of radom year intercepts:
    prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
    sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
    
    #random route effect
    for( m in 1:M ){#loop routes
    epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
    }#M
    
    #associated precision for random route effect
    precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
    sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
    
    
    # 
    # #derived parameters 
    #   for( k in 1:K ){
    #     for( j in 1:J ){
    #       #estimate yearly alpha diversity
    #       a.div[ j, k ] <- sum( z[ 1:S, j, k ] * JKmat[ j, k ] )
    #     } #J
    #   }#K
    # 
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "om3.txt"

#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             , 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             #, 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             #, 'delta' #random intercept for phylogenetic group
             #             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       det = ydf, #to avoid cyclic error
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
fm3 <- jags( win.data, inits = inits, params, modelname, 
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE ) 
fm3.time <- proc.time() - ptm
#fm3 with delta out




######### from om2: commenting out different aspects ##############################
########Model specification###############
sink( "om4.txt" )
cat("
    model{
    
    #ecological model
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #relate occupancy probability to random intercepts for route and year:
    logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
    eps.psi[ k ]
    for( s in 1:S ){ #loop species
    #modeling true occupancy
    #muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
    muz[ s, j, k ] <- psi[ j, k ] * bcr.occ[ s, bcr.id[j] ]
    z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
    }#S
    }#K
    }#J
    
    #Detection Model 
    for (s in 1:S){ #loop species
    for (j in 1:J){ #loop sites (rt_segment)
    for (k in 1:K){ #loop years 
    #relate p to predictors
    logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
    alpha[2] * Ord.ma[j,k] + 
    alpha[3] * TOD.ma[j, k] + 
    alpha[ 4 ] * Mass.scaled[ s ] +
    delta[ s ]
    mup[ s, j, k ] <- z[s, j, k] * p[s, j, k]
    ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
    }#K
    }#J
    }#S
    
    #  #defining data augmentation indicator prior
    # for( s in 1:S ){ #loop species
    #  for( j in 1:J ){ #loop segments
    #    for( k in 1:K ){ #loop years
    #      #restrict data augmentation only if undetected species is found in BCR      
    #      w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
    #      #indicator of whether undetected species should be added to j segment and k year
    #      w[ s, j, k ] ~ dbern(  omega.y )
    #      #only use omega for undetected species:
    #      omega.y <- ifelse( y[ s, j, k ] == 1, 1, omega )
    #     } #k
    #   } #j
    #  } #s
    #  
    #  #priors
    #  #data augmentation prior
    #  omega ~ dbeta(4,4)
    
    
    #priors for detection model 
    #Defining the intercept prior for alpha.0
    int.p <- log(mean.p / (1 - mean.p))
    mean.p ~ dbeta(4, 4)
    
    for(q in 1:Q){ #loop through p predictors
    alpha[q] ~ dnorm(0, 0.1)
    }
    
    # Priors for random species intercepts
    for(g in 1:G ){
    #phylo group mean response and error
    mu.spp[ g ] ~ dnorm(0, prec.spp[ g ] )
    prec.spp[ g ] <- 1 / ( sigma.spp[ g ] * sigma.spp[ g ] )
    sigma.spp[ g ] ~ dt(0, 1, 4) T(0, )
    } #g
    
    #loop through species
    for( s in 1:S ){
    #estimate species-specific response based on phylo group mean
    delta[ s ] ~ dnorm( mu.spp[ Phylo.V1.code[s] ],
    prec.spp[ Phylo.V1.code[s] ] )
    } #s
    
    #define intercepts for ecological model as mean probs:
    #logit main intercept
    int.psi <- log(mean.psi / (1 - mean.psi))
    #mean occupancy probability
    mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
    
    #random year intercepts
    for( k in 1:K ){
    eps.psi[k] ~ dnorm( 0, prec.psi ) T( -6, 6 )
    } #K loop
    
    #associated precision of radom year intercepts:
    prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
    sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
    
    #random route effect
    for( m in 1:M ){#loop routes
    epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
    }#M
    
    #associated precision for random route effect
    precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
    sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
    
    
    # 
    # #derived parameters 
    #   for( k in 1:K ){
    #     for( j in 1:J ){
    #       #estimate yearly alpha diversity
    #       a.div[ j, k ] <- sum( z[ 1:S, j, k ] * JKmat[ j, k ] )
    #     } #J
    #   }#K
    # 
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "om4.txt"

#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             , 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             #, 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             , 'delta' #random intercept for phylogenetic group
             #             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
fm4 <- jags( win.data, inits = inits, params, modelname, 
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE) 
fm4.time <- proc.time() - ptm

print(fm4)
#fm4 with data augmentation out





##############################################################################
######### different data augmentation set up ##############################
########Model specification###############
sink( "om5.txt" )
cat("
    model{
    
    #ecological model
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #relate occupancy probability to random intercepts for route and year:
    logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
    eps.psi[ k ]
    for( s in 1:S ){ #loop species
    #modeling true occupancy
    muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
    #only estimate z for surveyed sites and segments
    muzsurv[ s, j, k ] <- ifelse( det[ s, j, k ] >= 0, muz[ s, j, k ], 0 )
    z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
    }#S
    }#K
    }#J
    
    #Detection Model 
    for (s in 1:S){ #loop species
    for (j in 1:J){ #loop sites (rt_segment)
    for (k in 1:K){ #loop years 
    #relate p to predictors
    logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
    alpha[2] * Ord.ma[j,k] + 
    alpha[3] * TOD.ma[j, k] + 
    alpha[ 4 ] * Mass.scaled[ s ] +
    delta[ s ]
    #only use p for surveyed sites
    psurv[ s, j, k ] <- ifelse( det[ s, j, k ] >= 0, p[ s, j, k ], 0 )
    mup[ s, j, k ] <- z[s, j, k] * psurv[s, j, k]
    
    ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
    }#K
    }#J
    }#S
    
    #defining data augmentation indicator prior
    for( s in 1:S ){ #loop species
    for( j in 1:J ){ #loop segments
    for( k in 1:K ){ #loop years
    #restrict data augmentation only if undetected species is found in BCR      
    w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
    #indicator of whether undetected species should be added to j segment and k year
    w[ s, j, k ] ~ dbern(  omega.y )
    #only use omega for undetected species:
    omega.y <- ifelse( det[ s, j, k ] == 1, 1, omega )
    } #k
    } #j
    } #s
    
    #data augmentation prior
    omega ~ dbeta(4,4)
    
    #priors
    #priors for detection model 
    
    #Defining the intercept prior for alpha.0
    int.p <- log(mean.p / (1 - mean.p))
    mean.p ~ dbeta(4, 4)
    
    for(q in 1:Q){ #loop through predictors
    alpha[q] ~ dnorm(0, 0.1)
    }
    
    # Priors for random species intercepts
    for(g in 1:G ){
    #phylo group mean response and error
    mu.spp[ g ] ~ dnorm(0, prec.spp[ g ] )
    prec.spp[ g ] <- 1 / ( sigma.spp[ g ] * sigma.spp[ g ] )
    sigma.spp[ g ] ~ dt(0, 1, 4) T(0, )
    } #g
    
    #loop through species
    for( s in 1:S ){
    #estimate species-specific response based on phylo group mean
    delta[ s ] ~ dnorm( mu.spp[ Phylo.V1.code[s] ], 
    prec.spp[ Phylo.V1.code[s] ] )
    } #s
    
    
    #define intercepts for ecological model as mean probs:
    #logit main intercept
    int.psi <- log(mean.psi / (1 - mean.psi))
    #mean occupancy probability
    mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
    
    #random year intercepts
    for( k in 1:K ){
    eps.psi[k] ~ dnorm( 0, prec.psi ) T( -6, 6 )
    } #K loop
    
    #associated precision of radom year intercepts:
    prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
    sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
    
    #random route effect
    for( m in 1:M ){#loop routes
    epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
    }#M
    
    #associated precision for random route effect
    precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
    sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
    
    
    # 
    # #derived parameters 
    #   for( k in 1:K ){
    #     for( j in 1:J ){
    #       #estimate yearly alpha diversity
    #       a.div[ j, k ] <- sum( z[ 1:S, j, k ] * JKmat[ j, k ] )
    #     } #J
    #   }#K
    # 
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "om5.txt"

#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             , 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             #, 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             , 'delta' #random intercept for phylogenetic group
             #             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       det = ydf, #to avoid cyclic error
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
fm5 <- jags( win.data, inits = inits, params, modelname, 
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE) 
fm5.time <- proc.time() - ptm
#fm5 trying to avoid estimation of z and y for unsampled jXk

# #auto update the model
# upm2 <- autojags( win.data, inits = inits, params, modelname, 
#                   n.chains = 3, n.thin = 5, n.burnin = 0,
#                   iter.increment = 10, max.iter = 100,
#                   Rhat.limit = 1.1, save.all.iter=FALSE, parallel = TRUE )


################ end of script ##########################################