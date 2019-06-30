##########################################################
#########Single Visit Multispecies Occupancy Model########
######Script has code for JAGS specified model############
###############Created By: Spencer R Keyser###############
##########################################################

######Load Packages######
#rm(list = ls())
#install.packages("pacman")
library("pacman")
pacman::p_load(jagsUI, MCMCvis, here)

#########Package Loading Complete########
#########################################

#As per usual. Always checking github

###Load in the Data4OccModels Workspace###
##########################################

load(here::here('R Workspace/Phylo_Subs/Blackbirds.RData'))

##Hello github##


#########################################################################################
###########################################################################
####################### define MCMC settings ##############################

ni <- 100; nt <- 1; nb <- 0; nc <- 3 #iterations, thinning, burnin, chains
#ni <- 5000; nt <- 10; nb <- 1000; nc <- 3 #iterations, thinning, burnin, chains

#do these need to be on this script?
J <- length(unique(jdf$site.id))
M <- length(unique(jdf$rteno.id))
S <- length(unique(spp.occ$spp.id))
K <- 38

#Correcting the dimensions of the bcr.occ
bcr.occ <- bcr.occ[1:S, ]

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

########################################################################
### full model for single phylo group - #######
### but excluding delta[s]: species random intercept in detection model #
#############################################################################
########Model specification###############
# sink( "pm1.txt" )
# cat("
#     model{
#     
#     #ecological model
#     #for occupancy
#     for( j in 1:J ){ #loop segments
#       for( k in 1:K ){ #loop years
#       #relate occupancy probability to random intercepts for route and year:
#       logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
#                               eps.psi[ k ]
#         for( s in 1:S ){ #loop species
#         #modeling true occupancy
#         muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j, k ]
#         z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
#         }#S
#       }#K
#     }#J
#     
#     #Detection Model 
#     for (s in 1:S){ #loop species
#      for (j in 1:J){ #loop sites (rt_segment)
#         for (k in 1:K){ #loop years 
#           #relate detection probability to siteXyear, and species predictors
#           logit( p[s,j,k] ) <- int.p + alpha[1] * Obs.ma[j,k] +
#                                 alpha[2] * Ord.ma[j,k] + 
#                                 alpha[3] * TOD.ma[j, k] + 
#                                 alpha[ 4 ] * Mass.scaled[ s ] 
#           mup[ s, j, k ] <- z[s, j, k] * p[s, j, k]
#           #relate observed occupancy to detection prob and estimated occupancy.
#           ydf[s, j, k] ~ dbern( mup[ s, j, k ]  )
#     
#         }#K
#       }#J
#     }#S
#     
#     #defining data augmentation indicator prior
#     for( s in 1:S ){ #loop species
#       for( j in 1:J ){ #loop segments
#         for( k in 1:K ){ #loop years
#         #restrict data augmentation indicator to the corresponding BCR species      
#         w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, bcr.id[j] ]
#     
#         #indicator of whether species should be added to a particular segment
#         w[ s, j, k ] ~ dbern(  omega )
#         } #k
#       } #j
#     } #s
#     
#     
#     #priors
#     #data augmentation prior
#     omega ~ dbeta( 4, 4 )
#     
#     #priors for detection model 
#     
#     #Defining the intercept prior for alpha.0
#     int.p <- log(mean.p / (1 - mean.p))
#     mean.p ~ dbeta(4, 4)
#     
#     for(q in 1:Q){ #loop through predictors
#       alpha[ q ] ~ dnorm(0, 0.1)
#     }#Q
#     
#     #define intercepts for ecological model as mean probs:
#     #logit main intercept
#     int.psi <- log(mean.psi / (1 - mean.psi))
#     #mean occupancy probability
#     mean.psi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal
#     
#     #random year intercepts
#     for( k in 1:K ){
#       eps.psi[ k ] ~ dnorm( 0, prec.psi ) T( -6, 6 )
#     } #K loop
#     
#     #associated precision of radom year intercepts:
#     prec.psi <- 1 / (sigma.psi * sigma.psi) #variance of psi^2
#     sigma.psi ~ dt(0, 1, 4) T(0, ) #Variance of ps
#     
#     #random route effect
#     for( m in 1:M ){#loop routes
#       epsID.psi[ m ] ~ dnorm(0, precID.psi) T(-6, 6) #wide dist
#     }#M
#     
#     #associated precision for random route effect
#     precID.psi <- 1 / ( sigmaID.psi * sigmaID.psi )
#     sigmaID.psi ~ dt( 0, 1, 4 ) T(0, )
# 
#     for( s in 1:S ){ 
#       #for (k in 1:K){
#     #Pull out one species matrix * by non-surveyed sites and years
#     p.sp.mat[s] <- sum(p[ s, 1:J, 1:K ]  * JKsurv[ 1:J, 1:K ])
#     #p.sp.mat2 <- (means2$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
#     
#     #Sum across sites for each year
#     #p.sp.mat <- apply(p.sp.mat, 2, sum, na.rm = T)
#     #p.sp.mat2 <- apply(p.sp.mat2, 2, sum, na.rm = T)
#     
#     #divide by total number of sampled segments each year
#     p.sp[s] <-  p.sp.mat[ s ] / sum(surveyedJ[ 1:K ])
#     #p.mean[s] <- sum(p.sp[s, 1:K])
#     #p.sp2 <-  p.sp.mat2[ 1:K ] / surveyedJ[ 1:K ]
#     
#       #}#K
#     }#S
# 
#     } #model spec end
#     
#     ", fill = TRUE )
# 
# sink()
# 
# ################ end of model specification  #####################################
# modelname <- "pm1.txt"
# 
# #create initial values
# zst <- ydf
# zst[is.na( zst) ] <- 1
# inits <- function(){ list( z = zst ) }
# 
# #parameters monitored #only keep those relevant for model comparisons (with different variances)
# params <- c( 'int.psi' #intercept for occupancy model 
#              , 'eps.psi' #random year intercept
#              , 'epsID.psi' #random route intercept 
#              , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
#              , 'int.p' #intercept for detection submodel 
#              , 'alpha' #fixed coefficients for detection submodel 
#              , 'z' #only summary for estimated true occupancy
#              , 'p.sp' #only mean summary for detection probability
#              #, 'mean.a.div' #summary mean alpha diversity 
# )
# 
# #################################################################################
# ###### alternative  variable selection variances for occupancy model ############
# str( win.data <- list( ydf = ydf, #observed occupancy 
#                        J = J, K = K, S = S, Q =  4, M = M,
#                        #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
#                        JKsurv = JKsurv, #indicator of segments sampled (1) and not (0)
#                        bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
#                        rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
#                        bcr.occ = bcr.occ, #bcr indicator for each species
#                        Mass.scaled = spp.occ$Mass.scaled, #body mass
#                        TOD.ma = TOD.ma, #time of day
#                        Ord.ma = Ord.ma, #day of year
#                        Obs.ma = Obs.ma, #first observer year
#                        surveyedJ = surveyedJ
# ) )                
# #library( jagsUI )
# #call JAGS and summarize posteriors:
# ptm <- proc.time()
# #fm1 <- jags( win.data, inits = inits, params, modelname,
# #             n.chains = nc, n.thin = nt, n.iter = ni,
# #             n.burnin = nb, parallel = TRUE)
# 
# #auto update the model
# upm1 <- autojags( win.data, inits = inits, params, modelname,
#           n.chains = 3, n.thin = 10, n.burnin = 1000,
#           iter.increment = 1000, max.iter = 100000,
#           Rhat.limit = 1.15, save.all.iter=FALSE, parallel = TRUE )
# 
# fm1.time <- proc.time() - ptm
# 

##########################################
##### Multi-species occ model single visit ######
##### For species within a single phylogenetic group #######
#### Detection related to siteXyear fixed predictors: #######
##### including time of day, first-year observer effect #######
##### and species mass as well as random species intercept ##
## to borrow information across phylogenetic group ###
### Ecological model includes data augmentation parameter ###
### for species in BCR  and         #####
#####Random intercepts for year and site #######
##########################################

####Model 1 but with species random effects in detection model###

########Model specification###############
sink( "pm2.txt" )
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
    alpha[ q ] ~ dnorm(0, 0.1)
    }#Q
    
    #Random species effect delta[s]
    for(s in 1:S){ #loop through species
    delta[s] ~ dnorm( 0, prec.delta ) T( -6, 6 )
    } #S loop
    
    #Precision for random species effects
    prec.delta <- 1 / ( sigma.delta * sigma.delta )
    sigma.delta ~ dt( 0, 1, 4 ) T( 0, ) 
    
    
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
    
    for( s in 1:S ){ 
      #for (k in 1:K){
    #Pull out one species matrix * by non-surveyed sites and years
    p.sp.mat[s] <- sum(p[ s, 1:J, 1:K ]  * JKsurv[ 1:J, 1:K ])
    #p.sp.mat2 <- (means2$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
    
    #Sum across sites for each year
    #p.sp.mat <- apply(p.sp.mat, 2, sum, na.rm = T)
    #p.sp.mat2 <- apply(p.sp.mat2, 2, sum, na.rm = T)
    
    #divide by total number of sampled segments each year
    p.sp[s] <-  p.sp.mat[ s ] / sum(surveyedJ[ 1:K ])
    #p.mean[s] <- sum(p.sp[s, 1:K])
    #p.sp2 <-  p.sp.mat2[ 1:K ] / surveyedJ[ 1:K ]
    
      #}#K
    }#S
    
    } #model spec end
    
    ", fill = TRUE )

sink()

################ end of model specification  #####################################
modelname <- "pm2.txt" #changed model name to phylomodel2 
#create initial values
zst <- ydf
zst[is.na( zst) ] <- 1
inits <- function(){ list( z = zst ) }

#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'int.psi' #intercept for occupancy model 
             , 'eps.psi' #random year intercept in occupancy model
             , 'epsID.psi' #random route intercept in occupancy model
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
             , 'delta', 'sigma.delta' #random species intercept in detection model
             , 'z'
             #, 'p'
             , 'p.sp'
             #### need to add summary stats to the model ####
             # , 'z.prime' #only summary for estimated true occupancy
             # , 'mean.p.sp' #only mean summary for detection probability
             # , 'mean.a.div' #summary mean alpha diversity 
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       J = J, K = K, S = S, Q =  4, M = M,
                       #J = 50, K = 5, S = 233, G = 23, Q =  4, M = 84,
                       JKsurv = JKsurv, #indicator of segments sampled (1) and not (0)
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma, #first observer year
                       surveyedJ = surveyedJ
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
ptm <- proc.time()
# fm2 <- jags( win.data, inits = inits, params, modelname,
#              n.chains = nc, n.thin = nt, n.iter = ni,
#              n.burnin = nb, parallel = TRUE)


#auto update the model
upm1 <- autojags( win.data, inits = inits, params, modelname,
          n.chains = 3, n.thin = 10, n.burnin = 1000,
          iter.increment = 1000, max.iter = 50000,
          Rhat.limit = 1.15, save.all.iter=FALSE, parallel = TRUE )

fm1.time <- proc.time() - ptm


