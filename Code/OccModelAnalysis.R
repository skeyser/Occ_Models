##########################################################
#########Single Visit Multispecies Occupancy Model########
######Script has code for JAGS specified model############
###############Created By: Spencer R Keyser###############
##########################################################

######Load Packages######
rm(list = ls())

library(pacman)
#pacman::p_load(tidyverse, ggplot2, ggfortify, JagsUI)
install.packages("jagsUI")
library(jagsUI)
#########Package Loading Complete########
#########################################

###Load in the Data4OccModels Workspace###
##########################################

load('Data4OccModels_6_3.RData')


#########################################################################################
###########################################################################
####################### define MCMC settings ##############################

ni <- 2; nt <- 0; nb <- 0; nc <- 3 #iterations, thinning, burnin, chains

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

#####Playing with Beta distribution#######
p = seq(0, 1, length = 100)
plot(p, dbeta(p, 100, 100), ylab = "density", type = "l", col = 4)
lines(p, dbeta(p, 1, 10), type = "l", col = 1)
lines(p, dbeta(p, 1, 1), type = "l", col = 2)
lines(p, dbeta(p, 4, 4), type = "l", col = 3)
lines(p, dbeta(p, 100, 100), type = "l", col = 4)
lines(p, dbeta(p, 10, 1), type = "l", col = 5)
lines(p, dbeta(p, 0, 0), type = "l", col = 6)

###########################################
#Hello

########Model specification###############
sink( "om1.txt" )
cat("
    model{
    
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

    #defining data augmentation indicator prior
    for( s in 1:S ){ #loop species
       for( j in 1:J ){ #loop segments
      #   for( k in 1:K ){ #loop years
          omega[ s, j ] ~ dbeta(4,4)
          #restrict data augmentation indicator to the corresponding BCR species
          w[ s, j ] ~ dbern(  omega[ s, j ] )          
          w.bcr[ s, j ] <- w[ s, j ] * bcr.occ[ s, bcr.id[j] ]
      #   } #k
       } #j
    } #s

    #ecological model
    #for occupancy
    for( j in 1:J ){ #loop segments
       for( k in 1:K ){ #loop years
        #relate occupancy probability to random intercepts for route and year:
          logit( psi[ j, k ] ) <- int.psi + epsID.psi[ rteno.id[ j ] ] + 
                                  eps.psi[ k ]
        for( s in 1:S ){ #loop species
          #modeling true occupancy
          muz[ s, j, k ] <- psi[ j, k ] * w.bcr[ s, j ]
          z[ s, j, k ] ~ dbern( muz[ s, j, k ] )            
        }#S
      }#K
    }#J

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
             , 'w.bcr' #indicator variable of which species are added as augmented set
             , 'eps.psi' #random year intercept
             , 'epsID.psi' #random route intercept 
             , 'sigma.psi', 'sigmaID.psi' #std dev for random intercepts
             , 'z' #estimated true occupancy
             , 'p' #detection probability
             , 'int.p' #intercept for detection submodel 
             , 'alpha' #fixed coefficients for detection submodel 
 #            , 'delta' #random intercept for phylogenetic group
#             , 'a.div' #alpha diversity
)

#################################################################################
###### alternative  variable selection variances for occupancy model ############
str( win.data <- list( ydf = ydf, #observed occupancy 
                       J = J, K = K, S = S, G = G, Q =  4, M = M,
                       #J = 10, K = 5, S = S, G = 23, Q =  4, M = M,
                       bcr.id = jdf$bcr.id, #indicator of what BCR the segment belongs to
                       rteno.id = jdf$rteno.id, #indicator of what route the segment belongs to
                       bcr.occ = bcr.occ, #bcr indicator for each species
                       Mass.scaled = spp.occ$Mass.scaled, #body mass
#                       Phylo.V1.code = spp.occ$Phylo.V1.code, #phylo grouping
                       TOD.ma = TOD.ma, #time of day
                       Ord.ma = Ord.ma, #day of year
                       Obs.ma = Obs.ma #first observer year
                       
) )                
#library( jagsUI )
#call JAGS and summarize posteriors:
fm1 <- jags( win.data, inits = inits, params, modelname, #
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE ) 

#auto update the model
upm1 <- autojags( win.data, inits = inits, params, modelname, 
          n.chains = nc, n.thin = nt, n.burnin = 0,
          iter.increment = 100, max.iter = 1000,
          Rhat.limit = 1.1, save.all.iter=FALSE, parallel=TRUE )

################ end of script ##########################################