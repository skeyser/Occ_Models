##########################################################
#########Single Visit Multispecies Occupancy Model########
######Script has code for JAGS specified model############
###############Created By: Spencer R Keyser###############
##########################################################

######Load Packages######
rm(list = ls())

library(pacman)
pacman::p_load(tidyverse, ggplot2, ggfortify, R2jags)


#########Package Loading Complete########
#########################################

###Load in the Data4OccModels Workspace###
##########################################

load('workspace.RData')

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
Hello

########Model specification###############
sink("om1.txt")
cat("
    model{
    
    #prior
    #for occupancy model:
    psi1 ~ dunif(0, 1)
    #define intercepts for ecological model as mean probs:
    #Log-link function below
    int.phi <- log(mean.phi / (1 - mean.phi))
    mean.phi ~ dbeta(4,4) #mean occupancy, beta(4,4) similar to normal

    #random species intercepts
    for (s in 1:S){
    eps.phi[s] ~ dnorm(0, prec.phi) T(-10,10)
    } #S

    #associated precision of radomm species intercepts:
    prec.phi <- 1 / (sigma.phi * sigma.phi) #variance of phi^2
    sigma.phi ~ dt(0, 1, 4) T(0, ) #Variance of phi
    
    #random site effect
    for (m in 1:M){#loop sites
    epsID.phi[m] ~ dnorm(0, precID.phi) T(-10, 10) #wide dist
    }#M

    #associated precision for site random effect
    precID.phi <- 1 / (sigmaID.phi * sigmaID.phi)
    sigmaID.phi ~ dt(0, 1, 4) T(0, )

####Not doing random slopes?####
    #priors for beta predictors: 
    for(q in 1:Q){ #loop over number of predictors
    }


    #detection model (ObsChange, Date, Time)
    

    #ecological model
    #for occupancy
    for(s in 1:S){ #loop species
      for(j in 1:J){ #loop sites
      #estimated psi[,,1]:
        ydf[s, j, 1] ~ dbern(psi[s, j, 1])
        psi[s, j, 1] <- psi1
        
        for(k in 1:K){ #loop years 1980:2017
          #Detection corrected, p != 1
          ydf[s, j, k] ~ dbern(psi[s, j, k])
          #probability of occupancy starts at year 2
          psi[s, j, k] <- ydf[s, j, k-1] * phi[s, j, k-1] *
          p[s, j, ] #p is detection probability
}
        
}
}
    
    }")