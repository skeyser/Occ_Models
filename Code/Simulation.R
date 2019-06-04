##############################################################################
###                  Code created by Jen Cruz                         ########
#### This script contains methods for simulating data based on model  #######
#### results and reruning the model ( to check that it is valid)  ######
###   Model is a reduced model estimating nest persistence, colonization #####
####        and success in DynOcc-Suc_ReducedModels.R script            ######
##############################################################################

########### clean workspace and load required packages ####################
###########################################################################

#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
#set working dir() 
gc() #releases memory
#gcinfo( FALSE ) #views memory allocation

####### load relevant packages ###
#library( tidyr ) #dataframe combining.
library( dplyr ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #to use spread and other df functions
#library( ggplot2 ) #fancy plots.
#library( lubridate ) #easy date adjustments and calculations.
library( jagsUI ) #to run RJAGS

########## end of package loading ###########
getwd() #check working directory
#functions
##############
`expit` <-
  function(x){
    exp(x)/(1+exp(x))
  }

#plogis is the R equivalent to expit

`logit` <-
  function(x){
    log(x/(1-x))
  }
#### end of functions #####
#######################    import relevant data   

# #start with final model results:
load( "D:/RProjects/Osprey/ModelResults.RData" )
# Define model we want to extract output from 
mr <- fm1
#view model parameters
mr$mean$int.p
mr$mean$alpha
summary( mr )

##################################################################################################
######              Simulating data based on model results                ########
##################################################################################################

######## simulating the ecological relationships ###########################
# We use actual predictor values for each site and season but use the mean coefs that were #
#derived from the model to simulate data #

# Function to simulate nest occupancy and success data:
SimData <- function(){
  #Function that simulates nest occupancy and success observations based on model output
  # from rem1 which is a dynamic occupancy/success model for osprey
  
  #create temporary objects to hold intermediate output inside the function
  
  psi <- matrix( NA, nrow = J, ncol = K )
  z <- p <- array( NA, dim = c(S, J, K) )
  # prob of colonization [i, k ]

  for( i in 1:J ){ # loop over sites
    for( k in 2:K ){ #loop over 2:K years
    
      psi[ j, k ] <-  expit( mr$mean$int.psi + int.psi + epsID.psi[ jdf$rteno.id[ j ] ] + 
                           eps.psi[ k ] )
      for( s in 1:S ){ # loop over species
        #define data augmentation prob of inclusion
        w[ s, j, k ] <- rbinom( n = 1, size = 1, prob = 0.5 ) 
        w.bcr[ s, j, k ] <- w[ s, j, k ] * bcr.occ[ s, jdf$bcr.id[ j ] ]
        # occupancy state
        z[ s, j, k ] <- rbinom( n = 1, size = 1, prob = (psi[ j, k ] * w.bcr[s,j,k] ) )
    }# close k loop
    #simulate observations of nest occupancy for each site each season [i,k]
    ysim[ s,j,k ] <<- rbinom( n = K, size = 1, prob = p[ s, j, k ] * z[ s, j, k ] )
    for( k in 1:K ){ #loop over years
      # prob of nest success[ i, k ]
      zeta[ i, k ] <- expit( mr$mean$int.zeta + mr$mean$beta.zeta[2] * XK[ k, 'EglAbund' ] +
      #zeta[ i, k ] <- expit( mr$mean$int.zeta + mr$mean$beta.zeta * XK[ k, 'EglAbund' ] +
                             mr$mean$eps.zeta[ k ] + mr$mean$epsID.zeta[ i ] )
      #simulated observations of nests success only for sites that were occupied [i,k]
      sucsim[ i, k ] <<- rbinom( n = 1, size = 1, prob = zeta[ i, k ] )  * z[ i, k ]
      
    } #close k loop
  } #close i loop
  print( 'phi')
  print( phi[ 1:5, 1:6 ] )
  print( "psi" )
  print( psi[ 1:5, 1:6 ] )
  print( "zeta" )
  print( zeta[ 1:5, 1:6] )
  for( i in 1:M ){ #loop over sites
    
  }
  print( 'z' )
  print( z[1:5,1:10] )
  print( "ysim" )
  print( ysim[1:5, 1:10] )
  print( "sucsim" )
  print( sucsim[1:5, 1:10] )

} #close function

# Create matrices to store function output
ysim <- sucsim <- matrix( NA, nrow = M, ncol = K )
# Run function
SimData()
# View output
table( ysim ); table( sucsim )
colSums( ysim )
colSums( sucsim )
# Rename matrices so that you can run the function again:
ysim1 <- ysim
sucsim1 <- sucsim 
# # Run function
# SimData()
# # View output
# table( ysim ); table( sucsim )
# # Rename matrices so that you can run the function again:
# ysim2 <- ysim
# sucsim2 <- sucsim 


##################### rerun model ####################
#modelname <- "redm2.txt"
modelname <- "am1.txt"

#create initial values
zst <- ysim1
zst[] <- 1
head( zst )
inits <- function(){ list( z = zst ) }

#parameters monitored
params <- c( 'int.phi' #intercept for persistence model 
             , 'int.gam' #intercept for colonization model
             , 'int.zeta' #intercept for nest success model
             , 'beta.phi' #occupancy model coefficients included in the model (where indicator delta = 1 )
             , 'beta.gam' #coefficients for colonization model
             , 'beta.zeta' #coefs for nest success model
             , 'n_occ' #yearly occupancy 1986-2012
             , 'repsucc' #reproductive success: rate of nests that produced >0 chicks
             , 'psi.k' #yearly nest occupancy
             , 'phi.k' #mean yearly nest persistence
             , 'gam.k' #mean yearly nest colonization
             , 'zeta.k' #mean yearly prob of nest success
             , 'p' #yearly detection for occupancy model
             , 'int.p' #intercept for detection submodel of occupancy model
             , 'alpha.p' #coefficients for detection submodel 
 )
#################################################################################
###### alternative  variable selection variances for occupancy model ############
# str( win.data <- list( y_obs = ysim2, M = M, K = K , surv_ind = surv_ind, #indicator of whether nest was surveyed 
#                        suc_obs = sucsim2, 
#                        ForagArea = XM$ForagArea , inVNP = XM$inVNP, #indicator of whether nest was inside VNP 
#                        modXK = as.matrix( XK[ , c('EglAbund', #'CPUEyr',
#                                                   'FirstDay') ] ),
#                        #DeltaO = DeltaOcc_mat,
#                        c.var = 50, tau.var = 0.05 #variance values for beta priors
# ) ) 
str( win.data <- list( y_obs = ysim, M = M, K = K , surv_ind = surv_ind, #indicator of whether nest was surveyed 
                       suc_obs = sucsim, #failprob = c( 0, 0, 1 ),
                       #r_obs = r_obs, f_obs = f_obs, hp_obs = hp_obs, #binomial obs of successful, failed and highly prod nests
                       ForagArea = XM$ForagArea , inVNP = XM$inVNP, #indicator of whether nest was inside VNP 
                       modXK = as.matrix( XK[ , c('EglAbund', 'Fish',
                                                  'RainyJulian', 'Mean_minT', 'Mean_totRain',
                                                  'FirstDay') ] ), XKno = 5 ,#XKno is no. of K covs in ecological model only
                       MinDistO = MinDistOcc_mat, MinDistS = MinDistSucc_mat, 
                       DeltaO = DeltaOcc_mat, DeltaS = DeltaSucc_mat, 
                       UncontO = UncontOcc_mat, UncontS = UncontSucc_mat, #XKMno = 6,
                       CPUE = CPUE_mat, CPUEXArea = CPUEXArea_mat,
                       PondAge = XM$PondAge, #centered categorical of whether nest was in pond or not
                       Xtot = 15, #XKno + XKMno + 2: of which 1 is for 1 Mlevel cov & 1 for interaction bw foragarea and eglabund
                       c.var = 50, tau.var = 0.05 #variance values for beta priors
) )
#library( jagsUI )
#call JAGS and summarize posteriors:
sm1 <- jags( win.data, inits = inits, params, modelname, #
             n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb, parallel = TRUE ) 
# # using ysim2 and sucsim2 data:
# sm2 <- jags( win.data, inits = inits, params, modelname, #
#              n.chains = nc, n.thin = nt, n.iter = ni, 
#              n.burnin = nb, parallel = TRUE ) 

#####################    View results ######################################
#assign general name to simulation results
sm <- sm1

par( mfrow = c( 3, 3 ), ask = F , mar = c(3,2,2,2), cex = 1.3, tcl = -0.2,
     lwd = 2 )
plot( density( sm$sims.list$int.phi ), xlab ="", ylab = "", 
      main = expression( beta[phi][0] ), lty = 2  )
#lines( density( mr$sims.list$int.phi ), lwd = 2 )
abline( v = mr$mean$int.phi )#, lwd = 2 )
plot( density( sm$sims.list$int.gam ), xlab ="", ylab = "",
      main = expression( beta[gamma][0] ), lty = 2  )
#lines( density( mr$sims.list$int.gam ), lwd = 2 )
abline( v = mr$mean$int.gam )
plot( density( sm$sims.list$int.zeta ), xlab ="", ylab = "", 
      main = expression( beta[zeta][0] ), lty = 2  )
#lines( density( mr$sims.list$int.zeta ), lwd = 2 )
abline( v = mr$mean$int.zeta )
plot( density( sm$sims.list$beta.phi[,2] ), xlab ="", ylab = "",
        main = expression( beta[phi][2] ), lty = 2  )
#lines( density( mr$sims.list$beta.phi[,2] ), lwd = 2 )
abline( v = mr$mean$beta.phi[2] )

plot( density( sm$sims.list$beta.zeta[,2] ), xlab ="", ylab = "", 
      main = expression( beta[zeta][2] ), lty = 2  )
#lines( density( mr$sims.list$beta.zeta[,2] ), lwd = 2)
abline( v = mr$mean$beta.zeta[2] )
plot( 1,1, xaxt = "n", yaxt = "n", col = "white" , xlab ="", ylab = "")
plot( density( sm$sims.list$int.p ),xlab ="", ylab = "", 
      main = expression( alpha[0] ), lty = 2  )
#lines( density( mr$sims.list$int.p ), lwd = 2 )
abline( v = mr$mean$int.p, lwd = 2 )
plot( density( sm$sims.list$alpha.p[,1] ), xlab ="", ylab = "", 
      main = expression( alpha[1] ), lty = 2  )
#lines( density( mr$sims.list$alpha.p[,1] ), lwd = 2)
abline( v = mr$mean$alpha.p[2], lwd = 2 )
plot( density( sm$sims.list$alpha.p[,2] ), xlab ="", ylab = "", 
      main = expression( alpha[2] ), lty = 2  )
#lines( density( mr$sims.list$alpha.p[,1] ), lwd = 2)
abline( v = mr$mean$alpha.p[2], lwd = 2 )


##########  timeseries plots ##########
yrs <- length( yrrange )-1
rowid <- c( 0, yrs, ( yrs * 2 ), ( yrs * 3 ) )
rndlength <- rep( yrrange[1:yrs], 3 )
idcolname <- 'Year'
rndeffs <- c( "phi.k", "gam.k", "zeta.k" )
rnddfnames <- c( 'var', 'real', 'Estimate', 'lower', 'higher', "RE" )
rnddf <- data.frame( matrix( NA, ncol = length( rnddfnames ), 
                             nrow = length( rndlength ) ) )
colnames( rnddf ) <- rnddfnames
rnddf$var <- factor( rnddf$var, levels = rndeffs ) 
rnddf[ ,6 ] <- rndlength

for( c in 1:length( rndeffs ) ){
  #for( k in 1:M ){ 
  for( k in 1:yrs ){  
    rowval <- rowid[ c ] + k
    rnddf$var[ rowval ] <- rndeffs[ c ]
    varname <- paste0( "sm$sims.list$", rndeffs[ c ], "[,", k, "]" )
    realvar <- paste0( "mr$sims.list$", rndeffs[ c ], "[,", k, "]" )
    rnddf[ rowval, 2 ] <- mean( eval( parse( text = realvar ) ), na.rm = TRUE )
    rnddf[ rowval, 3 ] <- mean( eval( parse( text = varname ) ), na.rm = TRUE )
    rnddf[ rowval, 4:5 ] <- as.numeric( quantile( eval( parse( text = varname ) ),
                                                  c( 0.025, 0.975 ), na.rm = TRUE ) )
    #print( c( rowval, varname) )#, moddf[ rowval, 3], moddf[ rowval, 4:5] ) )
  }}

head(rnddf)

#### create ggplot elements:
# varlabels <- c( expression( paste( phi[k] ) ), expression( paste( gamma[k] ) ),
#                  expression( paste( zeta[k] ) ) )
# 
# rndgpnames <- c( 'phi.k' = label_bquote( phi[k] ), 
#                  'gam.k' = label_bquote( gamma[k] ),
#                  'zeta.k' = label_bquote( zeta[k] )
# )

ggplot( rnddf, aes( x = RE, y = Estimate, group = var ) ) +
  theme_bw() + 
  theme( legend.position = "none", 
         panel.border = element_rect( size = 1.2 ),
         axis.line = element_line( size = 1.3 ), 
         axis.text = element_text( size = 15 ),
         axis.title = element_text( size = 18, face = "bold" ),
         plot.title = element_text( size = 18, face = "bold" ),
         strip.text = element_blank(),#text( size = 20, face = "bold" ) ,
         strip.background = element_blank()) + 
  xlab( idcolname ) + 
  facet_grid( var ~. , scales = "free_y" ) +#, labeller = as_labeller( rndgpnames ) ) + # ) +
  geom_point( size = 3 ) +
  geom_line( aes( y = real ), size = 1.4 ) +
  geom_ribbon( aes( ymin = lower, ymax = higher, group = var ),
                 alpha = 0.4 ) +
  scale_x_continuous( breaks = round( seq( min(rndlength), max(rndlength), by = 6 ), 0 ) )

# ### using base:
# par( mfrow = c( 4, 1 ), ask = F , mar = c( 2, 5, 1, 1 ), cex = 1.3, tcl = -0.2,
#      lwd = 2 )
# plot( yrrange[1:K-1], mr$mean$phi.k, ylim = c(0,1), 
#       ylab = expression( phi[k] ) )
# mr$summary
# lines( yrrange[1:K-1], sm$mean$phi.k )
# plot( yrrange[1:K-1], mr$mean$gam.k, ylim = c(0,1), 
#       ylab = expression( gamma[k] ) )
# lines( yrrange[1:K-1], sm$mean$gam.k )
# plot( yrrange, mr$mean$zeta.k, ylim = c(0,1), 
#       ylab = expression( zeta[k] ) )
# lines( yrrange, sm$mean$zeta.k )
# plot( yrrange, mr$mean$p, ylim = c(0,1), 
#       ylab = expression( p[k] ), xlab = "Year" )
# lines( yrrange, sm$mean$p )
# 
# par( mfrow = c( 2, 1 ), ask = F , mar = c(3,4,2,2), cex = 1.3, tcl = -0.2,
#      lwd = 2 )
# plot( x = yrrange, y = mr$mean$n_occ, ylab = "Osprey abundance" )#, ylim = c(5,30)  )
# lines( yrrange, sm$mean$n_occ, lty = 2 )
# whiskerplot( sm, parameters = c( 'n_occ' ) )
# lines( 1:K, mr$mean$n_occ, lwd = 2, lty = 1 )
# 
# plot( colSums(y_obs, na.rm = TRUE ), mr$mean$n_occ, lwd = 2 )
# lines( yrrange, colSums(y_obs, na.rm = TRUE ), lty = 2 )

######################    END OF SCRIPT             ###########################################
##################################################################################################