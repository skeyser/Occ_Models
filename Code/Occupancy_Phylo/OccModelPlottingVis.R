###################################################################
################Summary Stats and plotting#########################
#############for Multi species occpancy models#####################
###############Created by : Spencer R Keyser#######################
###################################################################

#rm("upm1")

#Load in packages
library(pacman)
pacman::p_load("MCMCvis", "LaplacesDemon", "truncnorm")

#Package loading complete 

J <- 274
K <- 38
S <- max(spp.occ$spp.id)

###################################
#This bit is in the model body now#
###################################

# #Mean detection across species 
# mean.p.sp <- data.frame(Species = 1:length(unique(spp.occ$Species)), "Mean p No RE" = NA,
#                         "Mean p RE" = NA)
# 
# mean.p.sp$Species <- spp.occ$Species
# 
# #Loop throug all species
# for( s in 1:length(unique(spp.occ$spp.id)) ){ 
#   #Pull out one species matrix * by non-surveyed sites and years
#   p.sp.mat <- (means.output$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
#   p.sp.mat2 <- (means2$p[ s, 1:J, 1:K ]  * JKmat[ 1:J, 1:K ])
#   
#   #Sum across sites for each year
#   p.sp.mat <- apply(p.sp.mat, 2, sum, na.rm = T)
#   p.sp.mat2 <- apply(p.sp.mat2, 2, sum, na.rm = T)
#   
#   #divide by total number of sampled segments each year
#   p.sp <-  p.sp.mat[ 1:K ] / surveyedJ[ 1:K ]
#   p.sp2 <-  p.sp.mat2[ 1:K ] / surveyedJ[ 1:K ]
#   
#   #Put values into empty DF
#   mean.p.sp [s, 2] <- mean(p.sp)
#   mean.p.sp [s, 3] <- mean(p.sp2)
# }#S

##################################

#Extract subsets of the full model for saving the workspace 
means.output <- upm1$mean
Rhats <- upm1$Rhat
model.sum <- upm1$summary


#View All non-converged parameters 
bad.params <- function(x) {
  rhats <- unlist(Rhats)
  buggers <- rhats[rhats > 1.15]
  buggers[complete.cases(buggers)]
}

non.converge <- bad.params(Rhats)

###summaries of model parameters ########
###here is an example if we actually want to estimate & plot #
# relationships between alpha and detection prob #
# we could also save mean and CIs for the random effects as a #
# separate dataframe exported elsewhere # 

#define range for important covariates in models
#intercept for covariate matrix 
int <- rep( 1, 100 )
#cbind( int, sclegl)
#annual eagle abundance
sclegl <- seq( min( XK[ , 'EglAbund'], na.rm = TRUE ), max( XK[ , 'EglAbund'], na.rm = TRUE ), 
               length.out = 100 )
egl <- seq( min( Kcovs[ , 'EglAbund' ], na.rm = TRUE ), 
            max( Kcovs[ , 'EglAbund' ], na.rm = TRUE ), length.out = 100 )
#estimate mean relationship with predictor while keeping all other covs constant at their mean
#multiply coefficient matrix (for all iterations) against transposed covariate matrix using matrix algebra
phi.egl <- cbind( mr$sims.list$int.phi, mr$sims.list$beta.phi[,2] ) %*% 
  t(cbind( int, sclegl ) ) 
#calculate its mean and get inverse logit
phi.eglm <- expit( apply( phi.egl, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
phi.eglCI <- expit( apply( phi.egl, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) ) )

### you can then save output into dataframe and save it so #
# we can then plot it elsewhere#
#combine predicted estimates into a dataframe
predrships <- data.frame( sclegl, egl,  #standardise and not standardised covariates
                          phi.eglm, t( phi.eglCI ) #predicted response
)
# to plot:
predp <- ggplot( data = predrships ) + theme_classic() +
  theme( legend.position = "none", 
         text = element_text( size = 18 ), 
         axis.line = element_line( size = 1.3 ) ) + ylim( 0.0,1.0 ) 

ap <- predp +  xlab( "Bald eagle abundance" ) +
  geom_line( aes( x = egl, y = phi.eglm ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = egl, ymin = X2.5. , ymax = X97.5. ) ) + #, 
  ylab( "Probability of persistence" )            
##### end #########
#################################################################
#Pull out just the z matrix
z.prime <- means.output$z
#z.prime2 <-means2$z

#Loop through each species, site, and year to set 0s and 1s 
#manually using cut-off values

for (s in 1:S){
  for (j in 1:J){
    for (k in 1:K){
      #Generate multiple "cut-offs" for sensitivity analyses
      z.prime.5[s, j, k] <- ifelse(means1$z[s, j, k] > 0.5, 1, 0)
      z.prime.65[s, j, k] <- ifelse(means1$z[s, j, k] > 0.65, 1, 0)
      z.prime.75[s, j, k] <- ifelse(means1$z[s, j, k] > 0.75, 1, 0)
      #z.prime2[s, j, k] <- ifelse(means2$z[s, j, k] > .5, 1, 0) 
    }#K
  }#J
}#S

#Multiply each species matrix by the NA df for species totals
for( i in 1:S ){
  z.prime.5[ i, , ] <- z.prime.5[ i, , ] * JKmat #JKsurv
  z.prime.65[ i, , ] <- z.prime.65[ i, , ] * JKmat#JKsurv
  z.prime.75[ i, , ] <- z.prime.75[ i, , ] * JKmat #JKsurv
}

#Check differences
z.prime[1,,]
z.prime.5[1,,]
z.prime.65[1,,]
z.prime.75[1,,]

#Quick function to calculate total observations of each species
total <- function(x, y){
  df <- x[y,,] 
  #did I just break your function by using JKmat instead?
  yr.sums <- colSums(df, 1:length(nrow), na.rm = TRUE)
  sum(yr.sums, na.rm = TRUE)
}

#Create empty DF for storing values 
#total.observed <- data.frame("Observed" = 1:max(spp.occ$spp.id), "Estimated RE" = NA, 
#                             "Estimated No Sp RE" = NA)

total.observed <- data.frame("Observed" = 1:max(spp.occ$spp.id), "Estimated 0.5" = NA, 
                             "Estimated 0.65" = NA, "Estimated 0.75" = NA)

for (i in 1:S){#max(spp.occ$spp.id)){
  total.observed[i, 1] <- total(ydf, i)
  total.observed[i, 2] <- total(z.prime.5, i)
  total.observed[i, 3] <- total(z.prime.65, i)
  total.observed[i, 4] <- total(z.prime.75, i)
}

#Find yearly alpha diversity per subgroup
#Generate empty matrices
group.a.div <- matrix(NA, J, K )#274, 38)
group.a.div.5 <- matrix(NA, J, K )# 274, 38)
group.a.div.65 <- matrix(NA, J, K )#274, 38)
group.a.div.75 <- matrix(NA, J, K )#274, 38)

for (j in 1:J){
  for ( k in  1:K){
    #zeros those z for unsampled segments
    group.a.div[ j, k ] <- sum( ydf[1:S, j, k ], na.rm = TRUE )
    ## we timed it by JKmat above so don't need to do it again right?
    group.a.div.5[ j, k ] <- sum( z.prime.5[ 1:S, j, k ], na.rm = TRUE)# * JKsurv[ j, k ] )
    group.a.div.65[ j, k ] <- sum( z.prime.65[ 1:S, j, k ], na.rm = TRUE)# * JKsurv[ j, k ] )
    group.a.div.75[ j, k ] <- sum( z.prime.75[ 1:S, j, k ], na.rm = TRUE)# * JKsurv[ j, k ] )
    #keep only average estimate
    #mean.a.div[ j, k ] <- mean(a.div[ j, k ])
  }#K
}#J

alpha.div <- list(group.a.div, group.a.div.5, group.a.div.65,
                  group.a.div.75)



###########################################################################
####################################Plotting###############################
###########################################################################

##Traceplots for relevant parameters

#Priors used for plotting
#Alpha params
PR.norm <- rnorm(1000, 0, 3.126)

#Random effects precision
sigma.dt <- rhalft(1000, scale = 1, nu = 4)
prec.df <- 1 / (sigma.dt * sigma.dt)

#Conditional normal for random effects
PR.norm.cond <- rtruncnorm(1000, a = -6, b = 6, 0, sigma.dt)

#Random intercepts prior 
PR.beta <- rbeta(1000, 4, 4)


#Plotting 
#can we extract part of the workspace name so that we can use that to #
# label the files saved...that way it can be done automatically instead of manually?
#
#path where i get the name of the workspace I am running:
extpath <- getwd()
# extract name of file to use: 
workname <-  dir( extpath, pattern = '\\.RData', full.names = FALSE )
workname <- gsub(".RData", "", workname )
#does that work? we may need to have separate folders for each phylogroup so that 
# or there may be a way easier way of doing it?

MCMCtrace(upm1, params = 'alpha', priors = PR.norm, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'sigma.delta', priors = sigma.dt, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'delta', priors = PR.norm.cond, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'int.p', priors = PR.beta, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'int.psi', priors = PR.beta, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'epsID.psi', priors = PR.norm.cond, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'eps.psi', priors = PR.norm.cond, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'sigmaID.psi', priors = sigma.dt, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))

MCMCtrace(upm1, params = 'sigma.psi', priors = sigma.dt, ind = T, Rhat = T, 
          n.eff = T, pdf = T, file = here::here("Figures and Tables/TraceSwallows"))



#Parameter Plots 
MCMCplot(upm1, params = c("alpha"), rank = T, main = "Swallows Phylo Group Alphas", labels = c("Observer Effect", "Time of Day",
                                                                                               "Day of Year", "Body Mass"))

MCMCplot(upm1, params = c("int.psi", "int.p"), rank = T, main = "Phylo Group Intercepts", labels = c("Intercept Occupancy", 
                                                                                                              "Intercept Detection")) 
MCMCplot(upm1, params = c("epsID.psi", "sigmaID.psi"), main = "Phylo Group Site Effects")
                                                                                                                
MCMCplot(upm1, params = c("eps.psi", "sigma.psi"), main = "Phylo Group Year Effects")

MCMCplot(upm1, params = c("delta", "sigma.delta"), main = "Species RE", labels = c("BASW", "CASW", "CLSW", 
                                                                                            "PUMA", "BASW", "NRWS", 
                                                                                            "Sigma Delta"))
rm(list = setdiff(ls(), c("means.output", "Rhats", "model.sum",
                          "z.prime", "z.prime.5", "z.prime.65",
                          "z.prime.75", "total.observed", 
                          "alpha.div", "non.converge")))

##Save workspace
what.dir <- "R Workspace/Output"
phylo.group <- as.character(unique(spp.occ$Phylo.V1))
phylo.group <- gsub("/", "_", phylo.group)
phylo.group <- paste0(phylo.group, "", ".RData")

#save.image(file = here::here(paste0(what.dir, "/", phylo.group)))
         