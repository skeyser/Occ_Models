# "metric". [gow.method]


# Outputs are provided in a matrix and include for each site (row) the site's
# ID as given in row names of the site-by-species matrix, the total number of 
# species occurring at that site, and functional richness.

MSTrichness<-function(abundances, traits, 
                      gow.w=NULL, gow.ord="metric"){
  
  # load required packages
  library(vegan) 	# for function spantree
  library(FD) 	# for function gowdis
  
  # check coherence of species in 'traits' and 'abundances'
  if(sum(is.element(dimnames(abundances)[[2]], dimnames(traits)[[1]]))!= dim(abundances)[2]){
    stop(" ERROR : not all species in the site-by-species matrix are represented in the traits matrix ")
  }	
  
  # check for absence of NA in 'traits'
  if (sum(is.na(traits))!=0) print(" WARNING : NAs in 'traits' matrix ")
  
  # if any traits are entirely NA, remove
  for(t in 1:ncol(traits)){
    if(sum(is.na(traits[,t]))==nrow(traits)){
      print(paste(" WARNING: Trait", dimnames(traits)[[2]][t], 
                  "is NA throughout, so will automatically be removed from analysis", sep=" "))
      traits<-traits[,-t]
      gow.w<-gow.w[-t]
    }
  }
  
  # replace NAs in 'abundances' by zero
  abundances[which(is.na(abundances))]<- 0
  site.species<-abundances
  
  # compute gower distances between species; this is done for species 
  # across communities rather than on a site-by-site basis, as gower 
  # distances between pairs of species may otherwise vary by site depending
  # on the range in ordinal or numeric traits represented at each site.
  if(is.null(gow.w)){
    spdist<-as.matrix(gowdis(x=traits, ord=gow.ord))
  }else{
    spdist<-as.matrix(gowdis(x=traits, w=gow.w, ord=gow.ord))
  }
  
  # Initiate storage vectors for species richness, and total tree length
  spec.rich<- numeric(length=nrow(site.species))
  total.mst<-numeric(length=nrow(site.species))
  
  # cycle through one site (s) at a time
  for(s in 1:nrow(site.species)){
    #s<-1   # for testing only
    
    # identify which species occur at site s & save species richness
    specs<-which(site.species[s,]>0)
    spec.rich[s]<-length(specs)
    
    # subsample distance matrix accordingly
    keep<-dimnames(site.species)[[2]][specs]
    site.spdist<-as.dist(spdist[keep,keep])
    
    # compute minimum spanning tree
    mst<-spantree(d=site.spdist, toolong=0)
    
    # determine the tree's total length 
    total.mst[s]<-sum(mst$dist)
  }
  outputs<-c("Site", "TotalSpecies", "FunctionalRichness")
  results<-matrix(data=0, nrow<-nrow(site.species), ncol=length(outputs),
                  dimnames=list(dimnames(site.species)[[1]], outputs))
  RMresults<-data.frame(results)
  RMresults$Site<-dimnames(site.species)[[1]]
  RMresults$TotalSpecies<-spec.rich
  RMresults$FunctionalRichness<-total.mst
  return(RMresults)
}