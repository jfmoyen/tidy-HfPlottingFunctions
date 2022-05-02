### Hf-in-zrc processing functions

# Constants declaration



# Useful value
LuHf <- list(basalt=0.021,CC=0.0113)

# Initialize
earth.model<-function(age.to=0,age.from=4500,age.by=10,
                      lambda=1.87e-11,chur_I=0.282785,chur_R=0.0336){
  #' Calculate a range of useful values, to be used in subsequent calculations
  #' @param age.from, age.to, age.by: define the ages (Ma since present) for which to calculate
  #' @param lambda: Hf decay constant
  #' @param CHUR_I: Present-day 176Hf/177Hf for CHUR
  #' @param CHUR_R: Present-day 176Lu/177Hf for CHUR
  
  # Define the age sequence (0 is present)
  age.rng<-seq(age.to,age.from,by=age.by)
  
  # Calculate the evolution of the CHUR
  chur<-chur_I-chur_R*(exp(lambda*age.rng*1e6)-1)
  names(chur)<-age.rng
  
  return(list(age.rng=age.rng,chur=chur))
}


# Plot evolution line
evol.line<-function(time0,R0,eps0,dt=10){
  #' Plot an evolution line in a eHf - time diagram
  #' @param time0 Age (Ma) of the evolution start
  #' @param R0 176Lu/177Hf of the evolving reservoir
  #' @param eps0 Starting epsilon value
  #' @param dt time step (default is 10 Ma)
  #' @param earth an Earth model containing at least age range (for which to calculate) and CHUR, as defined by earth.model()
  
  deltaT<-seq(0,time0,dt)
  I0<-(1e-4*eps0+1)*chur[as.character(time0)]
  
  I<-I0+R0*(exp(lambda*deltaT*1e6)-1)
  I<-rev(I)
  names(I)<-age.rng[1:length(deltaT)]
  
  eps<-(I/chur-1)*1e4
  from<-length(deltaT)+1
  to<-length(age.rng)
  
  eps[from:to]<-NA
  names(eps)<-age.rng
  
  return(eps)
}