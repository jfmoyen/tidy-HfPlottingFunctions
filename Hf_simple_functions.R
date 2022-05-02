# CHUR Bouvier et al 2008 : 176Lu/177Hf = 0.0336 and 176Hf/177Hf = 0.282785
# CHUR Iizuka et al 2015  : 176Lu/177Hf = 0.0338 and 176Hf/177Hf = 0.279781
# CHUR Blich.T et al 1997 : 176Lu/177Hf = 0.0332 and 176Hf/177Hf = 0.282772

# DM Naeraa 2012 : 176Lu/177Hf = 0.0375 and 176Hf/177Hf = 0.283120
# DM Griffin 2000: 176Lu/177Hf = 0.0384 and 176Hf/177Hf = 0.283250
# DM Griffin 2004: 176Lu/177Hf = 0.0384 and 176Hf/177Hf = 0.279718  ???

# Typical values 177Lu/176Hf are 0.0093 for upper CC, 
# 0.0113 for continental crust, 0.022 for mafic rocks

# Guestimate Guitreau 2012 Mantle source
# present day epsilon = +2 to +5
# R = 0.0336 (chondrite)
# Dans G et al '12, CHUR after BT 1997 --> 0.0332 et réservoir source -> 0.0338
# Avec un CHUR à .0336 il faut "guesstimer" le réservoir à 0.0342 ou à peu près

CHUR_t<-function(age,lambda=1.867e-11,chur_I0=0.282785,chur_R0=0.0336){
  #' Calculate the value of CHUR at a given age
  #' @param age: Age (Ma since present) for which to calculate
  #' @param lambda: Hf decay constant
  #' @param CHUR_I0: Present-day 176Hf/177Hf for CHUR
  #' @param CHUR_R0: Present-day 176Lu/177Hf for CHUR 
  #' Default CHUR values are from Bouvier et al. 2008
  #' lambda from Söderlund et al. 2004

   chur<-chur_I0-chur_R0*(exp(lambda*age*1e6)-1)
  
   return(chur)
}

I_to_epsilon<-function(I,age,...){
  #' Calculate the epsilon value for a given 177/176Hf at age
  #' @param age: Age (Ma since present) for which to calculate
  #' @param I: 176Hf/177Hf ratio at time t.
  #' @param ...: further parameters to be passed to CHUR_t (possibly, alternative CHUR models)
  
  eps <- (I/CHUR_t(age,...)-1)*1e4
  
  return(eps)
}

epsilon_to_I<-function(epsilon,age,...){
  #' Calculate the 177/176Hf value for a given epsilon at age
  #' @param age: Age (Ma since present) for which to calculate
  #' @param epsilon: epsilon at time t.
  #' @param ...: further parameters to be passed to CHUR_t (possibly, alternative CHUR models)
  
  I <- CHUR_t(age,...)*(1+epsilon*1e-4)
 
  return(I)
}

# Note to self
##############
#
# "back" and "forward" calculations have similar but not identical equations !
# It = I0 - R0 * (e^lambda*t - 1)
# I0 = It + Rt * (e^lambda*t - 1)
# There does not appear to be a generalized form so we write two functions


Backcalculate_I<-function(I0,R0,lambda=1.87e-11,age){
  #' Calculate the initial 176Hf/177Hf from present day values
  #' @param I0: present day 176Hf/177Hf
  #' @param R0: present day 176Lu/177Hf
  #' @param lambda: Hf decay constant
  #' @param age: Age (Ma since present) for which to calculate
   
  It <- I0 - R0 * (exp(lambda*age*1e6) - 1)
  
  return(It)
}

Forwardcalculate_I<-function(It,Rt,lambda=1.87e-11,age){
  #' Calculate the initial 176Hf/177Hf from present day values
  #' @param It: initial 176Hf/177Hf
  #' @param Rt: initial 176Lu/177Hf
  #' @param lambda: Hf decay constant
  #' @param age: Age (Ma after initial) for which to calculate
  
  I0 <- It + Rt * (exp(lambda*age*1e6) - 1)
  
  return(I0)
}

Evol_from_known_point<-Vectorize(function(age_start,
                              epsilon_start=0,I_start=epsilon_to_I(epsilon_start,age_start,...),
                              R_start,
                              age,
                              return_epsilon=T,
                              lambda=1.87e-11,...){
  #' Calculate the evolution (i.e. crust evolution line) from a point in time
  #' where we know, or assume, epsilon and R
  #' @param age_start: the starting point
  #' @param epsilon_start: epsilon value at starting time
  #' @param I_start: 176/177Hf at starting time. Calculated from epsilon of not supplied.
  #' @param R_start: 177Lu/176Hf at starting time
  #' Typical values are 0.0093 for upper CC, 0.0113 for continental crust, 0.022 for mafic rocks
  #' @param age: age of calculation
  #' @param return_epsilon: If true, the result is returned as epsilon, else as 176Hf/177Hf
  #' @param lambda: Hf decay constant
  #' @param ...: further parameters to be passed to CHUR_t 
  #' via epsilon_to_I and I_to_epsilon (possibly, alternative CHUR models)

  if(age > age_start){
    # Calculating before the start point is meaningless...
    retval <- NA
  }else{
    I <- Forwardcalculate_I(It=I_start,Rt=R_start,lambda=lambda,age=age_start - age)
    
    # If required, convert to epsilon values
    if(return_epsilon){
      retval <- I_to_epsilon(I=I,age=age,...)
    }else{
      retval <- I
    }
  }
  
  return(retval)
}
)