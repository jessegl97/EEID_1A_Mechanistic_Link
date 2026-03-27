##### Homogeneous #####

# homogeneous
pred.Exp<-function(theta, dose) {
  f<- 1-exp(-theta*dose)
  return(f)
}

# likelihood for homogeneous model
Lik.fun.hom=function(theta,dose) {# function for calculating the negative log likelihood
  pi_pred.c=pred.Exp(theta=theta,dose=dose)
  pi_obs.c<- pos.c/tot.c
  # Deviance calcs - control (Haas 8.33)
  Y1<- pos.c*log((pi_pred.c+1e-15)/(pi_obs.c+1e-15))
  Y2<- (tot.c-pos.c)*log((1-pi_pred.c+1e-15)/(1-pi_obs.c+1e-15))
  Y<- (-2)*(sum(Y1)+sum(Y2))
  return(Y)
}

# likelihood for homogeneous model
nll.hom <- function (par) {
  Lik.fun.hom(theta=par[1], dose=dose)
}

##### BETA - POISSON #####

# approximate beta Poisson
pred.betaPoisson<-function(alpha, betap, dose){
  f<-1-(1+(dose/betap))^(-alpha)
  return(f)
}

# likelihood for beta poisson
Lik.fun.beta=function(alpha,betap,dose) {# function for calculating the negative log likelihood
  # print(c(alpha,betap)) #print current values for debugging
  pi_pred=pred.betaPoisson(alpha=alpha,betap=betap,dose=dose)
  pi_obs<-pos.c/tot.c
  # Deviance calcs - control (Haas 8.33)
  Y1<- pos.c*log((pi_pred+1e-15)/(pi_obs+1e-15))
  Y2<- (tot.c-pos.c)*log((1-pi_pred+1e-15)/(1-pi_obs+1e-15))
  Y.c<- (-2)*(sum(Y1)+sum(Y2))
  return(Y.c)
}

# likelihood for beta model
nll.bet <- function (par) {
  Lik.fun.beta(alpha=par[1],betap=par[2], dose=dose)#
}

##### GAMMA #####

#gamma - in terms of shape and scale
# k is shape, omega is scale 
pred.inf.gamma<-function(theta,k,omega,dose) {
  ((theta^(k-1))*(exp(-theta/omega)))/((omega^k)*(gamma(k)))*(exp(-theta*dose))
}


#integration function for gamma - shape and scale
int2<-function(k,omega,dose){
  n <- length(dose)
  ig=NA;id=NA
  for (i in 1:n){
    ig=integrate(pred.inf.gamma,0,Inf,k=k,omega=omega,dose=dose[i],stop.on.error = FALSE)$value
    id=rbind(id,ig)
    pi_pred=1-as.vector(id)
    pi_pred=pi_pred[!is.na(pi_pred)]
  } 
  return(pi_pred)
}

Lik.fun.gamma=function(k,omega,dose) {# function for calculating the negative log likelihood
  # print(c(k,omega)) #print current values for debugging
  pi_pred=int2(k,omega,dose)
  pi_obs<- pos.c/tot.c
  # Deviance calcs - control (Haas 8.33)
  Y1<- pos.c*log((pi_pred+1e-15)/(pi_obs+1e-15))
  Y2<- (tot.c-pos.c)*log((1-pi_pred+1e-15)/(1-pi_obs+1e-15))
  Y<- (-2)*(sum(Y1)+sum(Y2))
  return(Y)
}

##gamma model with shape and scale
nll.gam <- function (par) {
  Lik.fun.gamma(k=par[1],omega=par[2], dose=dose)#
}
