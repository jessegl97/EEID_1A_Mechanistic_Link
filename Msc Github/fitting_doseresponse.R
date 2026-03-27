library(tidyverse)

rm(list=ls()) # clear workspace
#dev.off() # close graphical device (if necessary)
cat("\014") # clear console

source("likelihood_functions.R")

actual = read.csv("day_41_removed_hofi_small_for_fitting_secondary.csv")
actual$primary_dose_scale = actual$primary_dose*0.001
actual$secondary_dose_scale = actual$secondary_dose*0.001

# no prior
actual_noprior = actual %>%
  filter(primary_dose==0)
actual_noprior = arrange(actual_noprior,by=secondary_dose_scale)

dose = actual_noprior$secondary_dose_scale
pos.c = actual_noprior$infected
tot.c = actual_noprior$N

hom.c<- optim(par=.01, fn=nll.hom, gr=NULL, method="Brent", lower=0,upper=1)#, control=list(trace=10))
bet.c<- optim(par=c(1,1), fn=nll.bet, gr=NULL, method="Nelder-Mead", control=list(trace=10))
gam.c<- optim(par=c(1,1), fn=nll.gam, gr=NULL, method="Nelder-Mead", control=list(trace=10))

# low

actual_lowprior = actual %>%
  filter(primary_dose==750)
actual_lowprior = arrange(actual_lowprior,by=secondary_dose_scale)

dose = actual_lowprior$secondary_dose_scale
pos.c = actual_lowprior$infected
tot.c = actual_lowprior$N

hom.l<- optim(par=.01, fn=nll.hom, gr=NULL, method="Brent", lower=0,upper=1)#, control=list(trace=10))
bet.l<- optim(par=c(1,1), fn=nll.bet, gr=NULL, method="Nelder-Mead", control=list(trace=10))
gam.l<- optim(par=c(1,1), fn=nll.gam, gr=NULL, method="Nelder-Mead", control=list(trace=10))


# high

actual_highprior = actual %>%
  filter(primary_dose==30000)
actual_highprior = arrange(actual_highprior,by=secondary_dose_scale)

dose = actual_highprior$secondary_dose_scale
pos.c = actual_highprior$infected
tot.c = actual_highprior$N

hom.h<- optim(par=.01, fn=nll.hom, gr=NULL, method="Brent", lower=0,upper=1)#, control=list(trace=10))
bet.h<- optim(par=c(10,10), fn=nll.bet, gr=NULL, method="Nelder-Mead", control=list(trace=10))
gam.h<- optim(par=c(1,1), fn=nll.gam, gr=NULL, method="Nelder-Mead", control=list(trace=10))


AppendData = function(x,y,z,exposure){
  df = data.frame(par1.gamma =z$par[1],par2.gamma=z$par[2],dev.gamma=z$value, par1.beta=y$par[1],par2.beta=y$par[2],dev.beta=y$value,par1.hom=x$par,dev.hom=x$value,LRTbeta = y$value-x$value,LRTgamma=z$value-x$value, group = exposure)
  return(df)
}
df1 <- AppendData(hom.c,bet.c,gam.c,"0 va")
df1 <- rbind(df1,AppendData(hom.l,bet.l,gam.l,"750 va"))
df1 <- rbind(df1,AppendData(hom.h,bet.h,gam.h,"30000 va"))
df1

write.csv(df1,file="deviance_params_noday41positives_03JAN2024.csv")







