##### DOSE RESPONSE BOOTSTRAPPED CIs ##### 

# requires the following files:
#     day_41_removed_hofi_small_for_fitting_secondary.csv
#     deviance_params_noday41positives.csv
#     likelihood_functions.R

# produces the following files (if lines uncommented):
#     all_psuedosets.csv
#     va30000_gamma_distributionCIs_new.csv (line 295)
#     va750_gamma_distributionCIs_new.csv (line 375)
#     va0_gamma_distributionCIs_new.csv (line 441)

rm(list=ls()) # clears workspace
#dev.off() # close graphical device (if necessary)
cat("\014") # clear console

# load important packages
library(base)
library(dplyr)
library(ellipse)
library(ggplot2)
library(MASS)
library(Matrix)
library(reshape)
library(reshape2)
library(sp)
library(stats)
library(tidyr)
library(utils)

# Load data generated from fitting_doseresponse.R
store2 = read.csv("deviance_params_noday41positives.csv")
z3 = read.csv("day_41_removed_hofi_small_for_fitting_secondary.csv")
source("likelihood_functions.R")

z3 = arrange(z3,primary_dose,secondary_dose)
z3$dose_scale = z3$secondary_dose/1000
z3$frac = z3$infected/z3$N
# log dose
z3$ldose = log10(z3$secondary_dose+1)
# group dose
z3$pdose.f = as.factor(z3$primary_dose)
# create 'bird groups'
z3$bird.groups = paste(z3$pdose.f, z3$population)
unique.birds = unique(z3$bird.groups)


# Homogeneous model
# Exponential - va 0
va.0.preds=pred.Exp(store2$par1.hom[store2$group=="0 va"],
                    z3$dose_scale[z3$primary_dose=="0"])

va.0 = z3 %>%
  filter(primary_dose=="0")%>%
  select(primary_dose, infected, N, secondary_dose, dose_scale)%>%
  mutate(preds = va.0.preds)%>%
  #calculate chi-squared residual at each dose (Equation 9.10 in Haas)
  mutate(err = 
           (infected/N - preds)/sqrt(preds*(1-preds)/N))
va.0$err[va.0$preds==0]=0
print(va.0)

# Heterogeneous model
# Gamma -  va 750
va.mid.preds=int2(store2$par1.gamma[store2$group=="750 va"],
                  store2$par2.gamma[store2$group=="750 va"],
                  z3$dose_scale[z3$primary_dose=="750"])

va.mid = z3 %>%
  filter(primary_dose=="750")%>%
  select(primary_dose, infected, N, secondary_dose, dose_scale)%>%
  mutate(preds = va.mid.preds)%>%
  #calculate chi-squared residual at each dose (Equation 9.10 in Haas)
  mutate(err = 
           (infected/N - preds)/sqrt(preds*(1-preds)/N))

va.mid$err[va.mid$preds==0]=0
print(va.mid)


# Heterogeneous
# Gamma - va 30000
va.high.preds=int2(store2$par1.gamma[store2$group=="30000 va"],
                   store2$par2.gamma[store2$group=="30000 va"],
                   z3$dose_scale[z3$primary_dose=="30000"])

va.high = z3 %>%
  filter(primary_dose=="30000")%>%
  select(primary_dose, infected, N, secondary_dose, dose_scale)%>%
  mutate(preds = va.high.preds)%>%
  #calculate chi-squared residual at each dose (Equation 9.10 in Haas)
  mutate(err = 
           (infected/N - preds)/sqrt(preds*(1-preds)/N))

va.high$err[va.high$preds==0]=0
print(va.high)

# Combine data
va.comb = bind_rows(va.high, va.0, va.mid)%>%
  drop_na(err)
va.comb$group.dose = paste(paste(va.comb$primary_dose,"va",sep=" "), va.comb$dose_scale, sep="_")
print(va.comb)

# pre-allocate for prediction
R=1000
pi_m.v=NA;pi_m.v.st=NA;dose.st=NA;pi_m.c=NA;pi_m.c.st=NA;

for (i in 1:length(unique(va.comb$group.dose))){
  print(i)
  #for each dose within each group
  for (j in 1:R) {
    #at one dose, sample all the errors for that group 1000 times 
    #and calculate a new predicted value for that dose
    # Equation 9.11 in Haas
    # pi_i^(m) = pi_i^(T) + eps_q * sqrt((pi_i^(T) *(1-pi_i^(T)))/n_i)
    # where   pi_i^(T) are the best fit proportions to the original data, i.e. va.comb$preds
    #         eps_q is a random index of the errors in defined above, i.e. va.comb$err
    #         n_i is the number of samples for that condition, i.e. va.comb$N
    pi_m.v=va.comb$preds[i]+sample(va.comb$err,1)*sqrt((va.comb$preds[i]*(1-va.comb$preds[i]))/va.comb$N[i])
    pi_m.v[pi_m.v<=0]=0
    pi_m.v[pi_m.v>=1]=1
    pi_m.v.st=c(pi_m.v.st,pi_m.v)
  }
}
# Need to remove the NA, the first observation
pi_m.v.st <- pi_m.v.st[!is.na(pi_m.v.st)]

# Pre-tidyverse - everything needs to stay in order
group.dose.st = rep(va.comb$group.dose,each=1000)
N.st = rep(va.comb$N,each=1000)

# Separate the doses and groups for the next fitting
x=colsplit(group.dose.st, "_",c( "bird.groups.st","dose.st"))

# Combine the vectors together to make a new dataframe for modelling
mframe=data.frame(pi_m.v.st,N.st, x, group.dose.st)
head(mframe)

# pre-allocate for estimating positives
mframe$pos.c=NA

# Estimate number of positives-based on the probabilities
for(i in 1:length(mframe$pi_m.v.st))  {
  set1=rbinom(prob=mframe$pi_m.v.st[i],n=mframe$N.st[i],size=1)
  mframe$pos.c[i]=sum(set1)
}  
# Total number in format to use likelihood functions
mframe$tot.c=mframe$N.st

# Counts 16 times because there are 16 conditions and need to pick the same dose and the right number of values to fit to
# Note: it would be clearer to group by the dose.groups
mframe$bs=rep(1:1000,times=16)

# Paste together in order to get the groups and sequence to line up
mframe$count.me = paste(mframe$bird.groups.st, mframe$bs, sep = "_")

head(mframe);tail(mframe)


## Fitting models to pseudo-datasets

# pre-allocate for fitting
store.pseudo.preds=data.frame(
                par1.gamma=NA, par2.gamma=NA, dev.gamma=NA,
                 par1.beta=NA, par2.beta=NA, dev.beta=NA, 
                 par1.hom =NA, dev.hom = NA,
                 group=NA)
store.pseudo.preds2=store.pseudo.preds
unique.counts = unique(mframe$count.me)

# fitting
for (i in 1:length(unique(mframe$count.me))){
  print(i)
  pos.c = mframe$pos.c[mframe$count.me==unique.counts[i]]
  tot.c = mframe$tot.c[mframe$count.me==unique.counts[i]]
  dose = mframe$dose.st[mframe$count.me==unique.counts[i]]
  group = unique.counts[i]
  gamma.o<- optim(par=c(1,1), fn=nll.gam, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  beta.o<- optim(par=c(1,1), fn=nll.bet, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  hom.o<- optim(par=c(.01), fn=nll.hom, gr=NULL, method="Brent",lower=0,upper=1, control=list(trace=10))
  store.pseudo.preds = data.frame(
                    par1.gamma = gamma.o$par[1], 
                     par2.gamma = gamma.o$par[2], 
                     dev.gamma = gamma.o$val[1],
                     par1.beta = beta.o$par[1], 
                     par2.beta = beta.o$par[2], 
                     dev.beta = beta.o$val[1],
                     par1.hom = hom.o$par[1], 
                     dev.hom = hom.o$val[1],
                     group = group)
  store.pseudo.preds2=bind_rows(store.pseudo.preds2,store.pseudo.preds)
}

p=colsplit(store.pseudo.preds2$group, "_",c( "bird.groups","counter"))
store.pseudo.preds2 = cbind(store.pseudo.preds2, p)
unique(store.pseudo.preds2$bird.groups)

ps = store.pseudo.preds2
ps = ps %>%
  drop_na()

write.csv(ps, "all_pseudosets.csv",row.names = FALSE)

cis = ps %>%
  drop_na()%>%
  group_by(bird.groups)%>%
  summarise(ci_beta_par1_low = quantile(par1.beta,0.025 ),
            ci_beta_par1_high = quantile(par1.beta,0.975 ),
            ci_beta_par2_low = quantile(par2.beta,0.025 ),
            ci_beta_par2_high = quantile(par2.beta,0.975 ),
            ci_gam_par1_low = quantile(par1.gamma,0.025 ),
            ci_gam_par1_high = quantile(par1.gamma,0.975 ),
            ci_gam_par2_low = quantile(par2.gamma,0.025 ),
            ci_gam_par2_high = quantile(par2.gamma,0.975 ),
            mean_gam_low = quantile(par1.gamma,0.025 )*quantile(par2.gamma,0.025 ),
            mean_gam_high = quantile(par1.gamma,0.975 )*quantile(par2.gamma,0.975 ),
            ci_hom_low = quantile(par1.hom,0.025 ),
            ci_hom_high = quantile(par1.hom,0.975 )
  )
head(cis)            
# Note: beta bounds are so wide because alpha and beta depend on each other
# Use the ellipse method

#write.csv(ps, "figs/all_pseudosets.csv",row.names = FALSE)

# No prior exposure
# homogeneous model
hom.vals = ps %>%
  filter(bird.groups == "0 va")%>%
  select(par1.hom, bird.groups, dev.hom) %>%
  mutate(ci_hom_low = quantile(par1.hom,0.025 )) %>%
  mutate(ci_hom_high = quantile(par1.hom,0.975 ))%>%
  filter(par1.hom>=ci_hom_low&par1.hom<=ci_hom_high)

# If want to save new distributions for confidence intervals, uncomment next line
#write.csv(hom.vals, "va0_CIs_hom_pseudosets_new.csv",row.names = FALSE)

# No prior exposure
# heterogeneous model 
# gamma distribution
p0=ggplot(data=subset(ps, bird.groups=="0 va"),aes(x=par1.gamma,y=par2.gamma))+
  geom_jitter(size=1)+
  stat_ellipse(level=0.95)+
  ylab("beta")+
  xlab("alpha")+
  xlim(0,100)+
  ylim(0,100)+
  theme_bw() +
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=23),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=20,face="italic"),legend.title = element_blank(),legend.background = element_blank(),legend.key=element_rect(fill="white",color="white"))
p0

# Extract components
build <- ggplot_build(p0)$data
points <- build[[1]]
ell <- build[[2]]

# Determine which points are inside the ellipse and this info add to data frame
dat2 <- data.frame(points[1:2],
                   in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))

# Plot points, separating in and out of ellipse
ggplot(dat2, aes(x, y)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()

# Select parameter sets in ellipse
dat3=dat2[dat2$in.ell=="TRUE",]
nrow(dat3)

dat3$group = "0 va"
dat3 = dat3 %>%
  dplyr::rename(par1 = x, par2 = y)

# Generate gamma distribution data using selected parameter sets
# pre-allocate
num=NA
name=NA
hx2=matrix(NA,ncol=3, nrow=1000)
hx2names=c("num","name","x")
hx2=as.data.frame(hx2)
names(hx2)=hx2names
hx=NA;hx.i.true=NA;hx.i.est=NA;hx.true=NA;hx.est=NA;hx.name.i=NA;hx.name=NA;x.i=NA
hx.x=NA;hx.alpha.st.i=NA;hx.alpha.st=NA;hx.betap.st.i=NA;hx.betap.st=NA;

# Generate gamma distributions
for (i in 1:nrow(dat3)) {
  x.i=seq(0,5,length=300)
  hx.x=c(hx.x,x.i)
  hx.i.est=dgamma(x.i,dat3$par1[i],dat3$par2[i])
  hx.est=c(hx.est,hx.i.est)
}
hx.est <- hx.est[!is.na(hx.est)]
hx.x <- hx.x[!is.na(hx.x)]
name2=rep("boot",times=length(hx.est))
gen.dist=data.frame(hx.x,hx.est,name2)
names(gen.dist)=c("x","hx","name")

# Plot output of the distributions
ggplot(data=gen.dist, aes(x=hx.x,y=hx.est ))+
  geom_line()+
  coord_cartesian(ylim=c(0,100))

# If want to save new gamma distributions for confidence intervals, uncomment next line
write.csv(gen.dist, "va0_gamma_distributionCIs_new.csv",row.names = FALSE)

## Low dose prior exposure (VA 750)
# heterogeneous model 
# gamma distribution
p750=ggplot(data=subset(ps, bird.groups=="750 va"),aes(x=par1.gamma,y=par2.gamma))+
  geom_jitter(size=1)+
  stat_ellipse(level=0.95)+
  ylab("beta")+
  xlab("alpha")+
  xlim(0,100)+
  ylim(0,100)+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=23),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=20,face="italic"),legend.title = element_blank(),legend.background = element_blank(),legend.key=element_rect(fill="white",color="white"))
p750


# Extract components
build <- ggplot_build(p750)$data
points <- build[[1]]
ell <- build[[2]]

# Determine which points are inside the ellipse and this info add to data frame
dat2 <- data.frame(points[1:2], 
                  in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))

# Plot points, separating in and out of ellipse
ggplot(dat2, aes(x, y)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()

# Select parameter sets in ellipse
dat3=dat2[dat2$in.ell=="TRUE",]
nrow(dat3)
dat3$group = "750 va"
dat3 = dat3 %>%
  dplyr::rename(par1 = x, par2 = y)

# Generate gamma distribution data using selected parameter sets
# pre-allocate
num=NA
name=NA
hx2=matrix(NA,ncol=3, nrow=1000)
hx2names=c("num","name","x")
hx2=as.data.frame(hx2)
names(hx2)=hx2names
hx=NA;hx.i.true=NA;hx.i.est=NA;hx.true=NA;hx.est=NA;hx.name.i=NA;hx.name=NA;x.i=NA
hx.x=NA;hx.alpha.st.i=NA;hx.alpha.st=NA;hx.betap.st.i=NA;hx.betap.st=NA;

x=seq(0.0001,5,length=300)
group.labels = rep(unique.birds, each=length(x))
dist.dat = data.frame(x = rep(x, times=length(unique.birds)), 
                      group = group.labels)
dist.dat$par1.gamma = store2$par1.gamma[match(dist.dat$group, store2$group)]
dist.dat$par2.gamma = store2$par2.gamma[match(dist.dat$group, store2$group)]
dist.dat$hxa=dgamma(dist.dat$x,shape=dist.dat$par1.gamma,scale=dist.dat$par2.gamma)

# Calculate mean
dist.dat$mean_gamma = dist.dat$par1.gamma*dist.dat$par2.gamma
unique(dist.dat$mean_gamma)

# Generate gamma distributions
for (i in 1:nrow(dat3)) {
  x.i=seq(0,5,length=300)
  hx.x=c(hx.x,x.i)
  hx.i.est=dgamma(x.i,dat3$par1[i],dat3$par2[i])
  hx.est=c(hx.est,hx.i.est)
}
hx.est <- hx.est[!is.na(hx.est)]
hx.x <- hx.x[!is.na(hx.x)]
name2=rep("boot",times=length(hx.est))
gen.dist=data.frame(hx.x,hx.est,name2)
names(gen.dist)=c("x","hx","name")

# Plot output of distributions
ggplot(data=gen.dist, aes(x=hx.x,y=hx.est ))+
  geom_line()+
  coord_cartesian(ylim=c(0,100))

# If want to save new gamma distributions for confidence intervals, uncomment next line
write.csv(gen.dist, "va750_gamma_distributionCIs_new.csv",row.names = FALSE)

## High dose prior exposure (VA 30000)
# heterogeneous model 
# gamma distribution
phigh=ggplot(data=subset(ps, bird.groups=="30000 va"),aes(x=par1.gamma,y=par2.gamma))+
  geom_jitter(size=1)+
  stat_ellipse(level=0.95)+
  ylab("beta")+
  xlab("alpha")+
  xlim(0,100)+
  ylim(0,100)+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=23),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=20,face="italic"),legend.title = element_blank(),legend.background = element_blank(),legend.key=element_rect(fill="white",color="white"))
phigh

# Extract components
build <- ggplot_build(phigh)$data
points <- build[[1]]
ell <- build[[2]]

# Determine which points are inside the ellipse and this info add to data frame
dat2 <- data.frame(points[1:2], 
                   in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))

# Plot points, separating in and out of ellipse
ggplot(dat2, aes(x, y)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()

# Select parameter sets in ellipse
dat3=dat2[dat2$in.ell=="TRUE",]
nrow(dat3)
dat3$group = "30000 va"
dat3 = dat3 %>%
  rename(par1 = x, par2 = y)

# Generate gamma distribution
# pre-allocate
num=NA
name=NA
hx2=matrix(NA,ncol=3, nrow=1000)
hx2names=c("num","name","x")
hx2=as.data.frame(hx2)
names(hx2)=hx2names
hx=NA;hx.i.true=NA;hx.i.est=NA;hx.true=NA;hx.est=NA;hx.name.i=NA;hx.name=NA;x.i=NA
hx.x=NA;hx.alpha.st.i=NA;hx.alpha.st=NA;hx.betap.st.i=NA;hx.betap.st=NA;

# Generate distributions
for (i in 1:nrow(dat3)) {
  x.i=seq(0,5,length=300)
  hx.x=c(hx.x,x.i)
  hx.i.est=dgamma(x.i,dat3$par1[i],dat3$par2[i])
  hx.est=c(hx.est,hx.i.est)
}
hx.est <- hx.est[!is.na(hx.est)]
hx.x <- hx.x[!is.na(hx.x)]
name2=rep("boot",times=length(hx.est))
gen.dist=data.frame(hx.x,hx.est,name2)
names(gen.dist)=c("x","hx","name")

# Plot output of distributions
ggplot(data=gen.dist, aes(x=hx.x,y=hx.est ))+
  geom_line()+
  coord_cartesian(ylim=c(0,100))

# If want to save new gamma distributions for confidence intervals, uncomment next line
write.csv(gen.dist, "va30000_gamma_distributionCIs_new.csv",row.names = FALSE)

################### Go to fitting_doseresponse.R to use confidence intervals in figures ##################
