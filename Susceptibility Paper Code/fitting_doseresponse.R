##### DOSE RESPONSE FITTING AND PLOTTING ##### 

# requires the following files:
#     day_41_removed_hofi_small_for_fitting_secondary.csv
#     likelihood_functions.R
#     va30000_gamma_distributionCIs.csv
#     va750_gamma_distributionCIs.csv
#     va0_gamma_distributionCIs.csv

# produces the following files:
#     deviance_params_noday41positives.csv
#     Fig/va_dose_response_fits_error_bars.png
#     Fig/heterogeneous_gamma_fits_dose-response.png
#     Fig/homogeneous_fits_dose-response.png


rm(list=ls()) # clears workspace
#dev.off() # close graphical device (if necessary)
cat("\014") # clear console

# load important packages
library(base)
library(binom)
library(deSolve)
library(dplyr)
library(forcats)
library(ggplot2)
library(lme4)
library(MASS)
library(Matrix)
library(reshape)
library(reshape2)
library(tibble)
library(tidyr)
library(utils)


#1. READ IN DATA, PREP, AND CREATE GROUPS
dat = read.csv("day_41_removed_hofi_small_for_fitting_secondary.csv")
z3 = dat


#2.0 DEFINING THE FUNCTIONS
source("likelihood_functions.R")


#2.1 FITTING THE MODELS

# scale the dose
z3$dose_scale = z3$secondary_dose*0.001
# rename infected and N for use likelihood functions
z3$pos = z3$infected
z3$tot = z3$N
# log dose
z3$ldose = log10(z3$secondary_dose+1)
# group dose
z3$pdose.f = as.factor(z3$primary_dose)
# create 'bird groups'
z3$bird.groups = paste(z3$pdose.f, z3$population)

unique.birds = unique(z3$bird.groups)
length.bird.groups = length(unique(z3$bird.groups))

# arrange data by groups and secondary dose
z3 = z3 %>%
  group_by(bird.groups)%>%
  arrange(secondary_dose)

# pre-allocate for dose fitting
store=data.frame(par1.gamma=NA, par2.gamma=NA, dev.gamma=NA,
                 par1.hom =NA, dev.hom = NA,
                 group=NA)
store2=store

# dose fitting
for (i in 1:length.bird.groups){
  pos.c = z3$pos[z3$bird.groups==unique.birds[i]]
  tot.c = z3$tot[z3$bird.groups==unique.birds[i]]
  dose = z3$dose_scale[z3$bird.groups==unique.birds[i]]
  group = unique.birds[i]
  gamma.o<- optim(par=c(1,1), fn=nll.gam, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  hom.o<- optim(par=c(.01), fn=nll.hom, gr=NULL, method="Brent",lower=0,upper=1, control=list(trace=10))
  store = data.frame(par1.gamma = gamma.o$par[1], 
                     par2.gamma = gamma.o$par[2], 
                     dev.gamma = gamma.o$val[1],
                     par1.hom = hom.o$par[1], 
                     dev.hom = hom.o$val[1],
                     group = group)
  store2=bind_rows(store2,store)
}

store2 = store2 %>%
  drop_na()

store2$gamma.mean = store2$par1.gamma*store2$par2.gamma
store2$gamma.cov = sqrt(store2$par1.gamma*store2$par2.gamma^2)/store2$gamma.mean

view(store2)

write.csv(store2, "deviance_params_noday41positives.csv", row.names = FALSE)

##3. PLOT MODEL OUTPUT

z3$frac = z3$pos/z3$tot

# add binonmial standard error - normal approximation
z3$se=(z3$frac*(1-z3$frac)/z3$N)^.5
# add 2 successes and 2 failures correction to the normal approximation when p = 0 and p = 1
z3$pos2 = z3$pos+2 # add 2 successes 
z3$tot2 = z3$tot+4 # add 2 successes and 2 failures
z3$frac2 = z3$pos2/z3$tot2
z3$se[z3$se==0|z3$se==1]=(z3$frac2*(1-z3$frac2)/z3$tot2)^.5
z3$se_lower = z3$frac-z3$se
z3$se_upper = z3$frac+z3$se



va750 = z3 %>%
  filter(bird.groups == "750 va")

# create a smooth dose sequence 
dose.seq=seq(from=min(z3$dose_scale),to=max(z3$dose_scale), by=.01)
ct=seq(from=1, to=length(dose.seq))
data2=data.frame(dose.seq,ct)
data2$dose=data2$dose.seq

data3 = data2 %>%
  dplyr::slice(rep(1:nrow(data2), times=length(unique.birds)) )
group.labels = rep(unique.birds, each=nrow(data2))
data3$group = group.labels  

# predict models
pred_gamma=NA; pred_hom=NA; pred_group=NA; preds_group_i=NA
current_dose=NA; current_dose_i=NA
for( i in 1:nrow(store2)){
  preds_1=int2(store2$par1.gamma[i],
               store2$par2.gamma[i],
               data3$dose[data3$group==store2$group[i]])
  
  pred_gamma = c(preds_1, pred_gamma)
  pred_group_i  = rep(store2$group[i], length(preds_1))
  pred_group  = c(pred_group_i, pred_group)
  current_dose_i = data3$dose[data3$group==store2$group[i]]
  current_dose = c(current_dose_i, current_dose)

  preds_3=pred.Exp(store2$par1.hom[i],
                   data3$dose[data3$group==store2$group[i]])
  pred_hom = c(preds_3, pred_hom)
}

dat.tog = data.frame(pred_group,pred_hom, pred_gamma, current_dose)

data4 = dat.tog %>%
  pivot_longer(c(pred_hom, pred_gamma), names_to = "fit_type", values_to = "preds")

z3$group = z3$bird.groups
data4$group = data4$pred_group

data4$group[data4$group=="NA"]=NA
data4 = data4 %>%
  drop_na()
unique(data4$group)

# just virginia
va.pred = data4 %>%
  filter(group=="0 va"|group=="750 va"|group=="30000 va")

va = z3 %>%
  filter(group=="0 va"|group=="750 va"|group=="30000 va")
va$group = as.factor(va$group)
va = va %>%
  mutate(group = fct_relevel(group, "0 va", "750 va", "30000 va"))
va.pred = va.pred %>%
  mutate(group = fct_relevel(group, "0 va", "750 va", "30000 va"))

##3.1 DOSE RESPONSE PLOT

# add exposure labels to predictions
head(va.pred)
va.pred$exposure = "no prior exposure"
va.pred$exposure[va.pred$pred_group=="30000 va"] = "high prior exposure"
va.pred$exposure[va.pred$pred_group=="750 va"] = "low prior exposure"

va.pred$model = "heterogeneous"
va.pred$model[va.pred$fit_type=="pred_hom"] = "homogeneous"

# add exposure labels to data
head(va)
va$exposure = "no prior exposure"
va$exposure[va$group=="30000 va"] = "high prior exposure"
va$exposure[va$group=="750 va"] = "low prior exposure"

# reorder by exposure
str(va$exposure)
va$exposure = as.factor(va$exposure)
va = va %>%
  mutate(exposure = fct_relevel(exposure, "no prior exposure", 
                                "low prior exposure", 
                                "high prior exposure"))

va.pred = va.pred %>%
  mutate(exposure = fct_relevel(exposure, "no prior exposure", 
                                "low prior exposure", 
                                "high prior exposure"))

x.axis1 = expression(log[10]~MG~Exposure~Concentration~"("~CCU~per~mu*l~")")

# plot dose response with se
r.log.va.final2=ggplot(data=va, aes(x=dose_scale,y=frac))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=se_lower, ymax=se_upper), width=.1)+
  geom_line(data=subset(va.pred),
            aes(x=current_dose, y=preds, color = model, linetype = model),size=1)+
  facet_wrap(~exposure)+
  ylab("Fraction Infected")+
  xlab(x.axis1)+
  coord_cartesian(ylim=c(0,1))+
  scale_linetype_manual(values=c("solid", "dotted"))+
  scale_colour_manual(values = c("#068DA9", "#7E1717"))+
  scale_x_log10()+
  theme_bw() +
  #for facets
  theme(strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15),
        panel.grid = element_blank(), axis.line=element_line(),
        legend.position="top",legend.title=element_blank(),
        legend.text = element_text(size=20))
r.log.va.final2

ggsave(file="Fig/va_dose_response_fits_error_bars.png",width=13,height=6,units="in",limitsize=FALSE)


##3.2 PLOT DISTRIBUTIONS

##--------------gamma--------------------##
x=seq(0.0001,5,length=300)
group.labels = rep(unique.birds, each=length(x))
dist.dat = data.frame(x = rep(x, times=length(unique.birds)), 
                      group = group.labels)

dist.dat$par1.gamma = store2$par1.gamma[match(dist.dat$group, store2$group)]
dist.dat$par2.gamma = store2$par2.gamma[match(dist.dat$group, store2$group)]

dist.dat$hxa=dgamma(dist.dat$x,shape=dist.dat$par1.gamma,scale=dist.dat$par2.gamma)

# calculate mean
dist.dat$mean_gamma = dist.dat$par1.gamma*dist.dat$par2.gamma
unique(dist.dat$mean_gamma)

# bring in simulated data for CI ribbons
hi.rib = read.csv("va30000_gamma_distributionCIs_new.csv")
hi.rib$exposure = "high prior exposure"

mid.rib = read.csv("va750_gamma_distributionCIs_new.csv")
mid.rib$exposure = "low prior exposure"

no.rib = read.csv("va0_gamma_distributionCIs_new.csv")
no.rib$exposure = "no prior exposure"

# combine data
ribs = bind_rows(hi.rib, mid.rib, no.rib)

rib.cis = ribs %>%
  group_by(exposure,x)%>%
  summarise(high.gam = max(hx, na.rm=T), low.gam = min(hx, na.rm=T))


# final distributions
va.dist =  dist.dat %>%
  filter(group=="0 va"|group=="750 va"|group=="30000 va") 
va.dist = va.dist %>%
  mutate(group = fct_relevel(group, "0 va", "750 va", "30000 va"))

va.dist$exposure = "no prior exposure"
va.dist$exposure[va.dist$group=="30000 va"] = "high prior exposure"
va.dist$exposure[va.dist$group=="750 va"] = "low prior exposure"

va.dist$exposure = as.factor(va.dist$exposure)
va.dist = va.dist %>%
  mutate(exposure = fct_relevel(exposure, "no prior exposure", 
                                "low prior exposure", 
                                "high prior exposure"))

# round for left join
rib.cis$x = round(rib.cis$x, 3)
va.dist$x = round(va.dist$x,3)

rib.cis = rib.cis %>%
  filter(high.gam!=Inf)

va.dist.comb = left_join(va.dist, rib.cis)

va.dist.comb = va.dist.comb %>%
  mutate(exposure = fct_relevel(exposure, "no prior exposure", "low prior exposure", "high prior exposure"))

# plot gamma distributions
rgamma=ggplot(data = va.dist.comb %>% drop_na(high.gam), aes(x=x,y=hxa))+ 
  geom_ribbon(aes(ymin=low.gam,ymax=high.gam),col="#D1DBBD",fill="#D1DBBD",size=3)+
  geom_line(size=1, color ="#21918c" )+
  facet_wrap(~exposure)+
  coord_cartesian(ylim=c(0,1.2))+
  geom_vline(aes(xintercept = mean_gamma), color="gray")+
  ylab(expression("f(x)"))+
  xlab(expression("x, susceptibility"))+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=20),axis.text.y=element_text(size=15),axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.title=element_blank(),legend.text = element_text(size=20,face="italic"))
rgamma

ggsave(file="Fig/heterogeneous_gamma_fits_dose-response.png",width=12,height=6,units="in",limitsize=FALSE)
 
##--------------homogeneous--------------------##

x=seq(0.0001,1.01,length=300)
group.labels = rep(unique.birds, each=length(x))
dist.dat = data.frame(x = rep(x, times=length(unique.birds)), 
                      group = group.labels)

dist.dat$par1.hom = store2$par1.hom[match(dist.dat$group, store2$group)]

va.dist =  dist.dat %>%
  filter(group=="0 va"|group=="750 va"|group=="30000 va") 
va.dist = va.dist %>%
  mutate(group = fct_relevel(group, "0 va", "750 va", "30000 va"))

va.dist$exposure = "no prior exposure"
va.dist$exposure[va.dist$group=="30000 va"] = "high prior exposure"
va.dist$exposure[va.dist$group=="750 va"] = "low prior exposure"

va.dist$exposure = as.factor(va.dist$exposure)
va.dist = va.dist %>%
  mutate(exposure = fct_relevel(exposure, "no prior exposure", 
                                "low prior exposure", 
                                "high prior exposure"))

# plot homogeneous model
rhom=ggplot(data=subset(va.dist, exposure=="no prior exposure"), aes(x=x,y=0, color="homogeneous"))+ 
  geom_line(size=1)+
  facet_wrap(~exposure, scales = "free")+
  scale_colour_manual(values=c("#7E1717"))+
  geom_vline(aes(xintercept = 0.49), color = "#D1DBBD", size=1, linetype="dashed")+
  geom_vline(aes(xintercept = 1), color = "#D1DBBD", size=1, linetype="dashed")+
  geom_vline(aes(xintercept = par1.hom), color = "#7E1717", size=3)+
  coord_cartesian(ylim=c(0,3))+
  ylab(expression("f(p)"))+
  xlab(expression("p, infection probability"))+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,size=15),
        panel.grid = element_blank(), axis.line=element_line(),
        legend.position="top",legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))
rhom

ggsave(file="Fig/homogeneous_fits_dose-response.png",width=5,height=6,units="in",limitsize=FALSE)

##4.0 MODEL COMPARISONS

model.comp = store2 %>%
  pivot_longer(
    cols = -group)

distr = colsplit(model.comp$name, "\\.", c("extra_name","dist"))
model.comp$dist = distr$dist
model.comp$extra_name = distr$extra_name
model.comp$npars = 2
model.comp$npars[model.comp$dist=="hom"]=1

len.dose = length(unique(z3$secondary_dose))
no.birds = 12
no.tot.birds = 72
model.comp$sort.col = 1
model.comp$sort.col[model.comp$dist=="gamma"] = 3

model.comp.gamma.hom = model.comp %>%
  filter(extra_name=="dev") %>%
  drop_na()%>%
  group_by(group)%>%
  arrange(sort.col)%>%
  mutate(diff.in.lik = lag(value)-value) %>% 
  mutate (diff.in.pars = npars-lag(npars))%>%
  mutate(LRT.p.value = dchisq(diff.in.lik, diff.in.pars))%>%
  filter(group=="0 va"|group=="750 va"|group=="30000 va")%>%
  mutate(LRT.p.value = round(LRT.p.value, 5))
model.comp.gamma.hom$LRT.p.value[model.comp.gamma.hom$diff.in.lik<0]="not better"
View(model.comp.gamma.hom)

#Gamma Coefficient of variation = std/mean
#no 0.8992116 
#low 1.6298956 
#high 2.5106963

#Gamma Variance parameter = 1/shape parameter of gamma (from Ben-Ami)
#no 0.8085816 
#low  2.6565598
#high 6.3035959

