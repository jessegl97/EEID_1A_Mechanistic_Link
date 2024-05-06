###########################
#Clean Code
#EEID 1a Antibody Analysis
#Jesse Garrett-Larsen
#6May24
###########################
  
#load packages
rm(list=ls())
library(tidyverse)
library(dplyr)
library(glmmTMB)
library(effects)
library(AICcmodavg)
library(emmeans)
library(MASS)
library(DHARMa)
library(lme4)
#install.packages("remotes")
#remotes::install_github("T-Engel/CValternatives")
library(CValternatives)
library(patchwork)
library(gridExtra)

source("dataCleaning_EEID1A.R")

####Format####
#set colors for primary treatment
pri_colors <- c("#D95F02", "#7570B3", "#1B9E77")
#set color scheme secondary treatment
sec_colors <- c("#8C754B", "#77AB59", "#59A5D8", "#9F77D9", "#FA8072")
#set theme
theme_set(theme_minimal())


#sample sizes
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#threshold cutoffs
m.ab$threshold_cutoff = 50
m.ab$seropos_cutoff = 0.061
m.ab$pathology_cutoff = 0

####Infection Rules####
#generate infected column - if copy number > 50, infected.
#for each timepoint only - can become infected or uninfected at next timepoint
m.ab <- m.ab %>%
  mutate(infected = ifelse(quantity>threshold_cutoff, 1, 0))

#generate infection data; if path load > 50 copies = infected
#1 = bird was successfully infected at any point during the priming phase of the experiment
#coalesce function replaces any NA with 0 here. 
#To be conservative, I am always assuming there is no path load unless it's measured with qPCR
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(infected_prim = ifelse(dpi < 42 & any(coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#generate infection data; if path load > 50 copies any time after secondary infection -> 1
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate (infected_sec = ifelse(dpi > 42 & any(coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#generate seropositivity data; if elisa OD > 0.061 = seropositive
#can become seropos or not at any point
m.ab <- m.ab %>%
  mutate(seropos = ifelse(elisa_od>seropos_cutoff, 1, 0))

#Seropos at any time during the priming phase
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(seropos_prim = ifelse(dpi < 42 & any(coalesce(seropos, 0) == 1), 1, 0)) %>%
  ungroup()

#Seropos at any time during the secondary phase
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(seropos_sec = ifelse(dpi > 42 & any(coalesce(seropos, 0) == 1), 1, 0)) %>%
  ungroup()

#generate eye score data; if eye score > 0 = diseased
#1 = bird had an eye score at that time point
m.ab <- m.ab %>%
  mutate(diseased = ifelse(tes>pathology_cutoff, 1, 0))

#ever diseased during primary infection
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(diseased_prim = ifelse(dpi < 42 & any(coalesce(diseased, 0) == 1), 1, 0)) %>%
  ungroup()

#ever diseased during secondary infection
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(diseased_sec = ifelse(dpi > 42 & any(coalesce(diseased, 0) == 1), 1, 0)) %>%
  ungroup()

#Rules for infection: We defined a bird as infected if it had eye score > 0, pathogen load > 50 copies, or both
m.ab <- m.ab %>%
  mutate(inf = ifelse((coalesce(quantity, 0) > threshold_cutoff | coalesce(tes, 0) > pathology_cutoff), 1, 0))%>%
  ungroup()


#ever inf during primary infection
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(inf_prim = ifelse(dpi < 42 & any(coalesce(inf, 0) == 1), 1, 0)) %>%
  ungroup()

#ever inf during secondary infection
m.ab <- m.ab %>%
  group_by(band_number) %>%
  mutate(
    # Check if any infection occurred after dpi 42
    inf_after_dpi42 = any(dpi > 42 & inf == 1),
    # Propagate the inf_after_dpi42 value to all instances within the group
    inf_sec = ifelse(inf_after_dpi42, 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-inf_after_dpi42)  # Remove the intermediate column

#Rules for infection: We defined a bird as infected if it had eye score > 0, pathogen load > 50 copies, or both
ggplot(m.ab, aes(x=dpi, y=tes, color=primary_treatment))+
  geom_jitter(height=0)+
  geom_jitter(data=m.ab, aes(x=dpi, y=log10(quantity1)), shape=17)+
  geom_hline(yintercept=log10(50))+
  labs(x="DPI", y= "TES/Log10(quantity1)")+
  facet_wrap(~inf)

#eye score and path load
ggplot(m.ab, aes(x=tes, y=log10(quantity1), color=primary_treatment))+
  geom_point(height=0)+
  facet_wrap(~inf_sec)

####Omit Birds####
#omit 2505 from data set for analysis as it was positive at quarantine
m.ab <- m.ab %>%
  filter(band_number != 2505)

#check birds that were not recovered by secondary infection: 2274, 2514, 2469, 2520
not_recovered <- m.ab %>% 
  filter(band_number %in% c(2274, 2514, 2469, 2520, 2494))

#confirm quantity > cutoff at dpi 41
ggplot(not_recovered, aes(x=dpi, y=quantity, color=band_number))+
  geom_point()+
  geom_path()+
  facet_wrap(~band_number)+
  ylim(min=0, max=700)


####Primary Challenge####
#data frame with only primary infection

p.ab <- m.ab %>%
  filter(dpi <=41)

#Baseline antibody levels did not differ
lm0 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma())
summary(lm0)
plot(allEffects(lm0))

emm_r <- emmeans(lm0, ~primary_treatment)
pairs(emm_r)

####Antibody Levels Across Primary Infection####
#I treat days post inoculation as a factor in these models and use band_number as a random effect
p.abt <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41))
p.abt$dpi.f <- as.factor(p.abt$dpi)
hist(p.abt$elisa_od)

#model comparison 
p1 <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2 <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2.5 <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), data=p.abt, family=Gamma(log))
p3<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p4<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p5 <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))

aictab(cand.set=list(p1, p2, p2.5, p3, p4, p5), modnames=c("p1", "p2", "p2.5", "p3", "p4", "p5"))

#p2 best model
lm1 <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log)) 
#fully interactive does not converge with inverse link- only log link
summary(lm1)
car::Anova(lm1, type = "III")
simulateResiduals(lm1, plot=T) #but residuals look terrible with log link

emm_results <- emmeans(lm1, ~ primary_treatment*dpi.f, scale="response")
pairs(emm_results)

#sex not significant
summary(p2.5)

#specify model to predict
mod.pred <- lm1
#Model predictions
dat.new=expand.grid(primary_treatment=unique(p.abt$primary_treatment),
                    elisa_od = unique(p.abt$elisa_od),
                    dpi.f = unique(p.abt$dpi.f))
#generate yhat values
dat.new$yhat = predict(mod.pred, type="response", newdata=dat.new, re.form=NA)
#generate SEM
preds = predict(mod.pred, type="link", newdata=dat.new, se.fit=TRUE, re.form=NA)
dat.new = cbind(dat.new, preds)

ilink <- family(mod.pred)$linkinv
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))
head(dat.new)


####Variability in Antibodies####
#add small constant to total eye score for calculation of variability
p.ab$tes.new <- p.ab$tes + .01
#new df with only individuals that were infected any time during primary challenge
p.ab$dpi.f <- as.factor(p.ab$dpi)
p.abi <- p.ab %>%
  filter(inf_prim == 1)%>%
  dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1)

#antibody variability calculations
a.cv <- p.abi %>% 
  group_by(dpi.f, primary_treatment)%>%
  drop_na(elisa_od)%>%
  dplyr::summarise(mean_od = mean(elisa_od),
                   bird_cv = sd(elisa_od)/mean(elisa_od),
                   bird_sd = sd(elisa_od),
                   bird_pv = PV(elisa_od),
                   elisa_od = elisa_od)
#add to df
p.abi.m <- left_join(p.abi, a.cv)

#overlay PV and elisa_od
#FIGURE 2
ggplot(data=p.abt, aes(x=dpi.f, y= elisa_od,
                       color=fct_rev(primary_treatment)))+
  geom_hline(yintercept = 0.061, linetype="dashed", alpha=0.25)+
  geom_jitter(size=2, alpha=0.5, width=0.2, height=0)+
   geom_path(data=dat.new, aes(x=dpi.f, y=yhat, group=primary_treatment, color=fct_rev(primary_treatment)))+
  geom_errorbar(data=dat.new, aes(ymin=Lower, ymax = Upper, x=dpi.f, y=yhat), color="black", width=0.1)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat, color=fct_rev(primary_treatment)), size=2, alpha=0.75, shape =16)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat), color="black", size=2, alpha=1, shape=1, stroke=0.1)+
  geom_point(data=a.cv, aes(x=dpi.f, y=bird_pv, fill=fct_rev(primary_treatment)), color="black", shape=25, size=3)+
  geom_path(data=a.cv, aes(x=dpi.f, y=bird_pv, group=primary_treatment, color=fct_rev(primary_treatment)),
            linetype="dotdash", size=0.5, alpha=0.5)+
  labs(y="Antibody Levels", 
       x= "Days Post Infection", shape="Primary Treatment", color="Primary Treatment", fill="Primary Treatment")+
  scale_color_manual(values=c(pri_colors))+
  scale_fill_manual(values=c(pri_colors))+
  theme_minimal()+
  scale_y_continuous(
    name= "Antibody Levels",
    sec.axi=sec_axis(~.*1, name="PV")
  )+
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "black", size=11, face="bold"),
    axis.text.y = element_text(color="black", face="bold"),
    axis.text.y.right = element_text(color="black", face="bold"),
    axis.title.y.right = element_text(color = "black", size=11, face = "bold")
  ) 

#ggsave('ab_all_prim.jpeg', width=6, height=4, units='in')

#eye score
e.cv <- p.abi %>% 
  group_by(dpi, primary_treatment)%>%
  summarise(bird_cv = sd(tes.new)/mean(tes.new),
            bird_sd = sd(tes.new),
            bird_pv = PV(tes.new))
#path load
q.cv <- p.abi %>% 
  group_by(dpi, primary_treatment)%>%
  drop_na(quantity1)%>%
  summarise(bird_cv = sd(quantity1)/mean(quantity1),
            bird_sd = sd(quantity1),
            bird_pv = PV(quantity1))
#antibodies
a.cv <- p.abi %>% 
  group_by(dpi, primary_treatment)%>%
  drop_na(elisa_od)%>%
  dplyr::summarise(mean_od = mean(elisa_od),
                   bird_cv = sd(elisa_od)/mean(elisa_od),
                   bird_sd = sd(elisa_od),
                   bird_pv = PV(elisa_od),
                   elisa_od = elisa_od)
#path load
qpv <- ggplot(q.cv, aes(x=dpi, y=bird_pv))+
  geom_point(aes(color=fct_rev(primary_treatment)))+
  geom_path(aes(group=primary_treatment, color=fct_rev(primary_treatment)))
qpv

#antibody pv
apv <- ggplot(a.cv, aes(x=dpi, y=bird_pv))+
  geom_point(aes(color=fct_rev(primary_treatment)))+
  geom_path(aes(group=primary_treatment, color=fct_rev(primary_treatment)))
apv

#antibody cv
acv <- ggplot(a.cv, aes(x=dpi, y=bird_cv))+
  geom_point(aes(color=fct_rev(primary_treatment)))+
  geom_path(aes(group=primary_treatment, color=fct_rev(primary_treatment)))
acv

#antibody pv
ggplot(a.cv , aes(x=dpi, y=bird_pv))+
  geom_point(aes(color=fct_rev(primary_treatment), alpha="PV"), size=3)+
  geom_line(aes(color=fct_rev(primary_treatment)))+
  geom_point(aes(x=dpi, y=bird_cv, color=fct_rev(primary_treatment), alpha="CV"), size=3)+
  geom_path(aes(dpi, y=bird_cv, color=fct_rev(primary_treatment)), linetype="dashed")+
  scale_color_manual(values = c(pri_colors))+
  scale_alpha_manual(values = c(PV=1, CV=0.1))+
  labs(x="Days Post Primary Inoculation", y="Variability", color="Primary Treatment", alpha = "Variability Metric")


#ggsave('ab_var_prim.jpeg', width=6, height=4, units='in')

####Remove Infected Birds for Secondary Analysis####

#omit from dataset - were still infected (path load) on day 41 prior to reinfection day 42
m.ab <- m.ab %>%
  filter(!(band_number %in% c(2274, 2514, 2469, 2520, 2494)))

m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  tbl_summary(
    by=secondary_dose
  )

m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  tbl_summary(
    by=inf_sec
  )

####Do Antibody Levels Predict Susceptibility?####
#set theme for susceptibility graphs
sus_theme_pre <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_blank(),  # Remove the x-axis title
    axis.title.y = element_text(size = 15, face="bold"),
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    strip.text = element_text(size = 10)  # Adjust the size of facet labels
  )

sus_theme_14 <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1), 
  axis.title.y = element_text(size = 15, face="bold"),
  legend.position = "none",  # Remove the legend
  plot.title = element_blank(),  # Remove the plot title
  axis.title.x = element_blank(),  # Remove the x-axis title
  legend.text = element_blank(),  # Remove the legend text
  legend.title = element_blank(),  # Remove the legend title
  legend.background = element_blank(),  # Remove the legend background
  strip.text = element_text(size = 10)  # Adjust the size of facet labels
)

sus_theme_41 <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1), 
  axis.title.y = element_text(size = 15, face="bold"),
  legend.position = "none",  # Remove the legend
  plot.title = element_blank(),  # Remove the plot title
  axis.title.x = element_text(size=15, face="bold"), 
  legend.text = element_blank(),  # Remove the legend text
  legend.title = element_blank(),  # Remove the legend title
  legend.background = element_blank(),  # Remove the legend background
  strip.text = element_text(size = 10)  # Adjust the size of facet labels
)


m.ab.w <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41))

wibird <- m.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)

wibird$elisa_od_pre <- wibird$`elisa_od_-8`
#could try scaling antibody levels to help w/ effect sizes: value - minimum / max - min 0 to 1
min_val <- min(wibird$elisa_od_pre)
max_val <- max(wibird$elisa_od_pre)

# Scale elisa_od_14 to a range of 0 to 1
wibird$scaled_elisa_od_pre <- (wibird$elisa_od_pre - min_val) / (max_val - min_val)

wibird <- wibird %>% mutate(totalTES = sum())
wibird$sec_dose.n <- wibird$secondary_dose+1
wibird$log10.sec_dose <- round(log10(wibird$sec_dose.n), digits = 2)

wibird$log10.sec_dose <- as.numeric(wibird$log10.sec_dose)

wibird$secondary_dose_fct <- factor(wibird$secondary_dose, levels = c(0, 30, 100, 300, 7000))

sec_dose_names <- c(
  "0" = "0",
  "1.49" = "30",
  "2" = "100",
  "2.48" = "300",
  "3.85" = "7000"
)

glm.ab.pre <- glm(inf_sec ~ elisa_od_pre + log10.sec_dose, data=wibird, family=binomial())
summary(glm.ab.pre)
car::Anova(glm.ab.pre, type = 3)
simulateResiduals(glm.ab.pre, plot=T)

# Scale elisa_od_pre to a range of 0 to 1
wibird$scaled_elisa_od_pre <- (wibird$elisa_od_pre - min_val) / (max_val - min_val)

glm.ab.pre.scaled <- glm(inf_sec ~ scaled_elisa_od_pre + log10.sec_dose, data=wibird, family=binomial())

summary(glm.ab.pre.scaled)
car::Anova(glm.ab.pre.scaled, type=3)
mod <- glm.ab.pre.scaled
dat.new=expand.grid(log10.sec_dose=unique(wibird$log10.sec_dose),
                    inf_sec=unique(wibird$inf_sec),
                    scaled_elisa_od_pre = unique(wibird$scaled_elisa_od_pre))#new grid to put predictions into
dat.new$yhat=predict(mod, type="response", newdata = dat.new)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

#scaled
sus_pre <- ggplot(wibird, aes(x=(scaled_elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~log10.sec_dose, nrow=1, labeller = as_labeller(sec_dose_names))+
  sus_theme_pre
sus_pre
# not scaled
# ggplot(wibird, aes(x=(elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
#   geom_point(alpha=0.5)+
#   geom_line(data=dat.new, aes(x=(elisa_od_pre), y=yhat), alpha = 1, size =1)+
#   geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
#   scale_color_manual(values = c(sec_colors))+
#   scale_fill_manual(values = c(sec_colors))+
#   labs(x="Antibody Levels Day -8", y= "Susceptibility (Secondary Infection 0|1)", color="log10(Secondary Dose)", fill ="log10(Secondary Dose)")+
#   facet_wrap(~log10.sec_dose,labeller = as_labeller(sec_dose_names), nrow=1)+
#   scale_x_continuous(name = "Antibody Levels Day -8",
#                      labels = scales::number_format(accuracy = 0.01))+
#   sus_theme

#ggsave('sus_ab_pre.jpeg', width=6, height=4, units='in')

##Do Antibody Levels on Day 14 Predict Susceptibility?
m.ab.w <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41))

wibird <- m.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)
wibird14 <- wibird %>%
  drop_na(elisa_od_14)

#could try scaling antibody levels to help w/ effect sizes: value - minimum / max - min 0 to 1
min_val <- min(wibird14$elisa_od_14)
max_val <- max(wibird14$elisa_od_14)

# Scale elisa_od_14 to a range of 0 to 1
wibird14$scaled_elisa_od_14 <- (wibird14$elisa_od_14 - min_val) / (max_val - min_val)

wibird14 <- wibird14 %>% mutate(totalTES = sum())
wibird14$sec_dose.n <- wibird14$secondary_dose+1
wibird14$log10.sec_dose <- round(log10(wibird14$sec_dose.n), digits = 2)

wibird14$log10.sec_dose <- as.numeric(wibird14$log10.sec_dose)

#Are antibody levels on day 14 predictive of secondary infection?
glm.ab14 <- glm(inf_sec ~ elisa_od_14 + log10.sec_dose, data=wibird14, family=binomial())
glm.ab14.scaled <- glm(inf_sec ~ scaled_elisa_od_14 + log10.sec_dose, data=wibird14, family=binomial())


AIC(glm.ab14, glm.ab14.scaled)
summary(glm.ab14)
summary(glm.ab14.scaled)

car::Anova(glm.ab14, type = 3)
simulateResiduals(glm.ab14, plot=T)
car::Anova(glm.ab14.scaled, type = 3)
simulateResiduals(glm.ab14.scaled, plot=T)

mod <- glm.ab14.scaled
dat.new=expand.grid(log10.sec_dose=unique(wibird14$log10.sec_dose),
                    inf_sec=unique(wibird14$inf_sec),
                    scaled_elisa_od_14 = unique(wibird14$scaled_elisa_od_14))
                  # elisa_od_14 = unique(wibird14$elisa_od_14))
dat.new$yhat=predict(mod, type="response", newdata = dat.new)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))


#scaled
sus_14 <- ggplot(wibird14, aes(x=(scaled_elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Susceptibility (Infected 1|0)", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~log10.sec_dose, nrow=1, labeller = as_labeller(sec_dose_names))+
  sus_theme_14
sus_14


#ggsave('sus_ab_14.jpeg', width=6, height=4, units='in')

#not scaled
# ggplot(wibird, aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
#   geom_line(data=dat.new, aes(x=(elisa_od_14), y=yhat), alpha = 1, size =1)+
#   geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
#   geom_point(alpha=0.5, size=2)+
#   scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
#   scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
#   labs(x="Antibody Levels Day 14", y= "Susceptibility (Infected 1|0)", color="Secondary Dose", fill ="Secondary Dose")+
#   scale_x_continuous(name = "Antibody Levels Day 14",
#                      labels = scales::number_format(accuracy = 0.01))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         legend.position="top",
#         legend.background= element_rect(size=0.25, linetype="solid"))+
#   facet_wrap(~log10.sec_dose, nrow=1)

####Do Antibody Levels on Day 41 Predict Susceptibility?
wibird <- m.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)

wibird41 <- wibird %>%
  drop_na(elisa_od_41)

#scale
min_val <- min(wibird41$elisa_od_41)
max_val <- max(wibird41$elisa_od_41)
wibird41$scaled_elisa_od_41 <- (wibird41$elisa_od_41 - min_val) / (max_val - min_val)
wibird41 <- wibird41 %>% mutate(totalTES = sum())
wibird41$sec_dose.n <- wibird41$secondary_dose+1
wibird41$log10.sec_dose <- round(log10(wibird41$sec_dose.n), digits = 2)

wibird41$log10.sec_dose <- as.numeric(wibird41$log10.sec_dose)

glm.ab41 <- glm(inf_sec ~ elisa_od_41 + log10.sec_dose, data=wibird41, family=binomial())

summary(glm.ab41)
car::Anova(glm.ab41, type = 3)
simulateResiduals(glm.ab41, plot=T)


glm.ab41.scaled <- glm(inf_sec ~ scaled_elisa_od_41 + log10.sec_dose, data=wibird41, family=binomial())
summary(glm.ab41.scaled)

mod <- glm.ab41.scaled
dat.new=expand.grid(log10.sec_dose=unique(wibird41$log10.sec_dose),
                    inf_sec=unique(wibird41$inf_sec),
                    #elisa_od_14_new=unique(wibird41$elisa_od_14_new),
                    scaled_elisa_od_41 = unique(wibird41$scaled_elisa_od_41))
#elisa_od_41 = unique(wibird41$elisa_od_14))
dat.new$yhat=predict(mod, type="response", newdata = dat.new)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

#scaled
sus_41 <- ggplot(wibird41, aes(x=(scaled_elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Scaled Antibody Levels", y= " ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~log10.sec_dose, nrow=1, labeller = as_labeller(sec_dose_names))+
  sus_theme_41
sus_41
#not scaled
# ggplot(wibird41, aes(x=(elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
#   geom_point(alpha=0.5)+
#   geom_line(data=dat.new, aes(x=(elisa_od_41), y=yhat), alpha = 1, size =1)+
#   geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
#   scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
#   scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
#   labs(x="Antibody Levels Day 14", y= "Susceptibility (Infected 1|0)", color="Secondary Dose", fill ="Secondary Dose")+
#   scale_x_continuous(name = "Antibody Levels Day 41",
#                      labels = scales::number_format(accuracy = 0.01))+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#         legend.position="top",
#         legend.background= element_rect(size=0.25, linetype="solid"))+
#   facet_wrap(~log10.sec_dose, nrow=1, labeller = as_labeller(sec_dose_names))

#ggsave('sus_ab_41.jpeg', width=6, height=4, units='in')

#combine
#FIGURE 3
grid.arrange(sus_pre, sus_14, sus_41, ncol=1)


####Variability With Hawley et al., 2024####
dose <- c("Sham", "Low", "High")
CV <- as.numeric(c(0.899, 1.630, 2.511))
ab.pv <- as.numeric(c(0.04156784, 0.14716858, 0.18885516))
ab.cv <- as.numeric(c(0.03997555, 0.17202405, 0.26176249))
mean_sus <- as.numeric(c(1.181, 0.446, 0.192))
mean_ab <- as.numeric(c(0.045, 0.04786, 0.0579))
dpi <- 41



het.df <- data.frame(dose, CV, mean_sus, ab.pv, ab.cv, dpi, mean_ab)

het.df$mean_sus_inv <- 1/het.df$mean_sus

#hawley et al 2024
het.paper.only <- ggplot(het.df, aes(x=Metric))+
  geom_point(size=5,aes(x="CV", y=CV, color=dose), shape = 17)+
  geom_point(size=5,aes(x="CV", y=CV, shape = "CV"), color="black", shape=2, stroke = 1)+
  scale_color_manual(values=c(pri_colors))+
  scale_y_continuous(
    name= "Variability in Susceptibility Hawley et al. 2024"
  )+
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "black", size=12, face="bold"),
    legend.position = "none"
  ) +
  labs(color = "Primary Treatment", shape= "Variability Metric")
het.paper.only

ab.paper.only <- ggplot(het.df, x=Metric)+
  geom_point(size=5,aes(x="PV", y=ab.pv, color=dose, shape="PV"))+
  geom_point(size=5,aes(x="CV", y=ab.cv, color=dose, shape = "CV"))+
  #geom_point(size=3.6,aes(x="PV", y=ab.pv), color="black", shape=1, stroke = 2)+
  scale_color_manual(values=c(pri_colors))+
  scale_y_continuous(
    name="Variability in Antibody Levels"
  )+
  scale_x_discrete(
    name="Metric"
  )+
  scale_shape_manual(values=c(17,19))+
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "black", size=12, face = "bold"),
    legend.position="right",
    legend.background= element_rect(size=0.25, linetype="solid"))+
  labs(color = "Primary Treatment", shape= "Variability Metric")
ab.paper.only

#FIGURE 4
grid.arrange(het.paper.only, ab.paper.only, ncol=2)
