###Logistic Regression
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

n_eye_scores <- master %>%
  filter(l_eye_score != "na")
n_eye_scores

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


# a.ab <- m.ab %>%
#   dplyr::select(band_number, dpi, quantity, tes, threshold_cutoff, pathology_cutoff, primary_treatment, secondary_dose, experiment_notes)

#Rules for infection: We defined a bird as infected if it had eye score > 0.5, pathogen load > 50 copies, or both
m.ab <- m.ab %>%
  mutate(inf = ifelse((coalesce(quantity, 0) > threshold_cutoff | coalesce(tes, 0) > 0.5), 1, 0))%>%
  ungroup()

a.ab <- m.ab %>%
  mutate(inf_quant = ifelse((coalesce(quantity, 0) > threshold_cutoff), 1, 0))%>%
  ungroup()

a.ab <- a.ab %>%
  mutate(inf_final = ifelse((coalesce(quantity, 0) > threshold_cutoff | coalesce(tes, 0) > 0.5), 1, 0))%>%
  ungroup()

differing_band_numbers <- a.ab %>%
  filter(inf != inf_quant & primary_treatment == "Sham") %>%
  select(band_number, inf_quant, inf_final, dpi, primary_treatment, secondary_dose, tes, quantity, experiment_notes)

# Display the differing band_numbers
print(differing_band_numbers)
#including tes > 0.5 includes 8 birds during secondary that would not be included with path load alone

ggplot(a.ab %>% filter(inf_final ==1), aes(x=dpi, y=log10(quantity+1), color=fct_rev(primary_treatment)))+
  geom_point()+
  geom_line(aes(groups=band_number))+
  geom_hline(yintercept=log10(50))

#2419, 2434 0.5 on dpi 21 sham primary treatment
inf.sham <- m.ab %>% filter(band_number %in% c(2419, 2434, 2569))%>%
  dplyr::select(band_number, dpi, quantity, tes, elisa_od)


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
  geom_hline(yintercept=log10(50), color="red")+
  geom_hline(yintercept=0.5, aes(color="black"))+
  labs(x="DPI", y= "TES/Log10(quantity1)")+
  facet_wrap(~inf)

#eye score and path load
ggplot(m.ab, aes(x=tes, y=log10(quantity1), color=fct_rev(primary_treatment)))+
  geom_point(height=0)+
  facet_wrap(~inf_sec~dpi)

a<-ggplot(m.ab, aes(x=dpi, y=elisa_od, color=primary_treatment))+
  geom_point()+
  geom_hline(yintercept=0.061)
b<-ggplot(m.ab %>% filter(band_number %in% c(2452, 2562)), aes(x=dpi, y=quantity1, color=band_number))+
  geom_point()
a+b
####Omit Birds####
###PROBLEM BIRDS; n= 10###

prob_birds <- m.ab %>% filter(band_number %in% c(2274, 2514, 2469, 2520, 2494, 2452, 2562, 2419, 2434, 2505))%>%
  select(band_number, dpi, primary_treatment, secondary_dose, quantity, tes, l_eye_score, r_eye_score, elisa_od)

prob_birds %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(by=dpi)

omit_primary <- m.ab %>% filter(band_number %in% c(2505, 2452, 2562, 2419, 2434, 2514))

omit_primary %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(by=dpi)

#Sham = 51, Low = 52, High = 53
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#omit 2505 from data set for analysis as it was positive at quarantine (n=1)
m.ab <- m.ab %>%
  filter(band_number != 2505)

#omit 2452 and 2562 bc positive qpcr and they are shams (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2452, 2562)) 

#Sham = 49, Low = 51, High = 53
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#omit sham birds with eye scores (0.5) (2419, 2434) (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2419, 2434))

#Sham = 47, Low = 51, High = 53
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#omit sham bird with path load (2514) (n=1)
m.ab <- m.ab %>% filter(band_number != 2514)

#Sham = 46, Low = 51, High = 53
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#Total omitted from primary = 7; all from Sham treatment

#check birds that were not recovered by secondary infection: 2274, 2469, 2520, 2494 (2514 not recovered but omitted above)
not_recovered <- m.ab %>% 
  filter(band_number %in% c(2274, 2469, 2520, 2494))%>%
  dplyr::select(band_number, dpi, primary_treatment, secondary_dose, quantity, tes)

not_recovered %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

#confirm quantity > cutoff at dpi 41
ggplot(not_recovered, aes(x=dpi, y=quantity, color=band_number))+
  geom_point()+
  geom_path()+
  facet_wrap(~band_number)+
  ylim(min=0, max=700)


#omit from dataset - were still infected (path load) on day 41 prior to reinfection day 42
m.ab <- m.ab %>%
  filter(!(band_number %in% c(2274, 2469, 2520, 2494)))

#Final Secondary sample sizes
m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )

m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary()
m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  )
master %>% filter(primary_treatment == "sham" & secondary_dose == 0)
####Do Antibody Levels Predict Susceptibility?####
####set theme for susceptibility graphs####
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

####Wide Data -8####
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

####Model -8####
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
wibird$dpi <- "DPPI -8"
dpi_namespre <- c(
  "DPPI -8" = "DPPI -8"
)
sus_pre <- ggplot(wibird, aes(x=(scaled_elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi ~log10.sec_dose, labeller = as_labeller(c(dpi_namespre, sec_dose_names)))#+
sus_theme_pre
sus_pre

sus_pre_comb <- ggplot(wibird, aes(x=(scaled_elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_jitter(alpha=0.5, size=2, width=0, height=0.01)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~dpi, labeller=as_labeller(c(dpi_namespre, sec_dose_names)), strip.position = "right")

sus_pre_comb
#cumulative not scaled


# not scaled
mod <- glm.ab.pre
dat.new=expand.grid(log10.sec_dose=unique(wibird$log10.sec_dose),
                    inf_sec=unique(wibird$inf_sec),
                    elisa_od_pre = unique(wibird$elisa_od_pre))#new grid to put predictions into
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

wibird <- wibird %>%
  group_by(log10.sec_dose)%>%
  arrange(elisa_od_pre) %>%  # Sort by elisa_od_pre (or another predictor)
  mutate(cumulative_points = cumsum(inf_sec) / seq_along(inf_sec))  # Cumulative sum of inf_sec

#names
wibird$dpi <- "DPPI -8"
dpi_namespre <- c(
  "DPPI -8" = "DPPI -8"
)

sus_pre_raw<-ggplot(wibird, aes(x=(elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_point(data = wibird, aes(x = elisa_od_pre, y = cumulative_points), size = 1) +  # Cumulative points
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi ~log10.sec_dose, labeller = as_labeller(c(dpi_namespre, sec_dose_names)))
sus_pre_raw
#ggsave('sus_ab_pre.jpeg', width=6, height=4, units='in')

####Wide Data 14####
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

####Model 14####
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

wibird14$dpi <- "DPPI 14"
dpi_names14 <- c(
  "DPPI 14" = "DPPI 14"
)
#scaled
sus_14 <- ggplot(wibird14, aes(x=(scaled_elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Infection Status Upon Secondary Challenge (0|1)", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names14,sec_dose_names)))#+
sus_theme_14
sus_14

sus_14_comb <- ggplot(wibird14, aes(x=(scaled_elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Susceptibility (Infected 0|1)", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~dpi, labeller = as_labeller(dpi_names14), strip.position = "right")

sus_14_comb
#not scaled
mod <- glm.ab14
dat.new=expand.grid(log10.sec_dose=unique(wibird14$log10.sec_dose),
                    inf_sec=unique(wibird14$inf_sec),
                    elisa_od_14 = unique(wibird14$elisa_od_14))
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

wibird14$dpi <- "DPPI 14"
dpi_names14 <- c(
  "DPPI 14" = "DPPI 14"
)

wibird14 <- wibird14 %>%
  group_by(log10.sec_dose)%>%
  arrange(elisa_od_14) %>%  # Sort by elisa_od_pre (or another predictor)
  mutate(cumulative_points = cumsum(inf_sec) / seq_along(inf_sec))  # Cumulative sum of inf_sec

sus_14_raw<-ggplot(wibird14, aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_point(alpha=0.5, size=2)+
  geom_point(data = wibird14, aes(x = elisa_od_14, y = cumulative_points), size = 1) +  # Cumulative points
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Infection Status Upon Secondary Challenge (0|1)", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names14,sec_dose_names)))
sus_14_raw

####Wide Data 41####
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

####Model 41####
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

wibird41$dpi <- "DPPI 41"
dpi_names41 <- c(
  "DPPI 41" = "DPPI 41"
)
#scaled
sus_41 <- ggplot(wibird41, aes(x=(scaled_elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Scaled Antibody Levels", y= " ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names41,sec_dose_names)))#+
sus_theme_41
sus_41

sus_41_comb <- ggplot(wibird41, aes(x=(scaled_elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(scaled_elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= scaled_elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) + 
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Scaled Antibody Levels", y= " ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~dpi, labeller = as_labeller(dpi_names41), strip.position = "right")

sus_41_comb
#not scaled
mod <- glm.ab41
dat.new=expand.grid(log10.sec_dose=unique(wibird41$log10.sec_dose),
                    inf_sec=unique(wibird41$inf_sec),
                    #elisa_od_14_new=unique(wibird41$elisa_od_14_new),
                    elisa_od_41 = unique(wibird41$elisa_od_41))
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

wibird41$dpi <- "DPPI 41"
dpi_names41 <- c(
  "DPPI 41" = "DPPI 41"
)

wibird41 <- wibird41 %>%
  group_by(log10.sec_dose)%>%
  arrange(elisa_od_41) %>%  # Sort by elisa_od_pre (or another predictor)
  mutate(cumulative_points = cumsum(inf_sec) / seq_along(inf_sec))  # Cumulative sum of inf_sec

sus_41_raw<-ggplot(wibird41, aes(x=(elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_point(data = wibird41, aes(x = elisa_od_41, y = cumulative_points), size = 1) +  # Cumulative points
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Antibody Levels", y= " ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names41,sec_dose_names)))
sus_41_raw
#ggsave('sus_ab_41.jpeg', width=6, height=4, units='in')

####Combine Plots####
#FIGURE 3
library(cowplot)
p1 <- sus_pre + sus_theme_pre
p2 <- sus_14 + sus_theme_14
p3 <- sus_41 + sus_theme_41
legend <-get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

figure3 <- arrangeGrob(
  legend,
  arrangeGrob(p1, p2, p3, ncol = 1),
  ncol = 1,
  heights = c(1, 10)
)

library("grid")
grid.newpage()
grid.draw(figure3)

#combined
p1 <- sus_pre_comb + sus_theme_pre
p2 <- sus_14_comb + sus_theme_14
p3 <- sus_41_comb + sus_theme_41
legend <-get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

figure3_comb <- arrangeGrob(
  legend,
  arrangeGrob(p1, p2, p3, ncol = 1),
  ncol = 1,
  heights = c(1, 10)
)

grid.newpage()
grid.draw(figure3_comb)

fig3 <- grid.arrange(sus_pre, sus_14, sus_41, ncol=1)
g <- arrangeGrob(sus_pre, sus_14, sus_41, ncol=1) #generates g
print(g)
#ggsave(file="sus_all.jpeg", g) #saves g

#Raw - not scaled
r1 <- sus_pre_raw + sus_theme_pre
r2 <- sus_14_raw + sus_theme_14
r3 <- sus_41_raw + sus_theme_41
legend <-get_legend(r1)
r1 <- r1 + theme(legend.position = "none")

Raw_all <- arrangeGrob(
  legend,
  arrangeGrob(r1, r2, r3, ncol = 1),
  ncol = 1,
  heights = c(1, 10)
)

grid.newpage()
grid.draw(Raw_all)

#the magnitude of the immune response (IgY antibodies) predicts susceptibility upon secondary infection
max.ab <- m.ab %>%
  filter(dpi <42) %>%
  group_by(band_number, primary_treatment, secondary_dose)%>%
  reframe(max_ab = max(elisa_od, na.rm=TRUE),
            inf_sec = inf_sec,
          dpi = dpi)
max.ab <- max.ab %>%
  filter(dpi == 41)
max.ab$log10.sec_dose <-log10(max.ab$secondary_dose+1)


glm.abmax <- glm(inf_sec ~ max_ab+log10.sec_dose*primary_treatment, data=max.ab, family=binomial())
glm.abmax.sec <- glm(inf_sec ~ max_ab+log10.sec_dose, data=max.ab, family=binomial())
summary(glm.abmax)
plot(allEffects(glm.abmax))
simulateResiduals(glm.abmax, plot=T)

car::Anova(glm.abmax, type="III")
#birds with higher antibody levels, regardless of primary infection (mostly), were less susceptible to reinfection.

glm.inf <- glm(inf_sec ~ max_ab, data=max.ab, family=binomial())

#model selection
g0 <- glm(inf_sec ~ max_ab+primary_treatment, data=max.ab, family=binomial())
g1 <- glm(inf_sec ~ max_ab+log10.sec_dose, data=max.ab, family=binomial())
g2 <- glm(inf_sec ~ max_ab+log10.sec_dose + primary_treatment, data=max.ab, family=binomial())
g3 <- glm(inf_sec ~ max_ab+log10.sec_dose*primary_treatment, data=max.ab, family=binomial())
g4 <- glm(inf_sec ~ 1, data=max.ab, family=binomial())
g5 <- glm(inf_sec ~ log2(max_ab)+log10.sec_dose, data=max.ab, family=binomial())
g6 <- glm(inf_sec ~ log10(max_ab)+log10.sec_dose, data=max.ab, family=binomial())
g7 <- glm(inf_sec ~ log10(max_ab)+log10.sec_dose*primary_treatment, data=max.ab, family=binomial())
aictab(cand.set=list(g0, g1, g2, g3, g4, g5, g6, g7), modnames=c("g0", "g1", "g2", "g3", "g4", "g5", "g6", "g7"))

drop1(g1)
drop1(g2)
drop1(g3)
drop1(g5)

#Model g1 is the most supported (lowest AIC)
library(car)
#check for multicollinearity - VIF above 5 = multicollinearity
vif(g1)
vif(g5)
par(mfrow=c(2,2))
plot(g5)
plot(g1)
simulateResiduals(g5, plot=T)
hist(resid(g5))

mod <- g5
dat.new=expand.grid(log10.sec_dose=unique(max.ab$log10.sec_dose),
                    inf_sec=unique(max.ab$inf_sec),
                    primary_treatment=unique(max.ab$primary_treatment),
                    max_ab = unique(max.ab$max_ab))

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

max.ab <- max.ab%>%
  group_by(log10.sec_dose)%>%
  arrange(max_ab) %>%  # Sort by elisa_od_pre (or another predictor)
  mutate(cumulative_points = cumsum(inf_sec) / seq_along(inf_sec))  # Cumulative sum of inf_sec

ggplot(max.ab, aes(x=log10(max_ab), y=inf_sec, color=as.factor(log10.sec_dose)))+
  geom_jitter(height=0.0, alpha=0.5)+
  geom_point(aes(x=log10(max_ab), y=cumulative_points))+
  geom_line(data=dat.new, aes(x=log10(max_ab), y=yhat, group=log10.sec_dose))+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= log10(max_ab), y=yhat, group=log10.sec_dose, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Max Antibody Level", y="Probability of Reinfection", color="Log10(Secondary Dose+1)", fill="Log10(Secondary Dose+1)")+
  facet_grid(~log10.sec_dose)
  


#does the magnitude of eye score predict susceptibility?
max.tes <- m.ab %>%
  filter(dpi <42) %>%
  group_by(band_number, primary_treatment, secondary_dose)%>%
  reframe(max_tes = max(tes, na.rm=TRUE),
          inf_sec = inf_sec,
          dpi = dpi)
max.tes <- max.tes %>%
  filter(dpi == 41)
max.tes$log10.sec_dose <-log10(max.ab$secondary_dose+1)

glm.tesmax <- glm(inf_sec ~ max_tes+log10.sec_dose*primary_treatment, data=max.tes, family=binomial())
summary(glm.tesmax)
plot(allEffects(glm.tesmax))
car::Anova(glm.tesmax, type="III")
simulateResiduals(glm.tesmax, plot=T)

mod <- glm.tesmax
dat.new=expand.grid(log10.sec_dose=unique(max.tes$log10.sec_dose),
                    inf_sec=unique(max.tes$inf_sec),
                    primary_treatment=unique(max.tes$primary_treatment),
                    max_tes = unique(max.tes$max_tes))
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

max.tes <- max.tes%>%
  group_by(log10.sec_dose, primary_treatment)%>%
  arrange(max_tes) %>%  # Sort by elisa_od_pre (or another predictor)
  mutate(cumulative_points = cumsum(inf_sec) / seq_along(inf_sec))  # Cumulative sum of inf_sec

ggplot(max.tes, aes(x=max_tes, y=inf_sec, color=as.factor(log10.sec_dose)))+
  geom_jitter(height=0.01, width=0, alpha=0.5)+
  geom_point(aes(x=max_tes, y=cumulative_points), size=2)+
  geom_line(data=dat.new, aes(x=max_tes, y=yhat, group=log10.sec_dose))+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= max_tes, y=yhat, group=log10.sec_dose, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1)+
  scale_color_manual(values=c(sec_colors))+
  labs(x="Max Eye Score", y="Probability of Reinfection", color="Log10(Secondary Dose+1)", fill="Log10(Secondary Dose+1)")+
  facet_grid(~primary_treatment~log10.sec_dose)
   