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

setwd("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link")
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


####Primary Challenge####
#Final Primary sample sizes for analysis
m.ab %>%
  dplyr::select(dpi, primary_treatment)%>%
  tbl_summary(
    by=dpi
  )

#data frame with only primary infection

p.ab <- m.ab %>%
  filter(dpi <=41)

#Baseline antibody levels did not differ
lm0 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0), family=Gamma())
summary(lm0)
plot(allEffects(lm0))

emm_r <- emmeans(lm0, ~primary_treatment)
pairs(emm_r)

####Does body condition / mass predict max antibody response? No.
p.ab$mass <- as.numeric(p.ab$mass)
max.ab <- p.ab %>%
  group_by(band_number, primary_treatment) %>%
  summarise(max_od = max(elisa_od, na.rm = TRUE),
            mass = as.numeric(mass[dpi == -8]),
            dpi = dpi,
            sex = sex)

lmm <- glm(max_od ~ mass + primary_treatment + sex, data=max.ab %>%filter(dpi == -8), family=Gamma())
summary(lmm)
plot(allEffects(lmm))

ggplot(max.ab%>%filter(dpi == -8), aes(x=mass, y=max_od, color=primary_treatment))+
  geom_boxplot(aes(x=mass, y=0.1), width=0.01, color="black")+
  geom_point()+
  geom_smooth(method = "glm", method.args = list(family="Gamma"))+
  facet_wrap(~sex~primary_treatment, nrow=2)

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

library(VGAM)
drop1(p1, test="Chisq")
drop1(p2, test="Chisq")
drop1(p2.5, test="Chisq")
drop1(p3, test="Chisq")
drop1(p4, test="Chisq")
drop1(p5, test="Chisq")

#drop1 function?
#p2 best model
p.abt$primary_treatment <- relevel(p.abt$primary_treatment, ref = "Sham")
lm1 <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log)) 

#fully interactive does not converge with inverse link- only log link
summary(lm1)
car::Anova(lm1, type = "III")
simulateResiduals(lm1, plot=T) #but residuals look terrible with log link
hist(resid(lm1))
simulation_output <- simulateResiduals(fittedModel = lm1, n = 1000)

# Plot residuals
plotSimulatedResiduals(simulation_output, rank = TRUE)



emm_results <- emmeans(lm1, ~ primary_treatment|dpi.f, scale="response")
pairwise <- pairs(emm_results, adjust="tukey")
(summary(pairwise))
pairwise_df <- as.data.frame(pairwise)

ggplot(pairwise_df, aes(x = contrast, y = estimate, color=dpi.f)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Pairwise Comparisons of Primary Treatments within Each dpi.f",
       x = "Comparison",
       y = "Estimated Difference") +
  coord_flip() +  # Flip the coordinates for better readability
  theme_minimal()

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
  dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number)

#antibody variability calculations
a.cv <- p.abi %>% 
  group_by(dpi.f, primary_treatment)%>%
  drop_na(elisa_od)%>%
  dplyr::reframe(mean_od = mean(elisa_od),
                   bird_cv = sd(elisa_od)/mean(elisa_od),
                   bird_sd = sd(elisa_od),
                   bird_pv = PV(elisa_od),
                   elisa_od = elisa_od,
                  band_number = band_number,
                 n_individuals = n_distinct(band_number))%>%
  ungroup()

#add to df
p.abi.m <- left_join(p.abi, a.cv, by=c("dpi.f", "band_number"))
p.abi.m <- p.abi.m %>%
  filter(dpi.f %in% c(-8, 14, 41))
summary_tibble <- a.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    mean_bird_cv = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    mean_bird_sd = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    mean_bird_pv = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    avg_elisa_od = mean(elisa_od, na.rm = TRUE),# Calculate the average elisa_od for each group
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
unique(p.abi.m$band_number)

##Bootstrapping - resample underlying data (elisa_od) for each group then recalculate bird_cv and bird_pv for each bootstrap sample
library(boot)
set.seed(123)

bootstrap_fn <- function(data, indices) {
  resampled_data <- data[indices, ]
  
  bird_cv <- sd(resampled_data$elisa_od, na.rm = TRUE) / mean(resampled_data$elisa_od, na.rm = TRUE)
  bird_pv <- PV(resampled_data$elisa_od)
  
  return(c(bird_cv, bird_pv))
}

#apply bootstrapping to each group
boot_results <- a.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  do({
    boot_out <- boot(data = ., statistic = bootstrap_fn, R = 1000)
    
    ci_cv <- boot.ci(boot_out, index=1, type= "perc")$percent[4:5]
    ci_pv <- boot.ci(boot_out, index=2, type= "perc")$percent[4:5]
    
    tibble(
      dpi.f = unique(.$dpi.f),
      primary_treatment = unique(.$primary_treatment),
      mean_bird_cv = mean(.$bird_cv, na.rm =TRUE),
      mean_bird_pv = mean(.$bird_pv, na.rm =TRUE),
      ci_lower_cv = ci_cv[1],
      ci_upper_cv = ci_cv[2],
      ci_lower_pv = ci_pv[1],
      ci_upper_pv = ci_pv[2]
    )
  }) %>%
  ungroup()
boot_results

#overlay PV and elisa_od
#FIGURE 2
ggplot(data=p.abt, aes(x=dpi.f, y= elisa_od,
                       color=fct_rev(primary_treatment)))+
  #geom_hline(yintercept = 0.061, linetype="dashed", alpha=0.25)+
  #geom_jitter(size=2, alpha=0.5, width=0.2, height=0)+
  geom_path(data=dat.new, aes(x=dpi.f, y=yhat, group=primary_treatment, color=fct_rev(primary_treatment)))+
  geom_errorbar(data=dat.new, aes(ymin=Lower, ymax = Upper, x=dpi.f, y=yhat), color="black", width=0.1)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat, color=fct_rev(primary_treatment)), size=2, alpha=0.75, shape =16)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat), color="black", size=2, alpha=1, shape=1, stroke=0.1)+
  #geom_errorbar(data = boot_results, aes(x = dpi.f, ymin = ci_lower_pv, ymax = ci_upper_pv, group = primary_treatment), 
  #              color = "black", width = 0.2) +  # Error bars for bird_pv confidence intervals
  #geom_point(data = boot_results, aes(x = dpi.f, y = mean_bird_pv, fill = fct_rev(primary_treatment)), 
  #           size = 3, shape = 25) +  # Bird_pv points
  #geom_path(data = boot_results, aes(x = dpi.f, y = mean_bird_pv, group = primary_treatment, color = fct_rev(primary_treatment)),
  #          size = 0.5, alpha = 0.5) +
 #CV
   # geom_point(data=a.cv, aes(x=dpi.f, y=bird_cv, fill=fct_rev(primary_treatment)), shape=15, size=3, alpha=0.1)+
  # geom_path(data=a.cv, aes(x=dpi.f, y=bird_cv, group=primary_treatment, color=fct_rev(primary_treatment)),
  #            linetype="dashed", size=0.5, alpha=0.25)+
  labs(y="Antibody Levels", 
       x= "Days Post Primary Inoculation", shape="Primary Treatment", color="Primary Treatment", fill="Primary Treatment")+
  # Consolidating scale_y_continuous() settings
  scale_y_continuous(
    #limits = c(0.03, 0.31),  # Limits for the main y-axis
    name = "Antibody Levels", 
  ) +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  theme_minimal() +
  
  # Theme adjustments
  theme(
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.text.y.right = element_text(color = "black", face = "bold", size = 15),
    axis.title.y.right = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )
  
ggplot(data=p.abt, aes(x=dpi.f, y=elisa_od, color=fct_rev(primary_treatment))) +
  geom_hline(yintercept = 0.061, linetype = "dashed", alpha = 0.25) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2, height = 0) +
  
  # Path and error bars for yhat data
  geom_path(data = dat.new, aes(x = dpi.f, y = yhat, group = primary_treatment, color = fct_rev(primary_treatment))) +
  geom_errorbar(data = dat.new, aes(ymin = Lower, ymax = Upper, x = dpi.f, y = yhat), color = "black", width = 0.1) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat, color = fct_rev(primary_treatment)), size = 2, alpha = 0.75, shape = 16) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat), color = "black", size = 2, alpha = 1, shape = 1, stroke = 0.1) +
  
  # Add bird_pv points and their confidence intervals (from boot_results)
  # geom_errorbar(data = boot_results, aes(x = dpi.f, y=mean_bird_pv, ymin = ci_lower_pv, ymax = ci_upper_pv, group = primary_treatment), 
  #               width = 0.1) +  # Error bars for bird_pv confidence intervals
  # geom_point(data = boot_results, aes(x = dpi.f, y = mean_bird_pv, fill = fct_rev(primary_treatment)), 
  #            size = 3, shape = 25) +  # Bird_pv points
  # 
  # geom_path(data = boot_results, aes(x = dpi.f, y = mean_bird_pv, group = primary_treatment, color = fct_rev(primary_treatment)),
  #           size = 0.5, alpha = 0.5, linetype="dashed") +
  
  geom_errorbar(data = boot_results, aes(x = dpi.f, y=mean_bird_cv, ymin = ci_lower_cv, ymax = ci_upper_cv, group = primary_treatment), 
                width = 0.05) +  # Error bars for bird_pv confidence intervals
  geom_point(data = boot_results, aes(x = dpi.f, y = mean_bird_cv, fill = fct_rev(primary_treatment)), 
             size = 3, shape = 25) +  # Bird_pv points
  
  geom_path(data = boot_results, aes(x = dpi.f, y = mean_bird_cv, group = primary_treatment, color = fct_rev(primary_treatment)),
            size = 0.5, alpha = 0.5, linetype="dashed") +  

  labs(y = "Antibody Levels", 
       x = "Days Post Primary Inoculation", 
       shape = "Primary Treatment", 
       color = "Primary Treatment", 
       fill = "Primary Treatment") +
  
  # Consolidating scale_y_continuous() settings
  scale_y_continuous(
    #limits = c(0.03, 0.31),  # Limits for the main y-axis
    name = "ELISA Optical Density [OD]", 
    sec.axis = sec_axis(~.*1, name = "PV")  # Secondary axis with a scaling factor of 1
  ) +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  theme_minimal() +
  
  # Theme adjustments
  theme(
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.text.y.right = element_text(color = "black", face = "bold", size = 15),
    axis.title.y.right = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

#ggsave('ab_all_prim.jpeg', width=6, height=4, units='in')
ggplot(boot_results, aes(x=dpi.f, y=mean_bird_cv, color=fct_rev(primary_treatment)))+
  geom_errorbar(data = boot_results, aes(x = dpi.f, y=mean_bird_cv, ymin = ci_lower_cv, ymax = ci_upper_cv, group = primary_treatment), 
              width = 0.05) +  # Error bars for bird_pv confidence intervals
  geom_point(data = boot_results, aes(x = dpi.f, y = mean_bird_cv, fill = fct_rev(primary_treatment)), 
             size = 3, shape = 25) +  # Bird_pv points
  geom_path(data = boot_results, aes(x = dpi.f, y = mean_bird_cv, group = primary_treatment, color = fct_rev(primary_treatment)),
            size = 0.5, alpha = 0.5, linetype="dashed") +  
  scale_color_manual(values=c(pri_colors))+
  scale_fill_manual(values=c(pri_colors))+
  theme_minimal()+
  labs(y = "Variability in Antibody Levels (PV)", 
       x = "Days Post Primary Inoculation", 
       shape = "Primary Treatment", 
       color = "Primary Treatment", 
       fill = "Primary Treatment")+
  theme(
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

#with sex
a.cv.s <- p.abi %>% 
  group_by(dpi.f, primary_treatment, sex)%>%
  drop_na(elisa_od)%>%
  dplyr::reframe(mean_od = mean(elisa_od),
                 bird_cv = sd(elisa_od)/mean(elisa_od),
                 bird_sd = sd(elisa_od),
                 bird_pv = PV(elisa_od),
                 elisa_od = elisa_od)


ggplot(data=p.abt, aes(x=dpi.f, y= elisa_od,
                       color=fct_rev(primary_treatment)))+
  geom_hline(yintercept = 0.061, linetype="dashed", alpha=0.25)+
  geom_jitter(size=2, alpha=0.5, width=0.2, height=0, aes(shape=sex))+
  geom_path(data=dat.new, aes(x=dpi.f, y=yhat, group=primary_treatment, color=fct_rev(primary_treatment)))+
  geom_errorbar(data=dat.new, aes(ymin=Lower, ymax = Upper, x=dpi.f, y=yhat), color="black", width=0.1)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat, color=fct_rev(primary_treatment)), size=2, alpha=0.75, shape =16)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat), color="black", size=2, alpha=1, shape=1, stroke=0.1)+
  geom_point(data=a.cv.s, aes(x=dpi.f, y=bird_pv, fill=fct_rev(primary_treatment), shape=sex), size=3)+
  geom_path(data=a.cv.s, aes(x=dpi.f, y=bird_pv, group=interaction(primary_treatment,sex), linetype=(interaction(primary_treatment, sex)), color=fct_rev(primary_treatment)),
            size=0.5, alpha=0.5)+
  #CV
  # geom_point(data=a.cv.s, aes(x=dpi.f, y=bird_cv, fill=fct_rev(primary_treatment)), shape=15, size=3, alpha=0.1)+
  # geom_path(data=a.cv.s, aes(x=dpi.f, y=bird_cv, group=primary_treatment, color=fct_rev(primary_treatment)),
  #            linetype="dashed", size=0.5, alpha=0.25)+
  labs(y="Antibody Levels",
       x= "Days Post Primary Inoculation", shape="Primary Treatment", color="Primary Treatment", fill="Primary Treatment")+
  scale_color_manual(values=c(pri_colors))+
  scale_fill_manual(values=c(pri_colors))+
  theme_minimal()+
  scale_y_continuous(
    name= "ELISA Optical Density [OD]",
    sec.axi=sec_axis(~.*1, name="PV")
  )+
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "black", size = 11, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.y.right = element_text(color = "black", face = "bold"),
    axis.title.y.right = element_text(color = "black", size = 11, face = "bold"),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )+
  facet_wrap(~sex)

m.ab$tes
m.ab14 <- m.ab %>% filter(dpi == 14)
is.na(m.ab14$tes)

hmm <- m.ab14 %>%
  group_by(primary_treatment) %>%
  summarise(bird_cv = sd(tes)/mean(tes),
            bird_pv = PV(tes),
            bird_pv_plus = PV(tes+0.0001))

hmm

1.78/1.09
0.515/0.674
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
###Issue is that a.cv has lots of occurrences so alpha doesn't matter - need to make df with one sham, one low, one high point
ggplot(a.cv , aes(x=dpi, y=bird_pv))+
  geom_point(aes(color=fct_rev(primary_treatment), alpha="PV"), size=3)+
  geom_line(aes(color=fct_rev(primary_treatment)))+
  geom_point(aes(x=dpi, y=bird_cv, color=fct_rev(primary_treatment), alpha="CV"), size=3)+
  geom_path(aes(dpi, y=bird_cv, color=fct_rev(primary_treatment)), linetype="dashed", alpha=0.5)+
  scale_color_manual(values = c(pri_colors))+
  scale_alpha_manual(values = c(PV=1, CV=0.01))+
  labs(x="Days Post Primary Inoculation", y="Variability", color="Primary Treatment", alpha = "Variability Metric")


#ggsave('ab_var_prim.jpeg', width=6, height=4, units='in')

####Remove Infected Birds for Secondary Analysis####

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

#names
wibird$dpi <- "DPPI -8"
dpi_namespre <- c(
  "DPPI -8" = "DPPI -8"
)

sus_pre_raw<-ggplot(wibird, aes(x=(elisa_od_pre), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi ~log10.sec_dose, labeller = as_labeller(c(dpi_namespre, sec_dose_names)))
sus_pre_raw
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

sus_14_raw<-ggplot(wibird14, aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new, aes(x=(elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Infection Status Upon Secondary Challenge (0|1)", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names14,sec_dose_names)))
sus_14_raw

sus_14_raw7000<-ggplot(wibird14%>%filter(secondary_dose == "7000"), aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new%>%filter(log10.sec_dose == "3.85"), aes(x=(elisa_od_14), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new%>%filter(log10.sec_dose == "3.85"), aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_point(alpha=0.5, size=2)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Infection Status Upon Secondary Challenge (0|1)", color="Secondary Dose", fill ="Secondary Dose", x="Antibody Levels")+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=15),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )
  #scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  
  #facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names14,sec_dose_names)))
sus_14_raw7000


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


sus_41_raw<-ggplot(wibird41, aes(x=(elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="Antibody Levels", y= " ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names41,sec_dose_names)))
sus_41_raw
#ggsave('sus_ab_41.jpeg', width=6, height=4, units='in')

#combine
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

####Ab SID 14
s.ab <- m.ab %>% filter(dpi == 56) %>%
  dplyr::select(bird_ID, dpi, elisa_od, quantity1, primary_treatment, secondary_dose, inf_sec, inf_prim)

ggplot(s.ab, aes(x=as.factor(primary_treatment), y=elisa_od, color=fct_rev(as.factor(inf_sec))))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(alpha=0.75, height=0, aes(shape=as.factor(inf_sec)))+
  # xlim(xmin=0.001, xmax=10000)+
  # scale_x_log10()+
  facet_grid(~secondary_dose)
  
inf.ab <- m.ab %>%
  filter(inf_sec == 1 & dpi %in% c(14, 41, 56))
ggplot(inf.ab, aes(x=dpi, y=elisa_od, color=fct_rev(primary_treatment)))+
  geom_hline(yintercept=0.061, linetype = "dashed", alpha=0.5)+
  geom_point(alpha=0.5)+
  geom_line(aes(group=as.numeric(band_number)))+
  facet_wrap(~secondary_dose)

####Do antibody levels in susceptible birds predict pathology or eye score?####
m.ab.w <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41, 56))

wibird <- m.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)

m.ab.w.tes <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, tes, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41, 56))
wibird.tes <- m.ab.w.tes %>% pivot_wider(names_from = dpi,
                                     names_glue = "{.value}_{dpi}",
                                     values_from = tes)

m.ab.w.quant <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, quantity1, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41, 56))
wibird.quant <- m.ab.w.quant %>% pivot_wider(names_from = dpi,
                                         names_glue = "{.value}_{dpi}",
                                         values_from = quantity1)

wibird.all <- left_join(wibird, wibird.quant, by="band_number")
wibird.all <- left_join(wibird.all, wibird.tes, by="band_number")

ggplot(wibird.all, aes(x=elisa_od_14, y=log10(quantity1_56)))+
  geom_point()+
  geom_smooth(method=glm)
hist(log10(wibird.all$quantity1_56+1))
hist(wibird.all$elisa_od_14)
glm.ab.quant <- glm(log10(quantity1_56+1) ~ elisa_od_41 + primary_treatment, data=wibird.all, family=Gamma())
summary(glm.ab.quant)
simulateResiduals(glm.ab.quant, plot=T)
ggplot(wibird.all, aes(x=log10(quantity1_56+1), y=elisa_od_14, color=fct_rev(primary_treatment), shape=as.factor(inf_sec)))+
  geom_point(size=3)+
  facet_wrap(~primary_treatment)


glm.ab.tes <- glmmTMB(tes_56+0.01 ~ elisa_od_14 + primary_treatment, data=wibird.all, family=Gamma())
summary(glm.ab.tes)
simulateResiduals(glm.ab.tes, plot=T)
plot(allEffects(glm.ab.tes))
ggplot(wibird.all, aes(x=tes_56, y=elisa_od_14, color=fct_rev(primary_treatment), shape=as.factor(inf_sec)))+
  geom_point(size=3)+
  facet_wrap(~primary_treatment)

####Variability With Hawley et al., 2024####
summary_df <- a.cv %>%
  filter(dpi == 41) %>%
  group_by(primary_treatment) %>%
  summarize(
    mean_bird_pv = mean(bird_pv, na.rm = TRUE),
    mean_bird_cv = mean(bird_cv, na.rm = TRUE),
    mean_mean_od = mean(mean_od, na.rm = TRUE)
  )
boot_results
summary_df
dose <- c("Sham", "Low", "High")
CV <- as.numeric(c(0.899, 1.630, 2.511))
ab.pv <- as.numeric(c(0.04156784, 0.14716858, 0.18885516))
ab.cv <- as.numeric(c(0.03997555, 0.17202405, 0.26176249))
ci_lower_cv <- as.numeric(c(0.0230, 0.119, 0.165))
ci_upper_cv <- as.numeric(c(0.0409, 0.206, 0.339))
ci_lower_pv <- as.numeric(c(0.0246, 0.0998, 0.144))
ci_upper_pv <- as.numeric(c(0.0452, 0.178, 0.231))
mean_sus <- as.numeric(c(1.181, 0.446, 0.192))
mean_ab <- as.numeric(c(0.045, 0.04786, 0.0579))
dpi <- 41



het.df <- data.frame(dose, CV, mean_sus, ab.pv, ab.cv, dpi, mean_ab, ci_lower_cv, ci_upper_cv, ci_lower_pv, ci_upper_pv)

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
  geom_errorbar(aes(ymin=ci_lower_pv, ymax=ci_upper_pv, y=ab.pv, x="PV", color=dose), width=0.1)+
  geom_point(size=5,aes(x="CV", y=ab.cv, color=dose, shape = "CV"))+
  geom_errorbar(aes(ymin=ci_lower_cv, ymax=ci_upper_cv, y=ab.cv, x="CV", color=dose), width=0.1)+
  #geom_point(size=3.6,aes(x="PV", y=ab.pv), color="black", shape=1, stroke = 2)+
  scale_color_manual(values=c(pri_colors))+
  scale_y_continuous(
    name="Variability in Antibody Levels DPPI 41"
  )+
  scale_x_discrete(
    name="Metric"
  )+
  scale_shape_manual(values=c(17,19))+
  theme_minimal()+
  theme(
    axis.title.y = element_text(color = "black", size=12, face = "bold"),
    legend.position="top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.background= element_rect(size=0.25, linetype="solid"))+
  labs(color = "Primary Treatment", shape= "Variability Metric")
ab.paper.only

legend <- cowplot::get_legend(ab.paper.only)

# Remove legends from plots
het.paper.only <- het.paper.only + theme(legend.position = "none")
ab.paper.only <- ab.paper.only + theme(legend.position = "none")

#FIGURE 4
Fig4 <- grid.arrange(legend, 
             arrangeGrob(het.paper.only, ab.paper.only, ncol = 2),
             ncol = 1,
             heights = c(1, 4))
#ggsave("var_sus_abs.jpeg", Fig4, width=10, height = 8, units="in", dpi = 300)


####Variability in reinfected birds####
#2449 had 0.5 eye score once but no path load > remove :( # This is now considered uninfected bird and not removed from the analysis
sec.ab <- m.ab %>%
  filter(inf_sec == 1 & dpi > 41 & secondary_dose == 7000 & band_number!=2449)

sec.ab %>%
  filter(secondary_dose == 7000 & dpi ==56)%>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  group_by(secondary_dose, primary_treatment)%>%
  tbl_summary()

m.ab %>%
  filter(dpi == 56)%>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  group_by(secondary_dose, primary_treatment, inf_sec)%>%
  tbl_summary()

#max quantity
max.quant.sec <- sec.ab %>%
  filter(dpi > 41) %>%
  group_by(band_number, primary_treatment, secondary_dose) %>%
  summarise(max_quantity = max(quantity, na.rm = TRUE))

max.tes.sec <- sec.ab %>%
  group_by(band_number, primary_treatment, secondary_dose)%>%
  summarise(max_tes = max(tes, na.rm=TRUE))

max.sec<-left_join(max.quant.sec, max.tes.sec, by="band_number")

max.sec$primary_treatment <- max.sec$primary_treatment.x
max.sec$secondary_dose <- max.sec$secondary_dose.x

max.sec <- max.sec %>%
  dplyr::select(-c(primary_treatment.y, secondary_dose.y, primary_treatment.x, secondary_dose.x))
max.sec$max_quantity1 <- max.sec$max_quantity+1

max.pv <- max.sec %>% 
  group_by(primary_treatment)%>%
  filter(secondary_dose == "7000")%>%
  drop_na(max_quantity1, max_tes)%>%
  dplyr::reframe(mean_quant = mean(max_quantity1),
                 bird_cv_quant = sd(max_quantity1)/mean(max_quantity1),
                 bird_sd_quant = sd(max_quantity1),
                 bird_pv_quant = PV(max_quantity1),
                 max_quantity1 = max_quantity1,
                mean_tes = mean(max_tes),
                bird_cv_tes = sd(max_tes)/mean(max_tes),
                bird_sd_tes = sd(max_tes),
                bird_pv_tes = PV(max_tes+0.001),
                max_tes = max_tes,
                band_number = band_number)

bootstrap_fn_dis <- function(data, indices) {
  # Resample the data (max_quantity and max_tes)
  resampled_data <- data[indices, ]
  
  # Calculate bird_cv and bird_pv for max_quantity
  bird_cv_quant <- ifelse(mean(resampled_data$max_quantity, na.rm = TRUE) != 0, 
                          sd(resampled_data$max_quantity, na.rm = TRUE) / mean(resampled_data$max_quantity, na.rm = TRUE), 
                          NA)  # Return NA if mean is zero
  bird_pv_quant <- ifelse(any(!is.na(resampled_data$max_quantity)), PV(resampled_data$max_quantity), NA)
  
  # Calculate bird_cv and bird_pv for max_tes, with handling for zero values
  bird_cv_tes <- ifelse(mean(resampled_data$max_tes, na.rm = TRUE) > 0, 
                        sd(resampled_data$max_tes, na.rm = TRUE) / mean(resampled_data$max_tes, na.rm = TRUE), 
                        NA)  # Return NA if mean is zero or near-zero
  bird_pv_tes <- ifelse(any(!is.na(resampled_data$max_tes)) && !all(resampled_data$max_tes == 0), 
                        PV(resampled_data$max_tes + 0.001), NA)  # Return NA if all values are 0
  
  # Return all four statistics
  return(c(bird_cv_quant, bird_pv_quant, bird_cv_tes, bird_pv_tes))
}

# Apply bootstrapping to each group (primary_treatment)
boot_results_dis <- max.sec %>%
  group_by(primary_treatment) %>%
  do({
    # Perform bootstrapping (1000 iterations)
    boot_out <- boot(data = ., statistic = bootstrap_fn_dis, R = 1000)
    
    # Calculate 95% confidence intervals for bird_cv and bird_pv for both quantities and tes
    ci_cv_quant <- tryCatch(boot.ci(boot_out, index = 1, type = "perc")$percent[4:5], 
                            error = function(e) c(NA, NA))  # Handle errors if CV is NA
    ci_pv_quant <- tryCatch(boot.ci(boot_out, index = 2, type = "perc")$percent[4:5], 
                            error = function(e) c(NA, NA))
    ci_cv_tes <- tryCatch(boot.ci(boot_out, index = 3, type = "perc")$percent[4:5], 
                          error = function(e) c(NA, NA))
    ci_pv_tes <- tryCatch(boot.ci(boot_out, index = 4, type = "perc")$percent[4:5], 
                          error = function(e) c(NA, NA))
    
    tibble(
      primary_treatment = unique(.$primary_treatment),
      mean_bird_cv_quant = mean(.$bird_cv_quant, na.rm = TRUE),
      mean_bird_pv_quant = mean(.$bird_pv_quant, na.rm = TRUE),
      ci_lower_cv_quant = ci_cv_quant[1],  # Confidence interval for bird_cv_quant
      ci_upper_cv_quant = ci_cv_quant[2],
      ci_lower_pv_quant = ci_pv_quant[1],  # Confidence interval for bird_pv_quant
      ci_upper_pv_quant = ci_pv_quant[2],
      
      mean_bird_cv_tes = mean(.$bird_cv_tes, na.rm = TRUE),
      mean_bird_pv_tes = mean(.$bird_pv_tes, na.rm = TRUE),
      ci_lower_cv_tes = ci_cv_tes[1],  # Confidence interval for bird_cv_tes
      ci_upper_cv_tes = ci_cv_tes[2],
      ci_lower_pv_tes = ci_pv_tes[1],  # Confidence interval for bird_pv_tes
      ci_upper_pv_tes = ci_pv_tes[2]
    )
  }) %>%
  ungroup()

# View the bootstrapped results
boot_results_dis

ggplot(max.pv, aes(x=primary_treatment, y=max_tes, color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=mean_tes), size=10, shape= 95)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_tes*10), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_tes*10), shape=5, size=5, color="black")+
  #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_tes*10, ymin=ci_lower_pv_tes*10, ymax=ci_upper_pv_tes*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  labs(x="Primary Treatment", y="Maximum Eyescore", color="Primary Treatment")+
  scale_y_continuous(limits=c(0, 6.5))+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv, aes(x=max_tes, fill=fct_rev(primary_treatment)))+
  geom_histogram(bins = 15, alpha=0.75)+
  geom_density()+
  geom_vline(aes(group=primary_treatment, xintercept=bird_pv_tes*10), color="black", linetype="dashed")+
  geom_vline(aes(group=primary_treatment, xintercept=mean_tes), color="brown", linetype="dashed")+
  scale_fill_manual(values=c(pri_colors))+
  facet_wrap(~primary_treatment, ncol=1)+
  labs(y="Count", x="Max Eye Score", color="Primary Treatment")

ggplot(max.pv, aes(x=primary_treatment, y=log10(max_quantity1), color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=log10(mean_quant+1)), size=10, shape= 95)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*10), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*10), shape=5, size=5, color="black")+
  #geom_point(aes(x=primary_treatment, y=bird_cv_quant), shape=18, size=6)+
  #geom_point(aes(x=primary_treatment, y=bird_cv_quant), shape=5, size=5, color="black")+
 #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_quant*10, ymin=ci_lower_pv_quant*10, ymax=ci_upper_pv_quant*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  scale_y_continuous(limits =c(0, 7.5))+
  labs(x="Primary Treatment", y="log10(Maximum Pathogen Load+1)", color="Primary Treatment")+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv, aes(x=primary_treatment, y=max_quantity1, color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=mean_quant), shape=95, size=10)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*100000), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*100000), shape=5, size=5, color="black")+
  #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_quant*10, ymin=ci_lower_pv_quant*10, ymax=ci_upper_pv_quant*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  #scale_y_continuous(limits =c(0, 7.5))+
  labs(x="Primary Treatment", y="Maximum Pathogen Load+1)", color="Primary Treatment")+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv, aes(x=log10(max_quantity1), fill=fct_rev(primary_treatment)))+
  geom_histogram(bins = 15, alpha=0.75)+
  geom_density()+
  geom_vline(aes(group=primary_treatment, xintercept=max.pv$bird_pv_quant*10), color="black", linetype="dashed")+
  geom_vline(aes(group=primary_treatment, xintercept=log10(mean_quant+1)), color="brown", linetype="dashed")+
  scale_fill_manual(values=c(pri_colors))+
  facet_wrap(~primary_treatment, ncol=1)+
  labs(y="Count", x="Log10(Max Pathogen Load+1)", fill="Primary Treatment")

####max of all infected birds regardless of secondary treatment####
sec.ab.all <- m.ab %>%
  filter(inf_sec == 1 & dpi > 41 & band_number!=2449)

sec.ab.all %>%
  filter(dpi ==56)%>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  group_by(secondary_dose, primary_treatment)%>%
  tbl_summary()

m.ab %>%
  filter(dpi == 56)%>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec)%>%
  group_by(secondary_dose, primary_treatment, inf_sec)%>%
  tbl_summary()

#max quantity
max.quant.sec.all <- sec.ab.all %>%
  filter(dpi > 41) %>%
  group_by(band_number, primary_treatment, secondary_dose, inf_sec) %>%
  summarise(max_quantity = max(quantity, na.rm = TRUE))

max.tes.sec.all <- sec.ab.all %>%
  group_by(band_number, primary_treatment, secondary_dose, inf_sec)%>%
  summarise(max_tes = max(tes, na.rm=TRUE))

max.sec.all<-left_join(max.quant.sec.all, max.tes.sec.all, by="band_number")

max.sec.all$primary_treatment <- max.sec.all$primary_treatment.x
max.sec.all$secondary_dose <- max.sec.all$secondary_dose.x
max.sec.all$inf_sec <- max.sec.all$inf_sec.x

max.sec.all <- max.sec.all %>%
  dplyr::select(-c(primary_treatment.y, secondary_dose.y, primary_treatment.x, secondary_dose.x, inf_sec.y, inf_sec.x))
max.sec.all$max_quantity1 <- max.sec.all$max_quantity+1

max.pv.all <- max.sec.all %>% 
  group_by(primary_treatment, inf_sec)%>%
  drop_na(max_quantity1, max_tes)%>%
  dplyr::reframe(mean_quant = mean(max_quantity1),
                 bird_cv_quant = sd(max_quantity1)/mean(max_quantity1),
                 bird_sd_quant = sd(max_quantity1),
                 bird_pv_quant = PV(max_quantity1),
                 max_quantity1 = max_quantity1,
                 mean_tes = mean(max_tes),
                 bird_cv_tes = sd(max_tes)/mean(max_tes),
                 bird_sd_tes = sd(max_tes),
                 bird_pv_tes = PV(max_tes+0.001),
                 max_tes = max_tes,
                 band_number = band_number,
                 secondary_dose = secondary_dose)

ggplot(max.pv.all %>% filter(inf_sec==1), aes(x=primary_treatment, y=max_tes, color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=mean_tes), size=10, shape= 95)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_tes*10), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_tes*10), shape=5, size=5, color="black")+
  #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_tes*10, ymin=ci_lower_pv_tes*10, ymax=ci_upper_pv_tes*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  labs(x="Primary Treatment", y="Maximum Eyescore", color="Primary Treatment", title = "All susceptible birds regardless of priming dose")+
  #facet_wrap(~secondary_dose)+
  #scale_y_continuous(limits=c(0, 6.5))+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv.all, aes(x=max_tes, fill=fct_rev(primary_treatment)))+
  geom_histogram(bins = 12, alpha=0.75)+
  geom_density()+
  geom_vline(aes(group=primary_treatment, xintercept=bird_pv_tes*10), color="black", linetype="dashed")+
  geom_vline(aes(group=primary_treatment, xintercept=bird_cv_tes), color="pink", linetype="dashed")+
  geom_vline(aes(group=primary_treatment, xintercept=mean_tes), color="brown", linetype="dashed")+
  scale_fill_manual(values=c(pri_colors))+
  facet_wrap(~primary_treatment, ncol=1)+
  labs(y="Count", x="Max Eye Score", fill="Primary Treatment", title = "All susceptible birds regardless of priming dose")

ggplot(max.pv.all, aes(x=primary_treatment, y=log10(max_quantity1), color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=log10(mean_quant+1)), size=10, shape= 95)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*10), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*10), shape=5, size=5, color="black")+
  #geom_point(aes(x=primary_treatment, y=bird_cv_quant), shape=18, size=6)+
  #geom_point(aes(x=primary_treatment, y=bird_cv_quant), shape=5, size=5, color="black")+
  #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_quant*10, ymin=ci_lower_pv_quant*10, ymax=ci_upper_pv_quant*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  #scale_y_continuous(limits =c(0, 7.5))+
  labs(x="Primary Treatment", y="log10(Maximum Pathogen Load+1)", color="Primary Treatment", title = "All susceptible birds regardless of priming dose")+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv.all, aes(x=primary_treatment, y=max_quantity1, color=fct_rev(primary_treatment)))+
  #geom_boxplot(width=0.25, outlier.alpha = 0)+
  geom_point(aes(x=primary_treatment, y=mean_quant), shape=95, size=10)+
  geom_jitter(alpha=0.5, height=0, width=0.125)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*100000), shape=18, size=6)+
  geom_point(aes(x=primary_treatment, y=bird_pv_quant*100000), shape=5, size=5, color="black")+
  #geom_errorbar(data=boot_results_dis, aes(x=primary_treatment, y=mean_bird_pv_quant*10, ymin=ci_lower_pv_quant*10, ymax=ci_upper_pv_quant*10), width=0.1, color="black")+
  scale_color_manual(values=c(pri_colors))+
  #scale_y_continuous(limits =c(0, 7.5))+
  labs(x="Primary Treatment", y="Maximum Pathogen Load+1)", color="Primary Treatment", title = "All susceptible birds regardless of priming dose")+
  theme(
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16)
  )

ggplot(max.pv.all, aes(x=log10(max_quantity1), fill=fct_rev(primary_treatment)))+
  geom_histogram(bins = 50, alpha=0.75)+
  geom_density()+
  geom_vline(aes(group=primary_treatment, xintercept=bird_pv_quant*10), color="black", linetype="dashed")+
  geom_vline(aes(group=primary_treatment, xintercept=log10(mean_quant+1)), color="brown", linetype="dashed")+
  scale_fill_manual(values=c(pri_colors))+
  facet_wrap(~primary_treatment, ncol=1)+
  labs(y="Count", x="Log10(Max Pathogen Load+1)", fill="Primary Treatment", title = "All susceptible birds regardless of priming dose")
