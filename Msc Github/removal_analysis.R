#### Supplementary to Host_Response.Rmd ####
#Analysis to compare statitics and variability in antibody levels across priming phase of EEID 1A
#Fifteen plasma samples were lost from DPPI 14 sampling and one from DPPI 41 sampling, resulting in reduced sample sizes for antibody analyses.
  #Of the 15 samples lost from DPPI 14, 14 were sham inoculated and one was inoculated with a low dose. 
  #The one plasma sample lost from DPPI 41 was from a bird inoculated with a low dose.

#This analysis will compare the results from the first part of the Host_Response analysis looking specifically at the effect of removing birds whose plasma
  #samples were lost from each day compared to only removing the missing data points for the day they were missing.

####read-in and format data####
rm(list=ls())
library(tidyverse)
library(dplyr)
library(glmmTMB)
library(effects)
library(AICcmodavg)
library(emmeans)

source("dataCleaning_EEID1A.R")

#set colors for primary treatment
pri_colors <- c("#1B9E77", "#7570B3", "#D95F02")
sec_colors <- c("#8C754B", "#77AB59", "#59A5D8", "#9F77D9", "#FA8072")
#set theme
theme_set(theme_bw())

#threshold cutoffs
m.ab$threshold_cutoff = 50
m.ab$seropos_cutoff = 0.061
m.ab$pathology_cutoff = 0

#make infected criteria
m.ab <- m.ab %>%
  # Sort data by band_number and dpi
  arrange(band_number, dpi) %>%
  # Replace NA values for quantity only
  mutate(quantity_clean = coalesce(quantity, 0)) %>%
  # Flag as infected if quantity > 50 on any single day
  mutate(quantity_flag = ifelse(quantity_clean > 50, 1, 0)) %>%
  # Check consecutive days for tes > 0.5, ignoring NAs for tes
  group_by(band_number) %>%
  mutate(
    tes_flag = ifelse(tes > 0.5 & !is.na(tes), 1, 0),
    tes_consecutive = tes_flag + lag(tes_flag, default = 0),
    tes_single_high = ifelse(tes >= 1 & !is.na(tes), 1, 0)
  ) %>%
  # Combine the three criteria: quantity > 50 OR tes > 0.5 on consecutive days OR tes >= 1 on any day
  mutate(inf = ifelse(quantity_flag == 1 | tes_consecutive > 1 | tes_single_high == 1, 1, 0)) %>%
  # Calculate inf_pri and inf_sec based on dpi
  ungroup() %>%
  group_by(band_number) %>%
  mutate(
    inf_pri = ifelse(any(inf == 1 & dpi < 42), 1, 0),
    inf_sec = ifelse(any(inf == 1 & dpi > 42), 1, 0)
  ) %>%
  # Un-group the data if needed
  ungroup() %>%
  # Select relevant columns (optional)
  select(-quantity_clean, -quantity_flag, -tes_flag, -tes_consecutive, -tes_single_high)

#omit samples (not ELISA)
#omit 2505 from data set for analysis as it was positive at quarantine (n=1)
m.ab <- m.ab %>%
  filter(band_number != 2505)

#omit 2452 and 2562 bc positive qpcr and they are shams (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2452, 2562)) 

#omit sham birds with eye scores (0.5) (2419, 2434) (n=2)
m.ab <- m.ab %>% filter(!band_number %in% c(2419, 2434))

#omit sham bird with path load (2514) (n=1)
m.ab <- m.ab %>% filter(band_number != 2514)

m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  ) %>%
  modify_header(
    label ~ "**Sample Sizes for Primary Analysis**"
  )
#Final Primary Sample Sizes: Sham = 46; Low = 51; High = 53

####Begin Analysis####
#data frame with only primary infection
p.ab <- m.ab %>%
  filter(dpi <=41)

p.ab %>%
  dplyr::select(dpi, primary_treatment, elisa_od)%>%
  tbl_summary(
    by=dpi
  ) %>%
  modify_header(
    label ~ "**PRIMARY ANALYSIS FINAL SAMPLE SIZES**"
  )

#Which samples did not have plasma samples?
missing <- p.ab %>% 
  filter(dpi %in% c(-8, 14, 41) & is.na(elisa_od)) %>%
  dplyr::select(band_number, dpi, primary_treatment, secondary_dose, elisa_od)
missing
#remove birds without plasma samples for antibody analysis
# Identify birds missing plasma samples on any of the specified days
birds_missing_samples <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41) & is.na(elisa_od)) %>%
  pull(band_number)  # Extract bird IDs with missing samples

# Remove these birds from the dataset entirely
rem <- p.ab %>%
  filter(!(band_number %in% birds_missing_samples))

rem %>%
  dplyr::select(dpi, primary_treatment, elisa_od)%>%
  tbl_summary(
    by=dpi
  ) %>%
  modify_header(
    label ~ "**Full Removal**"
  )

# Remove elisa_od NAs 
p.ab <- p.ab %>% 
  filter(!(dpi %in% c(-8, 14, 41) & is.na(elisa_od)))

#sample sizes
p.ab %>%
  dplyr::select(dpi, primary_treatment, elisa_od)%>%
  tbl_summary(
    by=dpi
  ) %>%
  modify_header(
    label ~ "**No Full Removal**"
  )

#Baseline; Does not affect baseline differences
lm0 <- lm(elisa_od~primary_treatment, data=p.ab %>% filter(dpi <0))
lm0.r <- lm(elisa_od~primary_treatment, data=rem %>% filter(dpi <0))
summary(lm0)
summary(lm0.r)
plot(allEffects(lm0))
simulateResiduals(lm0, plot=T)

emm_r <- emmeans(lm0, ~primary_treatment)
pairs(emm_r)

summary(lm0.r)
simulateResiduals(lm0.r, plot=T)

#compare r-squared
summary(lm0)$adj.r.squared
summary(lm0.r)$adj.r.squared

#compare Residual Standard Error; lower RSE = model fits data better
summary(lm0)$sigma
summary(lm0.r)$sigma

par(mfrow=c(1,2))
plot(lm0$fitted.values, lm0$residuals, main="Model lm0 Residuals")
plot(lm0.r$fitted.values, lm0.r$residuals, main="Model lm0.r Residuals")

#full days
#I treat days post inoculation as a factor in these models and use band_number as a random effect
p.abt <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41))
p.abt$dpi.f <- as.factor(p.abt$dpi)

p.abt.r <- rem %>%
  filter(dpi %in% c(-8, 14, 41))
p.abt.r$dpi.f <- as.factor(p.abt.r$dpi)

#model comparison no removal
p1 <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2 <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2.5 <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), data=p.abt, family=Gamma(log))
p3<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p4<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p5 <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))

aictab(cand.set=list(p1, p2, p2.5, p3, p4, p5), modnames=c("p1", "p2", "p2.5", "p3", "p4", "p5"))

#model comparison withg removal
p1.r <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log))
p2.r <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log))
p2.5.r <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), data=p.abt.r, family=Gamma(log))
p3.r <- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log))
p4.r <- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log))
p5.r <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log))

aictab(cand.set=list(p1.r, p2.r, p2.5.r, p3.r, p4.r, p5.r), modnames=c("p1.r", "p2.r", "p2.5.r", "p3.r", "p4.r", "p5.r"))

#the model without removal has lower AICc values, indicating better model performance. Removal makes models perform worse.
aictab(cand.set=list(p2, p2.r, p2.5, p2.5.r), modnames = c("p2", "p2.r", "p2.5", "p2.5.r"))

library(VGAM)
drop1(p2, test="Chisq")
drop1(p2.r, test="Chisq")

#final model; does not affect significance or effect sizes
#make Sham the reference category
p.abt$primary_treatment <- relevel(p.abt$primary_treatment, ref = "Sham")
p.abt.r$primary_treatment <- relevel(p.abt.r$primary_treatment, ref = "Sham")

lm1 <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log)) 
lm1.r <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number), data=p.abt.r, family=Gamma(log)) 

summary(lm1)
summary(lm1.r)
#intercept > expected log-mean elisa_od (-3.074 ~ 0.046) looks good

simulateResiduals(lm1, plot=T)
simulateResiduals(lm1.r, plot=T)
hist(resid(lm1))
hist(resid(lm1.r))

emm_results <- emmeans(lm1, ~ primary_treatment|dpi.f, scale="response")
pairwise <- pairs(emm_results, adjust="tukey")
(summary(pairwise))

emm_results.r <- emmeans(lm1.r, ~ primary_treatment|dpi.f, scale="response")
pairwise.r <- pairs(emm_results.r, adjust="tukey")
(summary(pairwise.r))

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

#specify model to predict
mod.pred <- lm1.r
#Model predictions
dat.new.r=expand.grid(primary_treatment=unique(p.abt$primary_treatment),
                    elisa_od = unique(p.abt$elisa_od),
                    dpi.f = unique(p.abt$dpi.f))
#generate yhat values
dat.new.r$yhat = predict(mod.pred, type="response", newdata=dat.new.r, re.form=NA)
#generate SEM
preds = predict(mod.pred, type="link", newdata=dat.new.r, se.fit=TRUE, re.form=NA)
dat.new.r = cbind(dat.new.r, preds)

ilink <- family(mod.pred)$linkinv
dat.new.r <- transform(dat.new.r,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

#not removed
ggplot()+
  #geom_jitter(data=p.abt, aes(x=dpi.f, y=elisa_od, color=primary_treatment), alpha=0.5)+
  geom_point(data=dat.new, aes(x=dpi.f, y=yhat, color=primary_treatment), shape=18, size=3)+
  geom_errorbar(data=dat.new, aes(x=dpi.f, y=yhat, ymin=Lower, ymax=Upper, color=primary_treatment), width=0.1)

#removed
ggplot()+
  #geom_jitter(data=p.abt.r, aes(x=dpi.f, y=elisa_od, color=primary_treatment), alpha=0.5)+
  geom_point(data=dat.new.r, aes(x=dpi.f, y=yhat, color=primary_treatment), shape=18, size=3)+
  geom_errorbar(data=dat.new.r, aes(x=dpi.f, y=yhat, ymin=Lower, ymax=Upper, color=primary_treatment), width=0.1)

####Variability####
    #add small constant to total eye score for calculation of variability
    p.ab$tes.new <- p.ab$tes + .01
    #new df with only individuals that were infected any time during primary challenge
    p.ab$dpi.f <- as.factor(p.ab$dpi)
    p.abi <- p.ab %>%
      filter(inf_pri == 1)%>%
      dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)
    
    #df with all individuals during primary challenge
    p.aba <-  p.ab %>%
      dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)
    
    #add small constant to total eye score for calculation of variability
    rem$tes.new <- rem$tes + .01
    #new df with only individuals that were infected any time during primary challenge
    rem$dpi.f <- as.factor(rem$dpi)
    p.abi.r <- rem %>%
      filter(inf_pri == 1)%>%
      dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)
    
    #df with all individuals during primary challenge
    p.aba.r <-  rem %>%
      dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)

# Define functions for calculating CV, PV, and V2
calculate_cv <- function(data) {
  return(sd(data) / mean(data))
}

calculate_pv <- function(data) {
  return(PV(data))
}

calculate_v2 <- function(data) {
  mean_x <- mean(data)
  sd_x <- sd(data)
  v2 <- sd_x^2 / (sd_x^2 + mean_x^2)
  return(v2)
}

# Set the number of bootstrap replicates
n_boot <- 1000

# Pipeline to calculate CV, PV, V2, and their 95% confidence intervals
a.cv <- p.aba %>%
  group_by(dpi.f, primary_treatment) %>%
  drop_na(elisa_od) %>%
  dplyr::reframe(
    # Calculate means and standard deviations
    elisa_od = elisa_od,
    mean_od = mean(elisa_od),
    bird_sd = sd(elisa_od),
    band_number = band_number,
    bird_cv = calculate_cv(elisa_od),   # Calculate CV for the original data
    bird_pv = calculate_pv(elisa_od),   # Calculate PV for the original data
    bird_v2 = calculate_v2(elisa_od),   # Calculate V2 for the original data
    
    # Bootstrap for each metric
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(elisa_od, replace = TRUE)))),#1000 new CV calculations from sampling the data 10 times randomly
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(elisa_od, replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(elisa_od, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for CV, PV, and V2
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)), #lower confidence interval (2.5th percentile)
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)), #upper confidence interval (97.5th percentile)
    
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975)),
    
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Optionally, remove the bootstrap columns to clean up the dataset
  select(-cv_bootstrap, -pv_bootstrap, -v2_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.m <- left_join(p.aba, a.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.m <- p.aba.m %>%
  select(-elisa_od.y)
p.aba.m$elisa_od <- p.aba.m$elisa_od.x
# p.abi.m$primary_treatment <- p.abi.m$primary_treatment.x
# p.abi.m <- p.abi.m %>%
#   select(-elisa_od.x, -primary_treatment.x)

#just days samples were taken
p.aba.m <- p.aba.m %>%
  filter(dpi.f %in% c(-8, 14, 41))

summary_tibble <- a.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    mean_bird_cv = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    mean_bird_sd = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    mean_bird_pv = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    avg_elisa_od = mean(elisa_od, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble)

# Removed variability
a.cv.r <- p.aba.r %>%
  group_by(dpi.f, primary_treatment) %>%
  drop_na(elisa_od) %>%
  dplyr::reframe(
    # Calculate means and standard deviations
    elisa_od = elisa_od,
    mean_od = mean(elisa_od),
    bird_sd = sd(elisa_od),
    band_number = band_number,
    bird_cv = calculate_cv(elisa_od),   # Calculate CV for the original data
    bird_pv = calculate_pv(elisa_od),   # Calculate PV for the original data
    bird_v2 = calculate_v2(elisa_od),   # Calculate V2 for the original data
    
    # Bootstrap for each metric
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(elisa_od, replace = TRUE)))),#1000 new CV calculations from sampling the data 10 times randomly
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(elisa_od, replace = TRUE)))),
    v2_bootstrap = list(replicate(n_boot, calculate_v2(sample(elisa_od, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for CV, PV, and V2
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)), #lower confidence interval (2.5th percentile)
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)), #upper confidence interval (97.5th percentile)
    
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975)),
    
    v2_lower_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.025)),
    v2_upper_ci = map_dbl(v2_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Optionally, remove the bootstrap columns to clean up the dataset
  select(-cv_bootstrap, -pv_bootstrap, -v2_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.m.r <- left_join(p.aba.r, a.cv.r, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.m.r <- p.aba.m.r %>%
  select(-elisa_od.y)
p.aba.m.r$elisa_od <- p.aba.m.r$elisa_od.x
# p.abi.m$primary_treatment <- p.abi.m$primary_treatment.x
# p.abi.m <- p.abi.m %>%
#   select(-elisa_od.x, -primary_treatment.x)

#just days samples were taken
p.aba.m.r <- p.aba.m.r %>%
  filter(dpi.f %in% c(-8, 14, 41))

summary_tibble.r <- a.cv.r %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    mean_bird_cv = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    mean_bird_sd = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    mean_bird_pv = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    avg_elisa_od = mean(elisa_od, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
print(summary_tibble.r)

#write.csv(summary_tibble.r, 
#"/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/ab_var_removed.csv", row.names=FALSE)

#write.csv(summary_tibble, 
#"/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/ab_var_not_removed.csv", row.names=FALSE)

#differences in pv
ggplot()+
  geom_point(data=summary_tibble, aes(x=dpi.f, y=mean_bird_pv, color=primary_treatment))+
  geom_errorbar(data=summary_tibble, aes(x=dpi.f, y=mean_bird_pv, ymin=lower_ci_pv, ymax=upper_ci_pv, color=primary_treatment), width=0.2)+
  geom_point(data=summary_tibble.r, aes(x=dpi.f, y=mean_bird_pv, color=primary_treatment), shape=17)+
  geom_errorbar(data=summary_tibble.r, aes(x=dpi.f, y=mean_bird_pv, ymin=lower_ci_pv, ymax=upper_ci_pv, color=primary_treatment), width=0.2, lty="dashed")+
  scale_color_manual(values=c(pri_colors))

#differences in cv
ggplot()+
  geom_point(data=summary_tibble, aes(x=dpi.f, y=mean_bird_cv, color=primary_treatment))+
  geom_errorbar(data=summary_tibble, aes(x=dpi.f, y=mean_bird_cv, ymin=lower_ci_cv, ymax=upper_ci_cv, color=primary_treatment), width=0.2)+
  geom_point(data=summary_tibble.r, aes(x=dpi.f, y=mean_bird_cv, color=primary_treatment), shape=17)+
  geom_errorbar(data=summary_tibble.r, aes(x=dpi.f, y=mean_bird_cv, ymin=lower_ci_cv, ymax=upper_ci_cv, color=primary_treatment), width=0.2, lty="dashed")+
  scale_color_manual(values=c(pri_colors))


####Secondary Susceptibility####
#omit from dataset - were still infected (path load) on day 41 prior to reinfection day 42
m.ab <- m.ab %>%
  filter(!(band_number %in% c(2274, 2469, 2520, 2494)))

rem <- m.ab %>%
  filter(!(band_number %in% birds_missing_samples))
rem <- rem %>%
  filter(!(band_number %in% c(2274, 2469, 2520, 2494)))

m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES SECONDARY NO REMOVAL**"
  )

rem %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES SECONDARY WITH REMOVAL**"
  )

m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary()%>%
  modify_header(
    label ~ "**By Primary and Secondary Dose**"
  )
m.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  )

rem$log10.sec_dose <- log10(rem$secondary_dose+1)
rem$dpi.f <- as.factor(rem$dpi)
m.ab$log10.sec_dose <- log10(m.ab$secondary_dose+1)
m.ab$dpi.f <- as.factor(m.ab$dpi)
#Model with removal
remmod14 <- rem %>%
  filter(dpi.f == 14)%>%
  dplyr::select(inf_sec, elisa_od, dpi.f, dpi, log10.sec_dose, band_number)

remmod41 <- rem %>%
  filter(dpi.f == 41)%>%
  dplyr::select(inf_sec, elisa_od, dpi.f, dpi, log10.sec_dose, band_number)

m.ab14 <- m.ab %>%
  filter(dpi.f == 14)%>%
  dplyr::select(inf_sec, elisa_od, dpi.f, dpi, log10.sec_dose, band_number)

m.ab41 <- m.ab %>%
  filter(dpi.f == 41)%>%
  dplyr::select(inf_sec, elisa_od, dpi.f, dpi, log10.sec_dose, band_number)


#I want to see whether elisa_od on dpi -8, 14, or 41 predict inf_sec.
glm.ab14 <- glm(inf_sec ~ elisa_od + log10.sec_dose, data=remmod14, family=binomial())
glm.m.ab14 <- glm(inf_sec ~ elisa_od + log10.sec_dose, data=m.ab14, family=binomial())

glm.ab41 <- glm(inf_sec ~ elisa_od + log10.sec_dose, data=remmod41, family=binomial())
glm.m.ab41 <- glm(inf_sec ~ elisa_od + log10.sec_dose, data=m.ab41, family=binomial())

summary(glm.ab14)
summary(glm.m.ab14)
summary(glm.ab41)
summary(glm.m.ab41)
simulateResiduals(glm.ab14, plot=T)
simulateResiduals(glm.m.ab14, plot=T)
simulateResiduals(glm.ab41, plot=T)
simulateResiduals(glm.m.ab41, plot=T)


mod <- glm.ab14
dat.new=expand.grid(log10.sec_dose=unique(remmod$log10.sec_dose),
                    inf_sec=unique(remmod$inf_sec),
                    # elisa_od = seq(min(remmod$elisa_od, na.rm=TRUE), 
                    #                max(remmod$elisa_od, na.rm=TRUE), length.out = 10),
                    elisa_od = unique(remmod$elisa_od),
                    dpi.f = unique(remmod$dpi.f))
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

ggplot(remmod, aes(x=(elisa_od), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_point(alpha=0.5)+
  geom_line(data=dat.new, aes(x=(elisa_od), y=yhat, color=as.factor(log10.sec_dose)), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_wrap(~dpi.f~log10.sec_dose, ncol=5)

