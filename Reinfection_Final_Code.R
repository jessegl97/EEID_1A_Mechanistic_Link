####Priming Dose Augments Host Response to Reinfection####
##Final Analysis
## 27 Jan 2025
## JGL

####Load Packages + Read in Data####
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
library(purrr)
#install.packages("remotes")
#remotes::install_github("T-Engel/CValternatives")
library(CValternatives)
library(patchwork)
library(gridExtra)
library(car)
library(VGAM)
library(writexl)
library(dunn.test)

#import data
source("dataCleaning_EEID1A.R")

#DFs: 
  #ab = raw antibody data
  #master = master df
  #m.ab = merged + formatted working df
  #not_recovered = birds that were not recovered by secondary infection: 2274, 2469, 2520, 2494 

#Infected column: inf = if quantity (where any NA is replaced with 0) is greater than cutoff of 50 copies or 
  #total eye score is greater or equal to 1 for at least one day, return 1, else return 0

#Format theme
#set colors for primary treatment
pri_colors <- c("#1B9E77", "#7570B3", "#D95F02")
#set color scheme secondary treatment
sec_colors <- c("#8C754B", "#77AB59", "#59A5D8", "#9F77D9", "#FA8072")
#set theme
theme_set(theme_bw())


m.ab %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=primary_treatment
  ) %>%
  modify_header(
    label ~ "**Sample Sizes for Primary Analysis Primary Treatment**"
  )
#Final Primary Sample Sizes: Sham = 46; Low = 51; High = 53

####Analysis 1) Do antibody levels vary across dpi and primary treatment?####
  #>Fifteen plasma samples were lost from DPPI 14 sampling and one from DPPI 41 sampling, resulting in reduced sample sizes
  #> for antibody analyses. Of the 15 samples lost from DPPI 14, 14 were sham inoculated and one was inoculated with a low dose. 
  #> The one plasma sample lost from DPPI 41 was from a bird inoculated with a low dose.
    #>I went back and compared analyses with removal of these birds from the entire dataset (see *removal_analysis.R*) and found
    #> that *removing these 15 birds did not change variability, effect sizes, or significance.* 
    #> Additionally, *the models without full removal had lower AICc values*, indicating better model performance. 
    #> Removal makes models perform worse.

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
#n=16
#2385, 2395, 2407, 2417, 2439, 2450, 2460, 2489, 2498, 2511, 2537, 2546, 2563, 2565, 2496, 2375
p.ab.missing <- p.ab %>% 
  filter(dpi %in% c(-8, 14, 41) & is.na(elisa_od)) %>%
  dplyr::select(band_number, dpi, primary_treatment, secondary_dose, elisa_od)

p.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**Before Removal**"
  )

p.ab.missing %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**To Remove**"
  )

# Remove elisa_od NAs; the samples that are missing
p.ab <- p.ab %>% 
  filter(!(dpi %in% c(-8, 14, 41) & is.na(elisa_od)))

p.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**After Removal**"
  )

##Baseline antibody levels did not differ between treatments (pre-inoculation)
lm0 <- glmmTMB(elisa_od~primary_treatment, data=p.ab %>% filter(dpi == '-8'), family=Gamma())
summary(lm0)
simulateResiduals(lm0, plot=T)

emm_r <- emmeans(lm0, ~primary_treatment)
pairs(emm_r)

#Do antibody levels differ across sampling days and treatment?
#I treat days post inoculation as a factor in these models and use band_number as a random effect

#Antibodies only measured dpi -8, 14, and 41
p.abt <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41))
p.abt$dpi.f <- as.factor(p.abt$dpi)


#Sample sizes
p.abt %>%
  dplyr::select(dpi.f, primary_treatment)%>%
  tbl_summary(
    by=dpi.f
  )%>%
  modify_header(
    label ~ "**Variability Sample Size**"
  )

#model comparison 
p1 <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2 <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p3 <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), data=p.abt, family=Gamma(log))
p4<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p5<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p6 <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p1d <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p2d <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p3d <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p4d<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p5d<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p6d <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, dispformula = ~dpi.f, family=Gamma(log))

aictab(cand.set=list(p1, p2, p3, p4, p5, p6, p1d, p2d, p3d, p4d, p5d, p6d), modnames=c("p1", "p2", "p3",  "p4", "p5", "p6", "p1d", "p2d", "p3d",  "p4d", "p5d", "p6d"))

p2 <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2d <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), dispformula = ~dpi.f*primary_treatment, data=p.abt, family=Gamma(log))
p2dd <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
aictab(cand.set=list(p2, p2d, p2dd), modnames=c("p2", "p2d", "p2dd"))

# Interaction between primary_treatment and dpi.f without sex with dispersion sub-model is best model by AICc

# Model selection based on AICc:
#   
#   K     AICc Delta_AICc AICcWt Cum.Wt      LL
# p2d 13 -3030.65       0.00   0.71   0.71 1528.76
# p3d 14 -3028.91       1.74   0.29   1.00 1528.96
# p6d  7 -2881.45     149.21   0.00   1.00 1447.85
# p1d  9 -2880.85     149.80   0.00   1.00 1449.64
# p4d 10 -2879.11     151.55   0.00   1.00 1449.81
# p5d 12 -2875.03     155.62   0.00   1.00 1449.89
# p2  11 -2796.41     234.25   0.00   1.00 1409.52
# p3  12 -2795.81     234.84   0.00   1.00 1410.28
# p1   7 -2608.78     421.87   0.00   1.00 1311.52
# p4   8 -2607.81     422.85   0.00   1.00 1312.07
# p5  10 -2604.85     425.81   0.00   1.00 1312.68
# p6   5 -2531.83     498.83   0.00   1.00 1270.98

#model comparison VGAM
# drop1(p1, test="Chisq") #primary_treatment and dpi.f contribute uniquely to explaining variation
# 
# drop1(p2, test="Chisq") #The interaction between primary_treatment and dpi.f > the effect of primary treatment
#                         #on the outcome varies significantly by dpi levels
# 
# drop1(p3, test="Chisq") # sex is not significant
# 
# drop1(p4, test="Chisq") # sex is not significant but interaction b/t primary_treatment and dpi.f is better
# 
# drop1(p5, test="Chisq") # primary_treatment*sex is not significant
# 
# drop1(p6, test="Chisq") # better with primary_treatment

drop1(p1d, test="Chisq") #primary_treatment and dpi.f contribute uniquely to explaining variation

drop1(p2d, test="Chisq") #The interaction between primary_treatment and dpi.f > the effect of primary treatment
#on the outcome varies significantly by dpi levels

drop1(p3d, test="Chisq") # sex is not significant

drop1(p4d, test="Chisq") # sex is not significant but interaction b/t primary_treatment and dpi.f is better

drop1(p5d, test="Chisq") # primary_treatment*sex is not significant

drop1(p6d, test="Chisq") # better with primary_treatment

#make Sham the reference category
p.abt$primary_treatment <- relevel(p.abt$primary_treatment, ref = "Sham")

hist(log(p.abt$elisa_od))

#Final Model; Does the interaction between priming treatment and dpi predict antibody response while controlling for individual bird id?
  #Dispersion formula to account for gamma model allowing dispersion parameter to vary by dpi and primary treatment
    #In gamma models, the variance is proportional to the squared mean
lm1 <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number),
                          dispformula = ~dpi.f*primary_treatment, 
                          family=Gamma(link="log"),
                          data=p.abt)

plot(residuals(lm1, type="pearson"))
simulateResiduals(lm1, plot=T)
plot(allEffects(lm1))
hist(p.abt$elisa_od)
summary(lm1)


testDispersion(lm1)
car::Anova(lm1, type = "III")
#install.packages("broom.mixed")
library(broom.mixed)  # For extracting `glmmTMB` results

# Extract model summary as data frames
disp_summary <- summary(lm1)$coefficients$disp
dispersion_effects <- data.frame(
  Term = rownames(disp_summary),
  Estimate = disp_summary[, 1],
  Std_Error = disp_summary[, 2],
  z_value = disp_summary[, 3],
  p_value = disp_summary[, 4],
  row.names = NULL
)
random_effects <- tidy(lm1, effects = "ran_pars")  # Extract random effects



#Post-Hoc Tukey
emm_results <- emmeans(lm1, ~ primary_treatment|dpi.f, scale="response")
pairwise <- pairs(emm_results, adjust="tukey")
(summary(pairwise))

# Convert pairwise comparisons to a data frame
pairwise_df <- as.data.frame(summary(pairwise))

# Save to an Excel file with multiple sheets
# write_xlsx(
#   list(
#     Fixed_Effects = fixed_effects,
#     Dispersion_Model = dispersion_effects,
#     Random_Effects = random_effects,
#     Pairwise_Comparisons = pairwise_df
#   ),
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/lm1_model_output.xlsx"
# )



#Model Predictions
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

#FIGURE 2 raw data + predictions
fig2.raw<-ggplot(data = p.abt, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +
  # Original jittered data points for elisa_od
  geom_jitter(size = 1.5, alpha = 0.5, width = 0.15, height = 0, shape=16) +
  
  # Line and points with error bars for predicted values
  geom_path(data = dat.new, aes(x = dpi.f, y = yhat, group = primary_treatment,
                                color = primary_treatment), alpha=0.75) +
  geom_errorbar(data = dat.new, aes(ymin = Lower, ymax = Upper, x = dpi.f, y = yhat),
                color = "black", width = 0.01) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat, color = primary_treatment),
             size = 2, alpha = 0.75, shape = 16) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat),
             color = "black", size = 2, alpha = 1, shape = 1, stroke = 0.1) +
  
  # Labels and axis
  labs(y = "Anti-MG IgY Antibody Levels [OD]", 
       x = "Days Post Primary Inoculation", 
       shape = "Variability Metric", 
       color = "Primary Treatment", 
       fill = "Primary Treatment") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  # Theme settings
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

fig2.raw


#Bayesian Gamma Model
# install.packages("brms")
# library(brms)

#Model antibody levels as Gamma-distributed with mean predicted by primary treatment * dpi.f with dispersion modeled explicitly by dpi.f
gamma_model <- readRDS("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Models/gamma_model.rds")
# gamma_model <- brm(
#   bf(elisa_od ~ primary_treatment * dpi.f + (1|band_number),
#      shape ~ dpi.f),
#   family = Gamma(link="log"),
#   data = p.abt,
#   chains=4, cores=4,
#   iter=6000, warmup=2000,
#   control=list(adapt_delta=0.999, max_treedepth=20),
#   prior = c(
#     prior(normal(-3, 1), class="Intercept"), # baseline around log(0.05) ≈ -3, SD=1
#     prior(normal(0, 1), class="b"),          # differences among groups around 0, SD=1
#     prior(exponential(1), class="sd"),       # positive, not too large random effect
#     prior(normal(0, 1), dpar="shape")        # mild prior for shape (dispersion)
#   ),
#   save_pars = save_pars(all = TRUE)
# )

#save model
#saveRDS(gamma_model, file = "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Models/gamma_model.rds")
summary_output <- summary(gamma_model)
#capture.output(summary_output, file = "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/bayes_ab_gamma_model_summary.txt")

# Check convergence diagnostics:
plot(gamma_model)

# Check posterior predictive checks:
pp_check(gamma_model)


# Posterior predictive checks
pp_check(gamma_model)

summary(gamma_model)

#loo_gamma<-loo(gamma_model, moment_match = TRUE, reloo = TRUE)
#saveRDS(loo_gamma, file = "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Models/loo_gamma.rds")
#capture.output(print(loo_gamma), file = "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/loo_gamma_summary.txt")
loo_gamma <- readRDS("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Models/loo_gamma.rds")

loo_gamma
AIC(lm1)

#LOOIC = -3035.2, AIC = -3183.5

#install.packages("tidybayes")
library(tidybayes)

# Generate posterior predictions
preds <- conditional_effects(gamma_model, effects = "dpi.f:primary_treatment", method = "posterior_epred")

# Convert to dataframe
pred_df <- as.data.frame(preds$`dpi.f:primary_treatment`)

# Plot predicted antibody levels over dpi and treatments
ggplot(pred_df, aes(x = dpi.f, y = estimate__, color = primary_treatment, group = primary_treatment)) +
  geom_jitter(data = p.abt, aes(x = dpi.f, y = elisa_od, color = primary_treatment),
              size = 1.5, alpha = 0.5, width = 0.15, height = 0, shape=16) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.05, lty="dashed") +
  theme_minimal() +
  scale_color_manual(values=c(pri_colors))+
  labs(
    title = "Predicted ELISA OD by dpi and Treatment",
    x = "Days Post Infection (dpi)",
    y = "Predicted ELISA OD",
    color = "Treatment"
  )


####Variability in Antibody Levels across primary####
#add small constant to total eye score for calculation of variability
p.ab$tes.new <- p.ab$tes + .01
p.ab$dpi.f <- as.factor(p.ab$dpi)

#new df with only individuals that were infected any time during primary challenge; p.abi = p.ab infected
p.abi <- p.ab %>%
  filter(inf_pri == 1)%>%
  dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)

#df with all individuals during primary challenge regardless of reinfection; p.aba = p.ab all
p.aba <-  p.ab %>%
  dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)

#> I calculate variability the same way across analyses in this experiment. 
  #> I use the below functions calculate_cv, calculate_pv, and calculate_v2. 
  #> I then bootstrap the original data by resampling the data with replacement (each bootstrap sample is created by
  #> randomly re-sampling all points in a treatment group - some points may appear multiple times in the same bootstrap
  #> sample, while others may be excluded). The 95% confidence intervals are then calculated from the bootstrap distribution
  #> corresponding to the 2.5th and 97.5th percentiles.

#> These confidence intervals represent the confidence that we have that the observed variability metric is
#> representative of the population as a whole.

#In many of the later tests, the sample size is very low, therefore the confidence intervals are very wide.

#Sample sizes
p.aba %>%
  dplyr::select(dpi.f, primary_treatment)%>%
  tbl_summary(
    by=dpi.f
  )%>%
  modify_header(
    label ~ "**Variability Sample Size**"
  )

#### Define functions for calculating CV, PV, and V2 ####
source("Variability_Functions.R")

####Calculate Variability in Antibody Levels Primary####
# Set the number of bootstrap replicates
n_boot <- 1000

# Pipeline to calculate CV, PV, V2, and their 95% confidence intervals of antibody levels
# Using p.aba (**ALL individuals regardless of infection status**)
  # Change to p.abi to look at only individuals that became infected
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
  # Remove the bootstrap columns to clean up the dataset
  select(-cv_bootstrap, -pv_bootstrap, -v2_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.m <- left_join(p.aba, a.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.m$elisa_od <- p.aba.m$elisa_od.x
p.aba.m <- p.aba.m %>%
  select(-elisa_od.y, -elisa_od.x)

#just days samples were taken
p.aba.m <- p.aba.m %>%
  filter(dpi.f %in% c(-8, 14, 41))

#summary table
summary_tibble_priming_ab <- a.cv %>%
  group_by(dpi.f, primary_treatment) %>%
  summarize(
    CV = mean(bird_cv, na.rm = TRUE),    # Calculate the mean of bird_cv for each group
    SD = mean(bird_sd, na.rm = TRUE),    # Calculate the mean of bird_sd for each group
    PV = mean(bird_pv, na.rm = TRUE),    # Calculate the mean of bird_pv for each group
    elisa_od = mean(elisa_od, na.rm = TRUE),# Calculate the average elisa_od for each group
    lower_ci_pv = mean(pv_lower_ci),
    upper_ci_pv = mean(pv_upper_ci),
    lower_ci_cv = mean(cv_lower_ci),
    upper_ci_cv = mean(cv_upper_ci),
    lower_ci_v2 = mean(v2_lower_ci),
    upper_ci_v2 = mean(v2_upper_ci),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
#write.csv(summary_tibble_priming_ab, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/ELISA Variability Primary/ab_var_not_removed.csv", row.names=FALSE)

#overlay PV and elisa_od
#FIGURE 2 combined
fig2.comb <- ggplot(data = p.aba.m, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +
  # Original jittered data points for elisa_od
  geom_jitter(size = 1.5, alpha = 0.5, width = 0.15, height = 0, shape=16) +
  
  # Line and points with error bars for predicted values
  geom_path(data = dat.new, aes(x = dpi.f, y = yhat, group = primary_treatment,
                                color = primary_treatment), alpha=0.75) +
  geom_errorbar(data = dat.new, aes(ymin = Lower, ymax = Upper, x = dpi.f, y = yhat),
                color = "black", width = 0.01) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat, color = primary_treatment),
             size = 2, alpha = 0.75, shape = 16) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat),
             color = "black", size = 2, alpha = 1, shape = 1, stroke = 0.1) +
  
  # V2 (CV^2) points and confidence intervals
  geom_point(data = a.cv, aes(x = dpi.f, y = bird_pv, fill = primary_treatment, shape = "PV"),
             size = 3, alpha = 0.2) +
  geom_errorbar(data = a.cv, aes(x = dpi.f, ymin = pv_lower_ci, ymax = pv_upper_ci,
                                 group = primary_treatment, color = primary_treatment),
                width = 0.05, alpha = 0.02) +
  geom_path(data = a.cv, aes(x = dpi.f, y = bird_pv, color = primary_treatment,
                             group = primary_treatment), lty = "dashed", alpha=0.2) +
  
  # Labels and axis
  labs(y = "Anti-MG IgY Antibody Levels [OD]", 
       x = "Days Post Primary Inoculation", 
       shape = "Variability Metric", 
       color = "Primary Treatment", 
       fill = "Primary Treatment") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  # Custom shape scale for V2 (CV^2)
  scale_shape_manual(
    name = "Variability Metric",
    values = c("PV" = 17),  # Set shape type for CV^2
    labels = c("PV" = "PV")
  ) +
  
  scale_y_continuous(
    name = "Anti-MG IgY Antibody Levels [OD]", 
    sec.axis = sec_axis(~ ., name = "Variability in IgY Antibody Levels (PV)")
  ) +
  
  # Theme settings
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.text.y.right = element_text(color = "black", face = "bold", size = 15),
    axis.title.y.right = element_text(color = "black", size = 15, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

fig2.comb

#Figure 2 var only
fig2.pv<-ggplot(data = p.aba.m, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +
  # Original jittered data points for elisa_od
  #geom_jitter(aes(shape=primary_treatment), size = 1.5, alpha = 1, width = 0.15, height = 0) +
  
  # PV
  geom_point(data = summary_tibble_priming_ab, aes(x = dpi.f, y = PV, fill = primary_treatment, shape="PV"),
             size = 3, alpha = 1) +
  geom_errorbar(data = summary_tibble_priming_ab, aes(x = dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv,
                                 group = primary_treatment, color = primary_treatment),
                width = 0.05, alpha = 0.75) +
  geom_path(data = summary_tibble_priming_ab, aes(x = dpi.f, y = PV, color = primary_treatment,
                             group = primary_treatment), lty = "solid", alpha=0.75) +
  
  # CV
  geom_point(data = summary_tibble_priming_ab, aes(x = dpi.f, y = CV, fill = primary_treatment, shape="CV"),
             size = 3, alpha = 0.5) +
  geom_errorbar(data = summary_tibble_priming_ab, aes(x = dpi.f, ymin = lower_ci_cv, ymax = upper_ci_cv,
                                 group = primary_treatment, color = primary_treatment),
                width = 0.05, alpha = 0.5, lty="dashed") +
  geom_path(data = summary_tibble_priming_ab, aes(x = dpi.f, y = CV, color = primary_treatment,
                             group = primary_treatment), lty = "dashed", alpha=0.75) +
  
  # Labels and axis
  labs(y = "Antibody Variability", 
       x = "Days Post Primary Inoculation", 
       shape = "Variability Metric", 
       color = "Primary Treatment", 
       fill = "Primary Treatment",
       lty="Variability Metric") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  scale_shape_manual(
    name = "Variability Metric",
    values = c("PV" = 19, "CV" = 19),  # Set shape type
    labels = c("PV" = "PV", "CV" = "CV"))+
  
  # Theme settings
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 15),
    axis.text.y.right = element_text(color = "black", face = "bold", size = 15),
    axis.title.y.right = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

fig2.pv

####Do Individual Differences in Antibody Levels Predict Secondary Susceptibility####
####Analysis 2) Does individual variation in antibody response to infection predict reinfection probability?
#For this analysis, I am asking whether individual variation on any days that antibodies were measured predict whether that individual will become reinfected upon 
  #secondary challenge. 
#First, I look at whether individual variation in antibody levels (antibodies from now on) 8 days prior to initial infection are predictive
  #of reinfection with the prediction that they will not be. Because individuals were challenged with 5 different secondary doses[log10], secondary dose[log10] is used
  #as a fixed effect. This first baseline analysis shows the effect of secondary dose[log10] on reinfection probability, which is expected to increase as secondary
  #dose[log10] increases.
#Next, I look at whether antibodies on day 14 post-primary challenge are predictive of susceptibility.
#Finally, I look at whether antibodies on day 41 post-primary challenge are predictive of susceptibility. 
  #I look at these three days independently with three separate models with the objective of asking whether one measure of antibody levels at any point during priming
  #infection (before, during, or after) can be used to predict reinfection probability.

#Four birds were not recovered just prior to reinfection on DPPI 41. These will be omitted from this analysis.
#omit from dataset - were still infected (path load) on day 41 prior to reinfection day 42
unique(not_recovered$band_number)

#Do not omit 2449 - it did not surpass the 50 copy threshold and therefore is not considered infected
s.ab <- m.ab %>%
  filter(!(band_number %in% c(2274, 2469, 2520, 2494)))

s.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES SECONDARY**"
  )

s.ab.w <- s.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41))

wibird <- s.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)

#DPPI -8
wibird$elisa_od_pre <- wibird$`elisa_od_-8`

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


wibird %>%
  dplyr::select(log10.sec_dose, elisa_od_pre, inf_sec)%>%
  tbl_summary(
    by=inf_sec
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES Susceptibility**"
  )

#Does individual variation in antibody levels prior to primary challenge (DPPI -8) with MG predict susceptibility upon reinfection?
glm.ab.pre <- glm(inf_sec ~ elisa_od_pre + log10.sec_dose, data=wibird, family=binomial())
summary(glm.ab.pre)
car::Anova(glm.ab.pre, type = 3)
simulateResiduals(glm.ab.pre, plot=T)

# Get the summary of the model
summary_output <- summary(glm.ab.pre)
# Extract coefficients table
coef_table <- summary_output$coefficients

# Convert to a data frame for easier export
coef_df <- as.data.frame(coef_table)

# Create a summary data frame with additional model information
model_info <- data.frame(
  AIC = AIC(glm.ab.pre),
  Deviance = summary_output$deviance,
  Null_Deviance = summary_output$null.deviance,
  Residual_Df = summary_output$df.residual,
  Null_Df = summary_output$df.null
)

final_summary_pre <- cbind(coef_df, model_info)

#write.csv(final_summary_pre, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/glm_ab_pre_full_summary.csv", row.names = TRUE)

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
  geom_jitter(alpha=0.5, height=0.01, width=0)+
  geom_line(data=dat.new, aes(x=(elisa_od_pre), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_pre, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y=" ", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi ~log10.sec_dose, labeller = as_labeller(c(dpi_namespre, sec_dose_names)))
sus_pre_raw


##DPPI 14
#Fifteen plasma samples were lost from DPPI 14 sampling and one from DPPI 41 sampling, resulting in reduced sample sizes for antibody analyses.
  #Of the 15 samples lost from DPPI 14, 14 were sham inoculated and one was inoculated with a low dose. 
  #The one plasma sample lost from DPPI 41 was from a bird inoculated with a low dose.

#Which samples did not have plasma samples?
missing <- m.ab %>% 
  filter(dpi %in% c(-8, 14, 41) & is.na(elisa_od)) %>%
  dplyr::select(band_number, dpi, primary_treatment, secondary_dose, elisa_od)
missing

#omit individuals with no plasma sample on DPPI 14
wibird14 <- wibird %>%
  drop_na(elisa_od_14)

wibird14 %>%
  dplyr::select(log10.sec_dose, elisa_od_14, inf_sec)%>%
  tbl_summary(
    by=inf_sec
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES DPPI 14**"
  )

#Does individual variation in antibody levels at peak primary infection (DPPI 14) MG predict susceptibility upon reinfection?
glm.ab14 <- glm(inf_sec ~ elisa_od_14 + log10.sec_dose, data=wibird14, family=binomial())
summary(glm.ab14)

simulateResiduals(glm.ab14, plot=T)
car::Anova(glm.ab14, type = 3)

# Get the summary of the model
summary_output <- summary(glm.ab14)
# Extract coefficients table
coef_table <- summary_output$coefficients

# Convert to a data frame for easier export
coef_df <- as.data.frame(coef_table)

# Create a summary data frame with additional model information
model_info <- data.frame(
  AIC = AIC(glm.ab14),
  Deviance = summary_output$deviance,
  Null_Deviance = summary_output$null.deviance,
  Residual_Df = summary_output$df.residual,
  Null_Df = summary_output$df.null
)

final_summary_14 <- cbind(coef_df, model_info)

#write.csv(final_summary_14, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/glm_ab14_full_summary.csv", row.names = TRUE)

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
  geom_jitter(alpha=0.5, height=0.01, width=0)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(y= "Reinfection Status", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme(axis.title.x=element_blank())+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names14,sec_dose_names)))
sus_14_raw

#DPPI 41
#drop birds with no plasma sample
wibird41 <- wibird %>%
  drop_na(elisa_od_41)

wibird41 <- wibird41 %>% mutate(totalTES = sum())
wibird41$sec_dose.n <- wibird41$secondary_dose+1
wibird41$log10.sec_dose <- round(log10(wibird41$sec_dose.n), digits = 2)

wibird41$log10.sec_dose <- as.numeric(wibird41$log10.sec_dose)

wibird41 %>%
  dplyr::select(log10.sec_dose, elisa_od_pre, inf_sec)%>%
  tbl_summary(
    by=inf_sec
  )%>%
  modify_header(
    label ~ "**FINAL SAMPLE SIZES PRE**"
  )

#Does individual variation in antibody levels just before reinfection (DPPI 41) MG predict susceptibility upon reinfection?
glm.ab41 <- glm(inf_sec ~ elisa_od_41 + log10.sec_dose, data=wibird41, family=binomial())

summary(glm.ab41)
car::Anova(glm.ab41, type = 3)
simulateResiduals(glm.ab41, plot=T)
# Get the summary of the model
summary_output <- summary(glm.ab41)
# Extract coefficients table
coef_table <- summary_output$coefficients

# Convert to a data frame for easier export
coef_df <- as.data.frame(coef_table)

# Create a summary data frame with additional model information
model_info <- data.frame(
  AIC = AIC(glm.ab41),
  Deviance = summary_output$deviance,
  Null_Deviance = summary_output$null.deviance,
  Residual_Df = summary_output$df.residual,
  Null_Df = summary_output$df.null
)

final_summary_41 <- cbind(coef_df, model_info)

#write.csv(final_summary_41, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/glm_ab41_full_summary.csv", row.names = TRUE)


mod <- glm.ab41
dat.new=expand.grid(log10.sec_dose=unique(wibird41$log10.sec_dose),
                    inf_sec=unique(wibird41$inf_sec),
                    primary_treatment=unique(wibird41$primary_treatment),
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
  geom_jitter(alpha=0.5, height=0.01, width=0)+
  geom_line(data=dat.new, aes(x=(elisa_od_41), y=yhat), alpha = 1, size =1)+
  geom_ribbon(data = dat.new, aes(ymin=Lower, ymax = Upper, x= elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7000"))+
  labs(x="IgY Antibody Levels [OD]", y= "Reinfection Status", color="Secondary Dose", fill ="Secondary Dose")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names41,sec_dose_names)))
sus_41_raw


#Compare Models
aictab(cand.set=list(glm.ab.pre, glm.ab14, glm.ab41), modnames=c("glm.ab.pre","glm.ab14", "glm.ab41"))


#ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Figures/sus_ab_14.png",
#       plot=sus_ab_14, width=8, height=4, dpi = 300)



# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Figures/sus_ab_41.png",
#        plot=sus_ab_41, width=8, height=4, dpi = 300)

#Figure 3
Fig3 <- sus_14_raw / sus_41_raw

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Figures/Fig3.png",
#         plot=Fig3, width=10, height=5, dpi = 300)

####Variability With Hawley et al., 2024####
#extract variability metrics from summary_tibble_priming_ab a.cv calculatons
summary_tibble_priming_ab %>% filter(dpi.f == 41)

dose <- c("Sham", "Low", "High")
CV <- as.numeric(c(0.899, 1.630, 2.511)) #Susceptibility CV Hawley et al., 2024 (p5)
ab.pv <- as.numeric(c(0.0592, 0.126, 0.187))
ab.cv <- as.numeric(c(0.0710, 0.153, 0.261))
ci_lower_cv <- as.numeric(c(0.0377, 0.106, 0.170))
ci_upper_cv <- as.numeric(c(0.09970471, 0.18774258, 0.33572240))
ci_lower_pv <- as.numeric(c(0.04002442, 0.09085289, 0.14255949))
ci_upper_pv <- as.numeric(c(0.07988264, 0.15355463, 0.22592119))
mean_sus <- as.numeric(c(1.181, 0.446, 0.192))
mean_ab <- as.numeric(c(0.04506522, 0.04812000, 0.05767925))
dpi <- 41

het.df <- data.frame(dose, CV, mean_sus, ab.pv, ab.cv, dpi, mean_ab, ci_lower_cv, ci_upper_cv, ci_lower_pv, ci_upper_pv)

#hawley et al 2024
het.paper.only <- ggplot(het.df, aes(x=Metric))+
  geom_point(size=5,aes(x="CV", y=CV, color=fct_rev(dose)), shape = 17)+
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
  geom_point(size=5,aes(x="PV", y=ab.pv, color=fct_rev(dose), shape="PV"))+
  geom_errorbar(aes(ymin=ci_lower_pv, ymax=ci_upper_pv, y=ab.pv, x="PV", color=dose), width=0.05)+
  geom_point(size=5,aes(x="CV", y=ab.cv, color=dose, shape = "CV"))+
  geom_errorbar(aes(ymin=ci_lower_cv, ymax=ci_upper_cv, y=ab.cv, x="CV", color=dose), width=0.05)+
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
    legend.position="right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.background= element_rect(size=0.25, linetype="solid"))+
  labs(color = "Primary Treatment", shape= "Variability Metric")
ab.paper.only


# Remove legends from plots
het.paper.only <- het.paper.only + theme(legend.position = "none")
ab.paper.only <- ab.paper.only + theme(legend.position = "right")

#Figure 4
het.paper.only + ab.paper.only


####Multilevel hierarchical modeling: Are antibody levels dppi 41 correlated with heterogeneity in susceptibility?
mlm41 <- wibird41 %>% 
  dplyr::select(band_number, primary_treatment, sex, elisa_od_41)

mlm41$primary_treatment <- as.factor(mlm41$primary_treatment)

het.df$primary_treatment <- as.factor(het.df$dose)
het.df$primary_treatment <- factor(het.df$primary_treatment)

mlm <- merge(mlm41, het.df, by="primary_treatment", all.x=TRUE)
mlm$sus_cv <- mlm$CV

#model: antibody levels predictor primary_treatment is radom effect as each group has one CV value.
hmod <- lm(sus_cv ~ elisa_od_41, data= mlm)
plot(allEffects(hmod))

ggplot(mlm, aes(x=sus_cv, y=elisa_od_41, color=primary_treatment))+
  geom_point()

#simple model - is heterogeneity in susceptibility correlated with antibody variability dppi 41?
mod_group <- lm(sus_cv ~ ab.pv + primary_treatment, data=mlm)
summary(mod_group)
plot(allEffects(mod_group))

ggplot(mlm, aes(x=sus_cv, y=ab.pv, color=primary_treatment))+
  geom_point()

ggplot(het.df, aes(x = ab.pv, y = CV, color=primary_treatment)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = ci_lower_pv, xmax= ci_upper_pv, x=ab.pv), width=0.01)+
  #geom_point(aes(x=ab.cv, y=CV, color=primary_treatment), size=2, alpha=0.5)+
  #geom_errorbar(aes(xmin = ci_lower_cv, xmax= ci_upper_cv, x=ab.cv), width=0.01, lty="dashed")+
  #geom_errorbar(aes(ymin = ))
  #geom_text(aes(label = primary_treatment), vjust = -1) +
  scale_color_manual(values=pri_colors)+
  labs(x = "Antibody Variability (PV)", 
       y = "Variability in Susceptibility (CV) Hawley et al., 2024", 
       title = "Comparison of Variability in Antibodies and Susceptibility",
       color="Primary Treatment") +
  theme_minimal()

ggplot(het.df, aes(x = ab.pv, y = ab.cv, color=primary_treatment)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = ci_lower_pv, xmax= ci_upper_pv, x=ab.pv), width=0.001)+
  geom_errorbar(aes(ymin = ci_lower_cv, ymax= ci_upper_cv, y=ab.cv), width=0.001)+
  #geom_errorbar(aes(ymin = ))
  #geom_text(aes(label = primary_treatment), vjust = -1) +
  labs(x = "Antibody Variability (PV)", 
       y = "Antibody Variability (CV)", 
       title = "Comparison of Variability in Antibodies",
       color = "Primary Treatment") +
  scale_color_manual(values=c(pri_colors))+
  theme_minimal()

####Secondary Challenge####
#df with only birds high secondary dose after reinfection
sec.ab <- m.ab %>%
  filter(dpi > 41 )

#remove birds that did not recover
unique(not_recovered$band_number)

sec.ab <- sec.ab %>% 
  filter(!band_number %in% c(2274, 2469, 2494, 2520))

#Sample sizes 
sec.ab %>%
  #filter(dpi == 56) %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, inf_sec) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary(
  by = dpi
  )%>%
  modify_header(
    label ~ "**By Primary and Secondary Dose**",
  )

#missing swab sample from 2426 dppi 46
who <- sec.ab %>% filter(band_number == 2426) %>%
  dplyr::select(date, band_number, dpi, quantity, tes, inf_pri, inf_sec)
#for comparison; all secondary doses
m.ab %>%
  filter(dpi > 41) %>%
  dplyr::select(dpi, inf_sec, primary_treatment, secondary_dose) %>%
  group_by(dpi, inf_sec, secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = dpi # Grouping by infected secondary
    #statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

sec.ab %>%
  filter(dpi == 56) %>%
  dplyr::select(inf_sec, primary_treatment, secondary_dose, band_number) %>%
  group_by(inf_sec, secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = inf_sec, # Grouping by infected secondary
    statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

#Number infected secondary 7000
sec.ab %>%
  filter(dpi == 56 & secondary_dose == "7000") %>%
  dplyr::select(inf_sec, primary_treatment, secondary_dose, band_number) %>%
  group_by(inf_sec, secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = inf_sec, # Grouping by infected secondary
    statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

#2449 is an issue
ggplot(sec.ab %>% filter(band_number == 2449), aes(x=dpi, y=tes, color=secondary_dose))+
  geom_point()

##We want to ask: Of the individuals that became reinfected, how variable were their responses to infection?
  #Therefore, filter so that inf_sec == 1
# sec.ab <- sec.ab %>%
#   filter(inf_sec == 1)

#Updated: individuals who were not susceptible are an important part of heterogeneity in disease responses
  #Therefore, they will be left in > any inf_sec

#Calculate max eye score and pathogen load
max.quant.sec <- sec.ab %>%
  group_by(band_number, primary_treatment, secondary_dose) %>%
  reframe(max_quantity = max(quantity, na.rm = TRUE))

max.quant.sec %>%
  filter(secondary_dose == "7000") %>%
  dplyr::select(primary_treatment, secondary_dose, band_number) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = primary_treatment, # Grouping by infected secondary
    statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

max.tes.sec <- sec.ab %>%
  group_by(band_number, primary_treatment, secondary_dose)%>%
  reframe(max_tes = max(tes, na.rm=TRUE))

max.tes.sec %>%
  filter(secondary_dose == "7000") %>%
  dplyr::select(primary_treatment, secondary_dose, band_number) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = primary_treatment, # Grouping by infected secondary
    statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

max.sec<-merge(max.quant.sec, max.tes.sec, by="band_number")

max.sec$primary_treatment <- max.sec$primary_treatment.x
max.sec$secondary_dose <- max.sec$secondary_dose.x

max.sec <- max.sec %>%
  dplyr::select(-c(primary_treatment.y, secondary_dose.y, primary_treatment.x, secondary_dose.x))
max.sec$max_quantity1 <- max.sec$max_quantity+1

max.sec %>%
  filter(secondary_dose == "7000") %>%
  dplyr::select(primary_treatment, secondary_dose, band_number) %>%
  group_by(secondary_dose, primary_treatment) %>%
  tbl_summary(
    by = primary_treatment, # Grouping by infected secondary
    statistic = list(band_number ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 

max.sec <- max.sec %>%
  filter(secondary_dose == "7000")

#Visualize TES Secondary
ggplot(sec.ab, aes(x=dpi, y=tes, color=primary_treatment))+
  geom_point()+
  geom_line(aes(group= band_number), alpha=0.5)+
  stat_summary(aes(group= primary_treatment), geom="line", fun="mean", size=1)+
  scale_color_manual(values=pri_colors)

#Visualize max TES Secondary
tes.max.s <- ggplot(max.sec, aes(x=primary_treatment, y=max_tes, color=primary_treatment))+
  geom_boxplot(outlier.shape = 17, width=0.5)+
  #geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5, binwidth = 0.1, stackratio = 1, alpha=0.75, aes(fill=primary_treatment))+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  labs(x="Primary Treatment", y="Maximum Eyescore (7000 Dose)", color="Primary Treatment", fill= "Primary Treatment")
tes.max.s  



#Visualize Quantity1 Secondary
ggplot(sec.ab, aes(x=dpi, y=quantity1, color=primary_treatment))+
  geom_hline(yintercept = 50, alpha=0.5, lty="dashed")+
  geom_jitter(width=0, height = 0, alpha=0.5)+
  geom_line(aes(group= band_number), alpha=0.5)+
  stat_summary(aes(group= primary_treatment), geom="line", fun="mean", linewidth=1)+
  #scale_y_log10()+
  scale_color_manual(values=pri_colors)


quant.max.s.log10 <- ggplot(max.sec, aes(x=primary_treatment, y=max_quantity1, color=primary_treatment))+
  geom_hline(yintercept = 51, lty="dashed", alpha=0.5)+
  geom_boxplot(outlier.shape = 17, width=0.5)+
  #geom_dotplot(binaxis="y", stackdir="center", dotsize=0.1, binwidth = 1, stackratio = 1, alpha=0.75, aes(fill=primary_treatment))+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  scale_y_log10()+
  labs(x="Primary Treatment", y="Log10 Maximum Quantity (7000 Dose)", color="Primary Treatment", fill= "Primary Treatment")
quant.max.s.log10


#Variability secondary
max.sec$max_tes <- max.sec$max_tes+0.001
#max_quantity1 is max_quantity + 1 <- add 1 because calculations do not work with 0s
  #log10(1) = 0 which does not work for variability calculations so add small constant
max.sec$lmax_quantity <- log10(max.sec$max_quantity+1)+0.001

max.s.v <- max.sec %>% 
  group_by(primary_treatment) %>%
  reframe(
    # Metrics for max_tes
    max_tes = max_tes,
    mean_tes = mean(max_tes-0.001), #to get actual mean
    bird_cv_tes = calculate_cv(max_tes),
    bird_sd_tes = sd(max_tes - 0.001),
    bird_pv_tes = calculate_pv(max_tes),
    bird_v2_tes = calculate_v2(max_tes),
    bird_se_tes = SE(max_tes-0.001), #to get actual SE
    
    # Metrics for max_quantity
    max_quantity = max_quantity1,
    mean_quantity = mean(max_quantity1),
    bird_cv_quantity = calculate_cv(max_quantity1),
    bird_sd_quantity = sd(max_quantity1),
    bird_pv_quantity = calculate_pv(max_quantity1),
    bird_v2_quantity = calculate_v2(max_quantity1),
    bird_se_quantity = SE(max_quantity1),
    
    # Metrics for log10(max_quantity)
    lmax_quantity = lmax_quantity,
    mean_lmax_quantity = mean(lmax_quantity),
    bird_cv_lmax_quantity = calculate_cv(lmax_quantity),
    bird_sd_lmax_quantity = sd(lmax_quantity),
    bird_pv_lmax_quantity = calculate_pv(lmax_quantity),
    bird_v2_lmax_quantity = calculate_v2(lmax_quantity),
    bird_se_lmax_quantity = SE(lmax_quantity),
    
    #Other info
    band_number = band_number,
    primary_treatment = primary_treatment,
    
    # Bootstrap for confidence intervals (max_tes)
    cv_bootstrap_tes = list(replicate(n_boot, calculate_cv(sample(max_tes, replace = TRUE)))),
    v2_bootstrap_tes = list(replicate(n_boot, calculate_v2(sample(max_tes, replace = TRUE)))),
    pv_bootstrap_tes = list(replicate(n_boot, calculate_pv(sample(max_tes, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (max_quantity)
    cv_bootstrap_quantity = list(replicate(n_boot, calculate_cv(sample(max_quantity, replace = TRUE)))),
    v2_bootstrap_quantity = list(replicate(n_boot, calculate_v2(sample(max_quantity, replace = TRUE)))),
    pv_bootstrap_quantity = list(replicate(n_boot, calculate_pv(sample(max_quantity, replace = TRUE)))),

  # Bootstrap for confidence intervals (lmax_quantity)
  cv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_cv(sample(lmax_quantity, replace = TRUE)))),
  v2_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_v2(sample(lmax_quantity, replace = TRUE)))),
  pv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_pv(sample(lmax_quantity, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for both max_tes and max_quantity
  mutate(
    # max_tes CIs
    cv_lower_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.025)),
    cv_upper_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.975)),
    v2_lower_ci_tes = map_dbl(v2_bootstrap_tes, ~ quantile(.x, 0.025)),
    v2_upper_ci_tes = map_dbl(v2_bootstrap_tes, ~ quantile(.x, 0.975)),
    pv_lower_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.025)),
    pv_upper_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.975)),
    
    # max_quantity CIs
    cv_lower_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    v2_lower_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.025)),
    v2_upper_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.975)),

    # lmax_quantity CIs
    cv_lower_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975)),
    v2_lower_ci_lmax_quantity = map_dbl(v2_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    v2_upper_ci_lmax_quantity = map_dbl(v2_bootstrap_lmax_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  # select(-cv_bootstrap_tes, -v2_bootstrap_tes, -pv_bootstrap_tes,
  #        -cv_bootstrap_quantity, -v2_bootstrap_quantity, -pv_bootstrap_quantity) %>%
  ungroup()

# Add to df and format names
max.sec.v <- left_join(max.sec, max.s.v, by = c("primary_treatment", "band_number"))

# Summary table for both max_tes and max_quantity
summary_tibble.s.v <- max.s.v %>%
  group_by(primary_treatment) %>%
  summarize(
    # Summary statistics for max_tes
    CV_tes = mean(bird_cv_tes, na.rm = TRUE),
    SD_tes = mean(bird_sd_tes, na.rm = TRUE),
    V2_tes = mean(bird_v2_tes, na.rm = TRUE),
    PV_tes = mean(bird_pv_tes, na.rm = TRUE),
    SE_tes = mean(bird_se_tes),
    
    lower_ci_pv_tes = mean(pv_lower_ci_tes),
    upper_ci_pv_tes = mean(pv_upper_ci_tes),
    lower_ci_cv_tes = mean(cv_lower_ci_tes),
    upper_ci_cv_tes = mean(cv_upper_ci_tes),
    lower_ci_v2_tes = mean(v2_lower_ci_tes),
    upper_ci_v2_tes = mean(v2_upper_ci_tes),
    
    # Summary statistics for max_quantity
    CV_quantity = mean(bird_cv_quantity, na.rm = TRUE),
    SD_quantity = mean(bird_sd_quantity, na.rm = TRUE),
    V2_quantity = mean(bird_v2_quantity, na.rm = TRUE),
    PV_quantity = mean(bird_pv_quantity, na.rm = TRUE),
    SE_quantity = mean(bird_se_quantity),
    
    lower_ci_pv_quantity = mean(pv_lower_ci_quantity),
    upper_ci_pv_quantity = mean(pv_upper_ci_quantity),
    lower_ci_cv_quantity = mean(cv_lower_ci_quantity),
    upper_ci_cv_quantity = mean(cv_upper_ci_quantity),
    lower_ci_v2_quantity = mean(v2_lower_ci_quantity),
    upper_ci_v2_quantity = mean(v2_upper_ci_quantity),
    
    # Summary statistics for lmax_quantity
    CV_lmax_quantity = mean(bird_cv_lmax_quantity, na.rm = TRUE),
    SD_lmax_quantity = mean(bird_sd_lmax_quantity, na.rm = TRUE),
    V2_lmax_quantity = mean(bird_v2_lmax_quantity, na.rm = TRUE),
    PV_lmax_quantity = mean(bird_pv_lmax_quantity, na.rm = TRUE),
    SE_lmax_quantity = mean(bird_se_lmax_quantity),
    
    lower_ci_pv_lmax_quantity = mean(pv_lower_ci_lmax_quantity),
    upper_ci_pv_lmax_quantity = mean(pv_upper_ci_lmax_quantity),
    lower_ci_cv_lmax_quantity = mean(cv_lower_ci_lmax_quantity),
    upper_ci_cv_lmax_quantity = mean(cv_upper_ci_lmax_quantity),
    lower_ci_v2_lmax_quantity = mean(v2_lower_ci_lmax_quantity),
    upper_ci_v2_lmax_quantity = mean(v2_upper_ci_lmax_quantity),
    
    # Number of individuals
    n_individuals = n_distinct(band_number),
    mean_tes = mean(mean_tes),
    mean_quantity = mean(mean_quantity),
    mean_lmax_quantity = mean(mean_lmax_quantity)-0.001
  ) %>%
  as_tibble()

print(summary_tibble.s.v)

s.var <- summary_tibble.s.v %>%
  dplyr::select(primary_treatment, CV_tes, PV_tes, mean_tes, SE_tes, 
                CV_quantity, PV_quantity, mean_quantity, SE_quantity, 
                CV_lmax_quantity, PV_lmax_quantity, mean_lmax_quantity, SE_lmax_quantity, 
                n_individuals)

#Visualize difference in variability metrics Max Eye Score Secondary
max_sec_tes.v<- ggplot(summary_tibble.s.v) +
  # geom_point(aes(x = "CV", y = CV_tes, color = primary_treatment), position = position_dodge(width = 0.5)) +
  # geom_errorbar(aes(x = "CV", ymin = lower_ci_cv_tes, ymax = upper_ci_cv_tes, y = CV_tes, color = primary_treatment), 
  #               width = 0., position = position_dodge(width = 0.5)) +
  
  geom_point(aes(x = "PV", y = PV_tes, color = primary_treatment), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = "PV", ymin = lower_ci_pv_tes, ymax = upper_ci_pv_tes, y = PV_tes, color = primary_treatment), 
                width = 0., position = position_dodge(width = 0.5)) +
  
  geom_point(aes(x = "V2", y = V2_tes, color = primary_treatment), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = "V2", ymin = lower_ci_v2_tes, ymax = upper_ci_v2_tes, y = V2_tes, color = primary_treatment), 
                width = 0., position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = pri_colors) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("CV" = 17, "PV" = 19, "V2" = 18),
                     labels = c("CV" = "CV", "PV" = "PV", "V2" = "CV²")) +
  labs(title = "Secondary Eye Score Variability",
       y = "Variability",
       x = "Metric",
       color = "Primary Treatment") +
  scale_y_continuous(limits = c(0,1))+
  theme_bw()+
  theme(legend.position = "none")


#Visualize difference in variability metrics Max Pathogen Load Secondary
max_sec_path.v <- ggplot(summary_tibble.s.v) +
  # geom_point(aes(x = "CV", y = CV_quantity, color = primary_treatment), position = position_dodge(width = 0.5)) +
  # geom_errorbar(aes(x = "CV", ymin = lower_ci_cv_quantity, ymax = upper_ci_cv_quantity, y = CV_quantity, color = primary_treatment), 
  #               width = 0., position = position_dodge(width = 0.5)) +
  
  geom_point(aes(x = "PV", y = PV_quantity, color = primary_treatment), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = "PV", ymin = lower_ci_pv_quantity, ymax = upper_ci_pv_quantity, y = PV_quantity, color = primary_treatment), 
                width = 0., position = position_dodge(width = 0.5)) +
  
  geom_point(aes(x = "V2", y = V2_quantity, color = primary_treatment), position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = "V2", ymin = lower_ci_v2_quantity, ymax = upper_ci_v2_quantity, y = V2_quantity, color = primary_treatment), 
                width = 0., position = position_dodge(width = 0.5)) +
  
  scale_color_manual(values = pri_colors) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("CV" = 17, "PV" = 19, "V2" = 18),
                     labels = c("CV" = "CV", "PV" = "PV", "V2" = "V2")) +
  labs(title = "Secondary Pathogen Load Variability",
       y = "Variability",
       x = "Metric",
       color = "Primary Treatment") +
  scale_y_continuous(limits = c(0,1))+
  theme_bw()

#Variability in max tes and pathogen load secondary
max_sec_tes.v+ max_sec_path.v

####Overlay Raw values and PV of Quantity and Eye Score Primary and Secondary####

#Eye Score Secondary
#Visualize max TES Secondary
tes.max.s <- ggplot(max.sec, aes(x=primary_treatment, y=max_tes, color=primary_treatment))+
  geom_boxplot(outlier.shape = 17, width=0.5)+
  #geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5, binwidth = 0.1, stackratio = 1, alpha=0.75, aes(fill=primary_treatment))+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  labs(x="Primary Treatment", y="Maximum Eyescore (7000 Dose)", color="Primary Treatment", fill= "Primary Treatment")
tes.max.s  

# Define a scale transformation factor
  #Scale factor puts PV and CI on a scale between 0 and the max TES
scale_factor.ts <- max(max.tes.sec$max_tes, na.rm = TRUE) / max(summary_tibble.s.v$PV_tes, na.rm = TRUE)

#this scale factor scales PV and CI between 0 and 6 keeping PV on a scale of 0 to 1
scale_factor.es <- 6 / 1

tes.max.sec.comb <- ggplot() +
  # Primary y-axis: max TES
  geom_boxplot(data=max.tes.sec %>% filter(secondary_dose == "7000"), aes(x=primary_treatment, y=max_tes, color=primary_treatment), outlier.shape = 8, width=0.25) +
  geom_jitter(data=max.tes.sec %>% filter(secondary_dose == "7000"), aes(x=primary_treatment, y=max_tes, color=primary_treatment), width=0.25, height=0, alpha=0.5) +
  # Secondary y-axis: PV TES (rescaled)
  geom_point(data=summary_tibble.s.v, aes(x=primary_treatment, y=PV_tes * scale_factor.es, color=primary_treatment), shape = 17, size=3) +
  geom_point(data=summary_tibble.s.v, aes(x=primary_treatment, y=PV_tes * scale_factor.es), shape = 24, size=3) +
  geom_errorbar(data=summary_tibble.s.v,
                aes(x=primary_treatment, y=PV_tes * scale_factor.es, ymin = lower_ci_pv_tes * scale_factor.es,
                    ymax = upper_ci_pv_tes * scale_factor.ts), width=0.15, lty="dashed") +
  
  # Manual color and fill
  scale_color_manual(values=pri_colors) +
  scale_fill_manual(values=pri_colors) +
  
  # Adjust y-axis labels with secondary y-axis
  scale_y_continuous(
    name = "Maximum Eyescore Secondary",
    sec.axis = sec_axis(~ . / scale_factor.es, name = "PV Total Eye Score")
  ) +
  
  labs(x="Primary Treatment", color="Primary Treatment", fill="Primary Treatment")

tes.max.sec.comb


#Path Load Secondary
#log10(quantity)
# Define a scale transformation factor
#scale_factor.lqs <- max(log10(max.quant.sec$max_quantity+0.001), na.rm = TRUE) / max(summary_tibble.s.v$PV_lmax_quantity, na.rm = TRUE)
max(log10(max.quant.sec$max_quantity))
scale_factor.lpl <- 6 / 1



quant.lmax.sec.comb <- ggplot() +
  # Primary y-axis: max quantity
  geom_hline(yintercept = log10(51), lty="dashed", alpha=0.5)+
  geom_boxplot(data=max.quant.sec %>% filter(secondary_dose == "7000"), aes(x=primary_treatment, y=log10(max_quantity+1), color=primary_treatment), outlier.shape = 8, width=0.25) +
  geom_jitter(data=max.quant.sec %>% filter(secondary_dose == "7000"), aes(x=primary_treatment, y=log10(max_quantity+1), color=primary_treatment), width=0.25, height=0, alpha=0.5) +
  
  # Secondary y-axis: PV TES (rescaled)
  geom_errorbar(data=summary_tibble.s.v,
                aes(x=primary_treatment, y=PV_lmax_quantity * scale_factor.lpl, ymin = lower_ci_pv_lmax_quantity * scale_factor.lpl,
                    ymax = upper_ci_pv_lmax_quantity * scale_factor.lpl), width=0.15, lty="dashed") +
  geom_point(data=summary_tibble.s.v, aes(x=primary_treatment, y=PV_lmax_quantity * scale_factor.lpl, color=primary_treatment),
             shape = 17, size=3, alpha=1) +
  geom_point(data=summary_tibble.s.v, aes(x=primary_treatment, y=PV_lmax_quantity * scale_factor.lpl), shape = 24, size=3) +
 
  
  # Manual color and fill
  scale_color_manual(values=pri_colors) +
  scale_fill_manual(values=pri_colors) +
  
  # Adjust y-axis labels with secondary y-axis
  scale_y_continuous(
    name = "Log10(Maximum Pathogen Load) Secondary",
    sec.axis = sec_axis(~ . / scale_factor.lpl, name = "PV Pathogen Load")
  ) +
  labs(x="Primary Treatment", color="Primary Treatment", fill="Primary Treatment")

quant.lmax.sec.comb
# #Figure 6
# quant.max.pri.comb / quant.max.sec.comb


#Alternatively, just secondary
tes.max.sec.comb / quant.lmax.sec.comb 

final_plot <- tes.max.sec.comb / quant.lmax.sec.comb

# Print the final plot
print(final_plot)

#Table 2; Summary of results for eyescore and pathogen load for secondary variability
s.var
#write_xlsx(s.var, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/secondary_variability.xlsx")

####Does primary treatment affect maximum eyescore or pathogen load upon secondary inoculation?####
kruskal_tes<- kruskal.test(max_tes ~ primary_treatment, data=max.sec)
kruskal_tes
#X2 = 19.06, df = 2, p < 0.0001

dunn_tes <- dunn.test(max.sec$max_tes, max.sec$primary_treatment, method = "bonferroni")
dunn_tes

dunn_tes_df <- tibble(
  Comparison = dunn_tes$comparisons,
  Z_value = dunn_tes$Z,
  P_unadjusted = dunn_tes$P,
  P_adjusted = dunn_tes$P.adjusted
)

kruskal_tes_summary <- tibble(
  statistic = kruskal_tes$statistic,
  df = kruskal_tes$parameter,
  p_value = kruskal_tes$p.value
)




#Kruskal Wallis Tests
kruskal_quant<- kruskal.test(lmax_quantity ~ primary_treatment, data=max.sec)
kruskal_quant
#X2 = 13.604, df = 2, p = 0.00111

dunn_quant <- dunn.test(max.sec$lmax_quantity, max.sec$primary_treatment, method = "bonferroni")
dunn_quant

dunn_quant_df <- tibble(
  Comparison = dunn_quant$comparisons,
  Z_value = dunn_quant$Z,
  P_unadjusted = dunn_quant$P,
  P_adjusted = dunn_quant$P.adjusted
)

kruskal_quant_summary <- tibble(
  statistic = kruskal_quant$statistic,
  df = kruskal_quant$parameter,
  p_value = kruskal_quant$p.value
)

# Export multiple sheets
# write_xlsx(
#   list(KruskalWallis_TES = kruskal_tes_summary, DunnTest_TES = dunn_tes_df,
#      KruskalWallis_lQuant = kruskal_quant_summary, DunnTest_lQuant = dunn_quant_df),
#  "EEID_1A_Antibody_Analysis_files/Results/combined_tes_lquant_sec_kruskal_dunn.xlsx"
# )

# Export multiple sheets
# write_xlsx(
#   list(KruskalWallis = kruskal_tes_summary, DunnTest = dunn_tes_df),
#  "EEID_1A_Antibody_Analysis_files/Results/combined_quant_sec_kruskal_dunn.xlsx"
# )

       