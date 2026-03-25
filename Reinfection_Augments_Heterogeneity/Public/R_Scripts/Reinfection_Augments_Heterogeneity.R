####Does prior pathogen exposure augment population-level heterogeneity in transmission-relevant traits?####
##Final Analysis; Peer Review
## Jun 2025; Updated Jan 2026
## JGL
#Requires: reinfection_response.rds, Variability_Functions.R, CIs_Sus_4042.R

####Load Packages + Read in Data####
rm(list=ls())

# Automatically install & load if missing
if (!require("pacman")) install.packages("pacman")
library(pacman)

# Core tidyverse
p_load(tidyverse, patchwork, broom.mixed, purrr)

# Modeling
p_load(glmmTMB, emmeans, AICcmodavg, VGAM, MASS)

# Diagnostics
p_load(DHARMa, car, effects)

# Visualization extras
p_load(ggh4x, grid, gridExtra)

# Statistics / post-hoc tests
p_load(dunn.test, brglm2, CValternatives)

# Reporting & tables
p_load(gtsummary, writexl)

#import data
#m.ab = merged + formatted working df
m.ab <- readRDS("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Final Dataframes/reinfection_response.rds")
#Infected column: inf = if quantity (where any NA is replaced with 0) is greater than cutoff of 50 copies return 1, else return 0

####Format theme####
#set colors for primary treatment
pri_colors <- c("#1B9E77", "#7570B3", "#D95F02")
#set color scheme secondary treatment
sec_colors <- c("#8C754B", "#77AB59", "#59A5D8", "#9F77D9", "#FA8072")
#set theme
theme_set(
  theme_bw() +
    theme(
      axis.title.y = element_text(color = "black", size = 15),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.text.x = element_text(color = "black", size = 15),
      legend.background = element_rect(linewidth = 0.25, linetype = "solid"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size=15),
      legend.text = element_text(size=15)
    )
)


#Sample Sizes; Table S4
m.ab %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose)%>%
  tbl_summary(
    by=secondary_dose
  )%>%
  modify_header(
    label ~ "**Starting Sample Size**"
  )

####Analysis 1) Do antibody levels vary across dpi and primary treatment?####
#>Fifteen plasma samples were lost from DPPI 14 sampling and one from DPPI 41 sampling, resulting in reduced sample sizes
#> for antibody analyses. Of the 15 samples lost from DPPI 14, 14 were sham inoculated and one was inoculated with a low dose. 
#> The one plasma sample lost from DPPI 41 was from a bird inoculated with a low dose.

#data frame with only primary infection
p.ab <- m.ab %>%
  filter(dpi <=41)

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

#Do antibody levels differ across sampling days and treatment?
#I treat days post inoculation as a factor in these models and use band_number as a random effect
#Antibodies only measured dpi -8, 14, and 41
p.abt <- p.ab %>%
  filter(dpi %in% c(-8, 14, 41))
p.abt$dpi.f <- as.factor(p.abt$dpi)

####Analysis #1: Do priming treatment or DPI predict antibody responses while controlling for individual bird ID?####
#model comparison
p1 <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p2 <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p3 <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), data=p.abt, family=Gamma(log))
p4<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p5<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
p6 <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, family=Gamma(log))
#Add dispersion formula to account for gamma model allowing dispersion parameter to vary by dpi
p1d <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p2d <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p3d <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p4d<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p5d<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), dispformula = ~dpi.f, data=p.abt, family=Gamma(log))
p6d <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, dispformula = ~dpi.f, family=Gamma(log))
#Add dispersion formula to account for gamma model allowing dispersion parameter to vary by dpi and primary treatment
p1dd <- glmmTMB(elisa_od~primary_treatment + dpi.f + (1|band_number), dispformula = ~primary_treatment*dpi.f, data=p.abt, family=Gamma(log))
p2dd <- glmmTMB(elisa_od~primary_treatment*dpi.f + (1|band_number), dispformula = ~primary_treatment*dpi.f, data=p.abt, family=Gamma(log))
p3dd <- glmmTMB(elisa_od~primary_treatment*dpi.f + sex + (1|band_number), dispformula = ~primary_treatment*dpi.f, data=p.abt, family=Gamma(log))
p4dd<- glmmTMB(elisa_od~primary_treatment + sex + dpi.f + (1|band_number), dispformula = ~primary_treatment*dpi.f, data=p.abt, family=Gamma(log))
p5dd<- glmmTMB(elisa_od~primary_treatment * sex + dpi.f + (1|band_number), dispformula = ~primary_treatment*dpi.f, data=p.abt, family=Gamma(log))
p6dd <- glmmTMB(elisa_od~1 + dpi.f + (1|band_number), data=p.abt, dispformula = ~primary_treatment*dpi.f, family=Gamma(log))

#Compare models with AICc
aictab(cand.set=list(p1, p2, p3, p4, p5, p6, p1d, p2d, p3d, p4d, p5d, p6d, p1dd, p2dd, p3dd, p4dd, p5dd, p6dd), 
       modnames=c("p1", "p2", "p3",  "p4", "p5", "p6", "p1d", "p2d", "p3d",  "p4d", "p5d", "p6d", "p1dd", "p2dd", "p3dd",  "p4dd", "p5dd", "p6dd"))

#sex not significant
summary(p3dd)

#Final Model; Does the interaction between priming treatment and dpi predict antibody response while controlling for individual bird id?
#Dispersion formula to account for gamma model allowing dispersion parameter to vary by dpi and primary treatment
lm1 <- glmmTMB(elisa_od~ primary_treatment*dpi.f + (1|band_number),
               dispformula = ~primary_treatment*dpi.f, 
               family=Gamma(link="log"),
               data=p.abt)

plot(residuals(lm1, type="pearson"))
simulateResiduals(lm1, plot=T)
summary(lm1)
car::Anova(lm1, type = "III")

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
fixed_effects <- tidy(lm1, effects = "fixed") #Extract fixed effects


#Post-Hoc Tukey
emm_results <- emmeans(lm1, ~ primary_treatment|dpi.f, scale="response")
pairwise <- pairs(emm_results, adjust="tukey")
(summary(pairwise))

# Convert pairwise comparisons to a data frame
pairwise_df <- as.data.frame(summary(pairwise))

# Table S1
# write_xlsx(
#   list(
#     Fixed_Effects = fixed_effects,
#     Dispersion_Model = dispersion_effects,
#     Random_Effects = random_effects,
#     Pairwise_Comparisons = pairwise_df
#   ),
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Supplementary/TableS1.xlsx"
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


####Variability in Antibody Levels across primary####
#add small constant to total eye score for calculation of variability
p.ab$tes.new <- p.ab$tes + .01
p.ab$dpi.f <- as.factor(p.ab$dpi)

#df with all individuals during primary challenge regardless of reinfection; p.aba = p.ab all
p.aba <-  p.ab %>%
  dplyr::select(dpi.f, dpi, primary_treatment, tes.new, elisa_od, quantity1, sex, band_number, inf_pri)

#> I calculate variability the same way across analyses in this experiment. 
#> I use functions calculate_cv and calculate_pv found in "Variability_Functions.R". 
#> I then bootstrap the original data by resampling the data with replacement (each bootstrap sample is created by
#> randomly re-sampling all points in a treatment group - some points may appear multiple times in the same bootstrap
#> sample, while others may be excluded). The 95% confidence intervals are then calculated from the bootstrap distribution
#> corresponding to the 2.5th and 97.5th percentiles.

#Read in variability functions
source("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/R_Scripts/Variability_Functions.R")

#Calculate Variability in Antibody Levels Primary
# Set the number of bootstrap replicates
n_boot <- 1000

# Pipeline to calculate CV and PV, and their 95% confidence intervals of antibody levels
# Using p.aba (**ALL individuals regardless of infection status**)
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
    
    # Bootstrap for each metric
    cv_bootstrap = list(replicate(n_boot, calculate_cv(sample(elisa_od, replace = TRUE)))), #1000 new CV calculations from sampling the data 1000 times randomly
    pv_bootstrap = list(replicate(n_boot, calculate_pv(sample(elisa_od, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for CV and PV
  mutate(
    cv_lower_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.025)), #lower confidence interval (2.5th percentile)
    cv_upper_ci = map_dbl(cv_bootstrap, ~ quantile(.x, 0.975)), #upper confidence interval (97.5th percentile)
    
    pv_lower_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.025)),
    pv_upper_ci = map_dbl(pv_bootstrap, ~ quantile(.x, 0.975))
  ) %>%
  # Remove the bootstrap columns to clean up the dataset
  dplyr::select(-cv_bootstrap, -pv_bootstrap) %>%
  ungroup()

#add to df and format names
p.aba.m <- left_join(p.aba, a.cv, by=c("dpi.f", "primary_treatment", "band_number"))
p.aba.m$elisa_od <- p.aba.m$elisa_od.x
p.aba.m <- p.aba.m %>%
  dplyr::select(-elisa_od.y, -elisa_od.x)

#just days samples were taken
p.aba.m <- p.aba.m %>%
  filter(dpi.f %in% c(-8, 14, 41))

#summary table; each individual bird is assigned the overall variability metric. For clarity, a df showing just the variability metric for each dpi x primary_treatment combination
  #This is achieved by taking the mean of each group as each individual bird in each group will have the same quantification, so mean condenses this down to one value.
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
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  ) %>%
  as_tibble()
#write.csv(summary_tibble_priming_ab, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/Antibody_Variability_Primary.csv", row.names=FALSE)

#FIGURE 2 raw data + predictions
fig2.raw<-ggplot(data = p.abt, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +
  # Original jittered data points for elisa_od
  geom_jitter(size = 1.5, alpha = 0.5, width = 0.15, height = 0, shape=16) +
  
  # Line and points with error bars for predicted values
  geom_path(data = dat.new, aes(x = dpi.f, y = yhat, group = primary_treatment,
                                color = primary_treatment), alpha=0.75) +
  geom_errorbar(data = dat.new, aes(ymin = Lower, ymax = Upper, x = dpi.f, y = yhat),
                color = "black", width = 0.025) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat, color = primary_treatment),
             size = 2, alpha = 0.75, shape = 16) +
  geom_point(data = dat.new, aes(x = dpi.f, y = yhat),
             color = "black", size = 2, alpha = 0.5, shape = 1, stroke = 0.1) +
  
  # Labels and axis
  labs(y = "Anti-MG IgY Antibody Levels (OD)", 
       x = "Days Post Primary Inoculation", 
       shape = "Variability Metric", 
       color = "Primary Treatment", 
       fill = "Primary Treatment") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  coord_cartesian(ylim=c(0.03,0.2))+
  scale_y_continuous(breaks=c(0.025, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.2))+
  # Theme settings
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    legend.background = element_rect(linewidth = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size=12),
    legend.text = element_text(size=12)
  )

fig2.raw


#Figure 2 var only
fig2.pv<-ggplot(data = p.aba.m, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +
  # PV
  geom_point(data = summary_tibble_priming_ab, aes(x = dpi.f, y = PV, fill = primary_treatment), shape=17,
             size = 3, alpha = 1) +
  geom_errorbar(data = summary_tibble_priming_ab, aes(x = dpi.f, ymin = lower_ci_pv, ymax = upper_ci_pv,
                                                      group = primary_treatment, color = primary_treatment),
                width = 0.05, alpha = 0.75) +
  geom_path(data = summary_tibble_priming_ab, aes(x = dpi.f, y = PV, color = primary_treatment,
                                                  group = primary_treatment), lty = "solid", alpha=0.75) +

    # Labels and axis
  labs(y = "Variability in IgY Antibody Levels (PV)", 
       x = "Days Post Primary Inoculation", 
       color = "Primary Treatment", 
       fill = "Primary Treatment") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  coord_cartesian(ylim=c(0,0.31))+
  
  # Theme settings
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size=12),
    legend.text = element_text(size=12)
  )

fig2.pv

#Combine two panels of Figure 2 with tags
fig2 <- (fig2.raw + fig2.pv & 
           theme(legend.position = 'top',
                 plot.tag = element_text(size = 20, face="bold"))) + 
  patchwork::plot_annotation(tag_levels = 'a')

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/Fig2.png",
#       plot=fig2, width=10, height=6, dpi = 300)

# ggsave(filename ="/Users/jesse/Documents/Virginia Tech/Research/EEID Heterogeneity Grant/Expt1A (VA94 Dose)/Manuscript/Final Scientific Reports Submission/Review/Final Review/Final Submission Mar2026/Final Submission/Figures/Fig2.png",
#        plot=fig2, width=10, height=6, dpi = 300)

#Figure S1
figS1<-ggplot(data = p.aba.m, aes(x = dpi.f, y = elisa_od, color = primary_treatment)) +

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
       color = "Primary Treatment", 
       fill = "Primary Treatment",
       shape = "Variability Metric",
       lty="Variability Metric") +
  
  # Custom color and fill scales
  scale_color_manual(values = c(pri_colors)) +
  scale_fill_manual(values = c(pri_colors)) +
  
  scale_shape_manual(
    name = "Variability Metric",
    values = c("PV" = 19, "CV" = 19),  # Set shape type
    labels = c("PV" = "PV", "CV" = "CV"))+
  
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    lty = guide_legend(order = 2)
  )+

  # Theme settings
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.y.right = element_text(color = "black", size = 15),
    axis.title.y.right = element_text(color = "black", size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.background = element_rect(size = 0.25, linetype = "solid"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size=15),
    legend.text = element_text(size=15)
  )

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/FigS1.png",
#       plot=figS1, width=9, height=6, dpi = 300)

####Are antibody levels equally variable across dpi and primary_treatment?####
#Brown-Forsythe
leveneTest(elisa_od ~ interaction(primary_treatment, dpi.f), data = p.aba, center = median)

dat41 <- p.aba %>% filter(dpi.f == 41)

# Get all pairwise combinations of groups
pairwise_results_ab_41 <- combn(unique(dat41$primary_treatment), 2, simplify = FALSE) %>%
  map_df(~{
    subdat <- dat41 %>% filter(primary_treatment %in% .x)
    test <- leveneTest(elisa_od ~ primary_treatment, data = subdat, center = median)
    tibble(
      group1 = .x[1],
      group2 = .x[2],
      F_value = test[1, "F value"],
      df = test[1, "Df"],
      df_denom = test[2, "Df"],
      p_value = test[1, "Pr(>F)"]
    )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Adjust for multiple comparisons

# write_csv(pairwise_results_ab_41,
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/BF Results/BF_pairwise_results_ab_41.csv")

dat14 <- p.aba %>% filter(dpi.f == 14)

# Get all pairwise combinations of groups
pairwise_results_ab_14 <- combn(unique(dat14$primary_treatment), 2, simplify = FALSE) %>%
  map_df(~{
    subdat <- dat14 %>% filter(primary_treatment %in% .x)
    test <- leveneTest(elisa_od ~ primary_treatment, data = subdat, center = median)
    tibble(
      group1 = .x[1],
      group2 = .x[2],
      df = test[1, "Df"],
      df_denom = test[2, "Df"],
      F_value = test[1, "F value"],
      p_value = test[1, "Pr(>F)"]
    )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Adjust for multiple comparisons

# write_csv(pairwise_results_ab_14,
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/BF Results/BF_pairwise_results_ab_14.csv")

####Analysis #2 Does individual variation in antibody response to infection predict reinfection probability?####
####Do Individual Differences in Antibody Levels Predict Secondary Susceptibility?
#For this analysis, I am asking whether individual variation on any days that antibodies were measured predict whether that individual will become reinfected upon 
#secondary challenge. 
  #First, I look at whether individual variation in antibody levels (antibodies from now on) 8 days prior to initial infection are predictive
  #of reinfection with the prediction that they will not be. Because individuals were challenged with 5 different secondary doses[log10], secondary dose[log10] is used
  #as a fixed effect. This first baseline analysis shows the effect of secondary dose[log10] on reinfection probability, which is expected to increase as secondary
  #dose[log10] increases. It is not directly reported in the manuscript as antibody levels were not predictive of susceptibility as predicted.
  #Next, I look at whether antibodies on day 14 post-primary challenge are predictive of susceptibility.
  #Finally, I look at whether antibodies on day 41 post-primary challenge are predictive of susceptibility. 
  #I look at these three days independently with three separate models with the objective of asking whether one measure of antibody levels at any point during priming
  #infection (before, during, or after) can be used to predict reinfection probability.

#New df
s.ab <- m.ab 

s.ab %>%
  filter(dpi == -8)%>%
  dplyr::select(primary_treatment, secondary_dose, inf_sec)%>%
  tbl_summary(
    by=inf_sec
  )%>%
  modify_header(
    label ~ "**Sample Size Secondary**"
  )

#Format wide
s.ab.w <- s.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi %in% c(-8, 14, 41))

wibird <- s.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)

#DPPI -8
wibird$elisa_od_pre <- wibird$`elisa_od_-8`

wibird <- wibird %>% mutate(totalTES = sum()) #Total eye score = sum of all eye scores for each bird
wibird$sec_dose.n <- wibird$secondary_dose+1
wibird$log10.sec_dose <- round(log10(wibird$sec_dose.n), digits = 2)

wibird$log10.sec_dose <- as.numeric(wibird$log10.sec_dose)

wibird$secondary_dose_fct <- factor(wibird$secondary_dose, levels = c(0, 30, 100, 300, 7000))

sec_dose_names <- c(
  "0" = "0  CCU/mL",
  "1.49" = "30 CCU/mL",
  "2" = "100 CCU/mL",
  "2.48" = "300 CCU/mL",
  "3.85" = "7,000 CCU/mL"
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

#write.csv(final_summary_pre, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Antibodies and Susceptibility/glm_ab_pre_full_summary.csv", row.names = TRUE)

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

#Table S2 Part 1
#write.csv(final_summary_14, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Antibodies and Susceptibility/glm_ab14_full_summary.csv", row.names = TRUE)

mod <- glm.ab14
dat.new14=expand.grid(log10.sec_dose=unique(wibird14$log10.sec_dose),
                    inf_sec=unique(wibird14$inf_sec),
                    elisa_od_14 = unique(wibird14$elisa_od_14))
dat.new14$yhat=predict(mod, type="response", newdata = dat.new14)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new14, se.fit =T)
#bind se's and fitted points
dat.new14 = cbind(dat.new14, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new14 <- transform(dat.new14,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

wibird14$dpi <- "DPPI 14"
dpi_names14 <- c(
  "DPPI 14" = "DPPI 14"
)

sus_14_raw<-ggplot(wibird14, aes(x=(elisa_od_14), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_line(data=dat.new14, aes(x=(elisa_od_14), y=yhat), alpha = 1, linewidth =1)+
  geom_ribbon(data = dat.new14, aes(ymin=Lower, ymax = Upper, x= elisa_od_14, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  geom_jitter(alpha=0.5, height=0.01, width=0)+
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7,000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7,000"))+
  labs(y= "Reinfection Status", color="Secondary Dose (CCU/mL)", fill ="Secondary Dose (CCU/mL)")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust =1, size=12))+
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

#Table S2 Part 2
#write.csv(final_summary_41, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Antibodies and Susceptibility/glm_ab41_full_summary.csv", row.names = TRUE)


mod <- glm.ab41
dat.new41=expand.grid(log10.sec_dose=unique(wibird41$log10.sec_dose),
                    inf_sec=unique(wibird41$inf_sec),
                    primary_treatment=unique(wibird41$primary_treatment),
                    #elisa_od_14_new=unique(wibird41$elisa_od_14_new),
                    elisa_od_41 = unique(wibird41$elisa_od_41))
dat.new41$yhat=predict(mod, type="response", newdata = dat.new41)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new41, se.fit =T)
#bind se's and fitted points
dat.new41 = cbind(dat.new41, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new41 <- transform(dat.new41,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

wibird41$dpi <- "DPPI 41"
dpi_names41 <- c(
  "DPPI 41" = "DPPI 41"
)

sus_41_raw<-ggplot(wibird41, aes(x=(elisa_od_41), y=(inf_sec), color=as.factor(log10.sec_dose)))+
  geom_jitter(alpha=0.5, height=0.01, width=0)+
  geom_line(data=dat.new41, aes(x=(elisa_od_41), y=yhat), alpha = 1, linewidth =1)+
  geom_ribbon(data = dat.new41, aes(ymin=Lower, ymax = Upper, x= elisa_od_41, y=yhat, fill=as.factor(log10.sec_dose)), alpha = 0.05, linetype="dashed", size=0.1) +
  scale_color_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7,000"))+
  scale_fill_manual(values = c(sec_colors), labels=c("0", "30", "100", "300", "7,000"))+
  labs(x="IgY Antibody Levels [OD]", y= "Reinfection Status", color="Secondary Dose (CCU/mL)", fill ="Secondary Dose (CCU/mL)")+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01))+
  theme(axis.text.x = element_text(angle=45, hjust =1, size = 12),
        legend.position = "none")+
  facet_grid(~dpi~log10.sec_dose, labeller = as_labeller(c(dpi_names41,sec_dose_names)))
sus_41_raw

#Compare Models
aictab(cand.set=list(glm.ab.pre, glm.ab14, glm.ab41), modnames=c("glm.ab.pre","glm.ab14", "glm.ab41"))

#Figure 3
#Combine two panels of Figure 2 with tags
Fig3 <- (sus_14_raw / sus_41_raw & 
           theme(plot.tag = element_text(size = 20, face="bold"))) + 
  patchwork::plot_annotation(tag_levels = 'a')


# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/Fig3.png",
#         plot=Fig3, width=9, height=6, dpi = 300)


####Analysis #3 Variability in Antibody levels and Susceptibility from Hawley et al., 2024####
#extract variability metrics from summary_tibble_priming_ab a.cv calculations
summary_tibble_priming_ab %>% filter(dpi.f == 41)

source("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/R_Scripts/CIs_Sus_2024.R")


# Filter to just dpi.f == 41
summary_41 <- summary_tibble_priming_ab %>%
  filter(dpi.f == 41) %>%
  arrange(primary_treatment)

# Create het.df directly from summary_41
het.df <- summary_41 %>%
  transmute(
    dose = primary_treatment,
    ab_pv = PV,
    ab_cv = CV,
    ci_lower_pv = lower_ci_pv,
    ci_upper_pv = upper_ci_pv,
    ci_lower_cv = lower_ci_cv,
    ci_upper_cv = upper_ci_cv,
    dpi = dpi.f
  ) %>%
  #Susceptibility estimates (CV Hawley et al., 2024)
  mutate(
    CV_Hawley2024 = c(0.899, 1.630, 2.511),
    CV_Hawley2024_calc = CoV_calc,
    mean_CV_boot = cv_summary$mean_CV,
    ci_lower_boot = ci_lower_boot,
    ci_upper_boot = ci_upper_boot
  )



#Figure 4
Fig4 <- ggplot(het.df, aes(x = ab_pv, y = CV_Hawley2024, color=dose)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_lower_pv, xmax= ci_upper_pv, x=ab_pv), width=0)+
  geom_errorbar(aes(ymin = ci_lower_boot, ymax = ci_upper_boot, y = CV_Hawley2024), width= 0)+
  scale_color_manual(values=pri_colors)+
  coord_cartesian(ylim=c(0,4), xlim=(c(0,0.25)))+
  labs(x = "Variability in IgY Antibody Levels (PV)", 
       y = "Variability in Susceptibility (CV) Hawley et al., 2024", 
       color="Primary Treatment") 

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/Fig4.png",
#         plot=Fig4, width=10, height=8, dpi = 300)


####Analysis #4 Variability in Pathogen Load and Pathology During Secondary Challenge####
#df with only birds high secondary dose after reinfection
sec.ab <- m.ab %>%
  filter(dpi > 41)

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

#Calculate max pathogen load 7000 secondary
max.quant.sec <- sec.ab %>%
  filter(secondary_dose == "7000") %>%
  group_by(band_number, primary_treatment, secondary_dose, inf_sec) %>%
  reframe(max_quantity = max(quantity, na.rm = TRUE))

#Calculate max pathology 7000 secondary
max.tes.sec <- sec.ab %>%
  filter(secondary_dose == "7000") %>%
  group_by(band_number, primary_treatment, secondary_dose, inf_sec)%>%
  reframe(max_tes = max(tes, na.rm=TRUE))

#Join
max.sec <- max.quant.sec %>%
  left_join(max.tes.sec,  by = c("band_number","primary_treatment","secondary_dose","inf_sec")) 

max.sec$max_quantity1 <- max.sec$max_quantity+1

max.sec %>%
  dplyr::select(primary_treatment, secondary_dose, inf_sec) %>%
  group_by(secondary_dose, primary_treatment, inf_sec) %>%
  tbl_summary(
    by = primary_treatment, # Grouping by infected secondary
    statistic = list(inf_sec ~ "{n} ({p}%)") # Count occurrences and show percentages
  ) %>%
  modify_header(label ~ "**Infected Secondary**") 


max.sec <- max.sec %>%
  mutate(
    primary_treatment_names = factor(
      dplyr::recode(
        primary_treatment,
        "Sham" = "Sham Priming",
        "Low"  = "Low Priming",
        "High" = "High Priming"
      ),
      levels = c("Sham Priming", "Low Priming", "High Priming")
    )
  )

####Hurdle model: 1) Logistic infection status

#Complete separation: All sham = 0 inf_sec
with(max.sec, table(primary_treatment, inf_sec))

#Use Firth logistic regression
m_inf_firth <- glm(
  inf_sec ~ primary_treatment,
  data = max.sec,
  family = binomial("logit"),
  method = "brglmFit",
  type = "AS_mean"
)

simulateResiduals(m_inf_firth, plot=T)

tidy(m_inf_firth, exponentiate = TRUE, conf.int = TRUE)

summary(m_inf_firth)

emmeans(m_inf_firth, ~ primary_treatment, type = "response")


#Variability secondary
max.sec$max_tes <- max.sec$max_tes+0.001

#max_quantity1 is max_quantity + 1 <- add 1 because calculations do not work with 0s; added to all so does not affect variability
#log10(1) = 0 which does not work for variability calculations so add small constant
max.sec$lmax_quantity <- log10(max.sec$max_quantity+1)+0.001


#All 7000 dose birds
max.s.v <- max.sec %>% 
  group_by(primary_treatment_names) %>%
  reframe(
    # Metrics for max_tes
    max_tes = max_tes,
    mean_tes = mean(max_tes-0.001), #to get actual mean
    bird_cv_tes = calculate_cv(max_tes),
    bird_sd_tes = sd(max_tes - 0.001),
    bird_pv_tes = calculate_pv(max_tes),
    bird_se_tes = SE(max_tes-0.001), #to get actual SE
    
    # Metrics for max_quantity
    max_quantity = max_quantity1,
    mean_quantity = mean(max_quantity1)-1,
    bird_cv_quantity = calculate_cv(max_quantity1),
    bird_sd_quantity = sd(max_quantity1),
    bird_pv_quantity = calculate_pv(max_quantity1),
    bird_se_quantity = SE(max_quantity1),
    
    # Metrics for log10(max_quantity)
    lmax_quantity = lmax_quantity,
    mean_lmax_quantity = mean(lmax_quantity),
    bird_cv_lmax_quantity = calculate_cv(lmax_quantity),
    bird_sd_lmax_quantity = sd(lmax_quantity),
    bird_pv_lmax_quantity = calculate_pv(lmax_quantity),
    bird_se_lmax_quantity = SE(lmax_quantity),
    
    #Other info
    band_number = band_number,
    primary_treatment = primary_treatment,
    
    # Bootstrap for confidence intervals (max_tes)
    cv_bootstrap_tes = list(replicate(n_boot, calculate_cv(sample(max_tes, replace = TRUE)))),
    pv_bootstrap_tes = list(replicate(n_boot, calculate_pv(sample(max_tes, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (max_quantity)
    cv_bootstrap_quantity = list(replicate(n_boot, calculate_cv(sample(max_quantity, replace = TRUE)))),
    pv_bootstrap_quantity = list(replicate(n_boot, calculate_pv(sample(max_quantity, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (lmax_quantity)
    cv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_cv(sample(lmax_quantity, replace = TRUE)))),
    pv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_pv(sample(lmax_quantity, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for both max_tes and max_quantity
  mutate(
    # max_tes CIs
    cv_lower_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.025)),
    cv_upper_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.975)),
    pv_lower_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.025)),
    pv_upper_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.975)),
    
    # max_quantity CIs
    cv_lower_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    
    # lmax_quantity CIs
    cv_lower_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  # select(-cv_bootstrap_tes, -v2_bootstrap_tes, -pv_bootstrap_tes,
  #        -cv_bootstrap_quantity, -v2_bootstrap_quantity, -pv_bootstrap_quantity) %>%
  ungroup()

# Add to df and format names
max.sec.v <- left_join(max.sec, max.s.v, by = c("primary_treatment_names", "band_number"))

# Summary table for both max_tes and max_quantity
summary_tibble.s.v <- max.s.v %>%
  group_by(primary_treatment_names) %>%
  summarize(
    # Summary statistics for max_tes
    CV_tes = mean(bird_cv_tes, na.rm = TRUE),
    SD_tes = mean(bird_sd_tes, na.rm = TRUE),
    PV_tes = mean(bird_pv_tes, na.rm = TRUE),
    SE_tes = mean(bird_se_tes),
    
    lower_ci_pv_tes = mean(pv_lower_ci_tes),
    upper_ci_pv_tes = mean(pv_upper_ci_tes),
    lower_ci_cv_tes = mean(cv_lower_ci_tes),
    upper_ci_cv_tes = mean(cv_upper_ci_tes),
    
    # Summary statistics for max_quantity
    CV_quantity = mean(bird_cv_quantity, na.rm = TRUE),
    SD_quantity = mean(bird_sd_quantity, na.rm = TRUE),
    PV_quantity = mean(bird_pv_quantity, na.rm = TRUE),
    SE_quantity = mean(bird_se_quantity),
    
    lower_ci_pv_quantity = mean(pv_lower_ci_quantity),
    upper_ci_pv_quantity = mean(pv_upper_ci_quantity),
    lower_ci_cv_quantity = mean(cv_lower_ci_quantity),
    upper_ci_cv_quantity = mean(cv_upper_ci_quantity),
    
    # Summary statistics for lmax_quantity
    CV_lmax_quantity = mean(bird_cv_lmax_quantity, na.rm = TRUE),
    SD_lmax_quantity = mean(bird_sd_lmax_quantity, na.rm = TRUE),
    PV_lmax_quantity = mean(bird_pv_lmax_quantity, na.rm = TRUE),
    SE_lmax_quantity = mean(bird_se_lmax_quantity),
    
    lower_ci_pv_lmax_quantity = mean(pv_lower_ci_lmax_quantity),
    upper_ci_pv_lmax_quantity = mean(pv_upper_ci_lmax_quantity),
    lower_ci_cv_lmax_quantity = mean(cv_lower_ci_lmax_quantity),
    upper_ci_cv_lmax_quantity = mean(cv_upper_ci_lmax_quantity),
    
    # Number of individuals
    n_individuals = n_distinct(band_number),
    mean_tes = mean(mean_tes),
    mean_quantity = mean(mean_quantity),
    mean_lmax_quantity = mean(mean_lmax_quantity)-0.001
  ) %>%
  as_tibble()

#Variability metrics
#Table S3
s.var <- summary_tibble.s.v %>%
  dplyr::select(primary_treatment_names, CV_tes, PV_tes, mean_tes, SE_tes,
                CV_quantity, PV_quantity, mean_quantity, SE_quantity,
                CV_lmax_quantity, PV_lmax_quantity, mean_lmax_quantity, SE_lmax_quantity,
                n_individuals)

#write.csv(s.var, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Supplementary/TableS3.csv", row.names = TRUE)

#Visualize difference in variability metrics Max Eye Score Secondary
max_sec_tes.v<- ggplot(summary_tibble.s.v) +
  geom_point(aes(x = primary_treatment_names, y = PV_tes, color = primary_treatment_names), shape = 17, size=3) +
  geom_errorbar(aes(x = primary_treatment_names, ymin = lower_ci_pv_tes, ymax = upper_ci_pv_tes, y = PV_tes, color = primary_treatment_names), 
                width = 0, alpha=0.75) +
  scale_color_manual(values = pri_colors) +
  labs(
    y = "Variability [PV] \n Maximum Eyescore",
    x = NULL,
    color = "Primary Treatment") +
  scale_y_continuous(limits = c(0,1))+
  theme(legend.position = "none")
max_sec_tes.v


#Visualize difference in variability metrics Max Pathogen Load Secondary
max_sec_path.v <- ggplot(summary_tibble.s.v) +
  geom_point(aes(x = primary_treatment_names, y = PV_lmax_quantity, color = primary_treatment_names), shape = 17, size=3) +
  geom_errorbar(aes(x = primary_treatment_names, ymin = lower_ci_pv_lmax_quantity, ymax = upper_ci_pv_lmax_quantity, y = PV_lmax_quantity, color = primary_treatment_names), 
                width = 0., alpha = 0.75) +
  
  scale_color_manual(values = pri_colors) +
  labs(
    y = "Variability [PV] \n Log10(Maximum Pathogen Load)",
    x = NULL,
    color = "Primary Treatment") +
  scale_y_continuous(limits = c(0,1))
max_sec_path.v

#Visualize max TES Secondary
tes.max.s <- ggplot(max.sec, aes(x=primary_treatment_names, y=max_tes, color=primary_treatment_names))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  #geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5, binwidth = 0.1, stackratio = 1, alpha=0.75, aes(fill=primary_treatment_names))+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_y_continuous(limits=c(0,6), breaks=c(0, 2, 4, 6))+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  labs(x=NULL, y="Maximum Eyescore", color="Primary Treatment", fill= "Primary Treatment")+
  theme(legend.position = "none")
tes.max.s  


#Max log10(quantity) secondary
quant.max.s.log10 <- ggplot(max.sec, aes(x=primary_treatment_names, y=max_quantity1, color=primary_treatment_names))+
  geom_hline(yintercept = 50, color="black", alpha=0.5, lty="dashed")+
  geom_boxplot(outlier.shape = 8, width=0.5)+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  scale_y_log10(limits = c(1, 1e6))+
  labs(x=NULL, y="Log10(Maximum Pathogen Load)", color="Primary Treatment", fill= "Primary Treatment")
quant.max.s.log10


####Variability in Secondary Responses Reinfected Birds Only####
#All 7000 dose birds; reinfected birds only
max.s.v.i <- max.sec %>% 
  filter(inf_sec == 1)%>% #only birds that become reinfected
  group_by(primary_treatment_names) %>%
  reframe(
    # Metrics for max_tes
    max_tes = max_tes,
    mean_tes = mean(max_tes-0.001), #to get actual mean
    bird_cv_tes = calculate_cv(max_tes),
    bird_sd_tes = sd(max_tes - 0.001),
    bird_pv_tes = calculate_pv(max_tes),
    bird_se_tes = SE(max_tes-0.001), #to get actual SE
    
    # Metrics for max_quantity
    max_quantity = max_quantity1,
    mean_quantity = mean(max_quantity1)-1,
    bird_cv_quantity = calculate_cv(max_quantity1),
    bird_sd_quantity = sd(max_quantity1),
    bird_pv_quantity = calculate_pv(max_quantity1),
    bird_se_quantity = SE(max_quantity1),
    
    # Metrics for log10(max_quantity)
    lmax_quantity = lmax_quantity,
    mean_lmax_quantity = mean(lmax_quantity),
    bird_cv_lmax_quantity = calculate_cv(lmax_quantity),
    bird_sd_lmax_quantity = sd(lmax_quantity),
    bird_pv_lmax_quantity = calculate_pv(lmax_quantity),
    bird_se_lmax_quantity = SE(lmax_quantity),
    
    #Other info
    band_number = band_number,
    primary_treatment_names = primary_treatment_names,
    
    # Bootstrap for confidence intervals (max_tes)
    cv_bootstrap_tes = list(replicate(n_boot, calculate_cv(sample(max_tes, replace = TRUE)))),
    pv_bootstrap_tes = list(replicate(n_boot, calculate_pv(sample(max_tes, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (max_quantity)
    cv_bootstrap_quantity = list(replicate(n_boot, calculate_cv(sample(max_quantity, replace = TRUE)))),
    pv_bootstrap_quantity = list(replicate(n_boot, calculate_pv(sample(max_quantity, replace = TRUE)))),
    
    # Bootstrap for confidence intervals (lmax_quantity)
    cv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_cv(sample(lmax_quantity, replace = TRUE)))),
    pv_bootstrap_lmax_quantity = list(replicate(n_boot, calculate_pv(sample(lmax_quantity, replace = TRUE))))
  ) %>%
  # Calculate 95% confidence intervals for both max_tes and max_quantity
  mutate(
    # max_tes CIs
    cv_lower_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.025)),
    cv_upper_ci_tes = map_dbl(cv_bootstrap_tes, ~ quantile(.x, 0.975)),
    pv_lower_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.025)),
    pv_upper_ci_tes = map_dbl(pv_bootstrap_tes, ~ quantile(.x, 0.975)),
    
    # max_quantity CIs
    cv_lower_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.975)),
    
    # lmax_quantity CIs
    cv_lower_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    cv_upper_ci_lmax_quantity = map_dbl(cv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975)),
    pv_lower_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.025)),
    pv_upper_ci_lmax_quantity = map_dbl(pv_bootstrap_lmax_quantity, ~ quantile(.x, 0.975))
  ) %>%
  # Remove bootstrap columns
  # select(-cv_bootstrap_tes, -v2_bootstrap_tes, -pv_bootstrap_tes,
  #        -cv_bootstrap_quantity, -v2_bootstrap_quantity, -pv_bootstrap_quantity) %>%
  ungroup()

# Add to df and format names
max.sec.v.i <- left_join(max.sec, max.s.v.i, by = c("primary_treatment_names", "band_number"))

# Summary table for both max_tes and max_quantity
summary_tibble.s.v.i <- max.s.v.i %>%
  group_by(primary_treatment_names) %>%
  summarize(
    # Summary statistics for max_tes
    CV_tes = mean(bird_cv_tes, na.rm = TRUE),
    SD_tes = mean(bird_sd_tes, na.rm = TRUE),
    PV_tes = mean(bird_pv_tes, na.rm = TRUE),
    SE_tes = mean(bird_se_tes),
    
    lower_ci_pv_tes = mean(pv_lower_ci_tes),
    upper_ci_pv_tes = mean(pv_upper_ci_tes),
    lower_ci_cv_tes = mean(cv_lower_ci_tes),
    upper_ci_cv_tes = mean(cv_upper_ci_tes),
    
    # Summary statistics for max_quantity
    CV_quantity = mean(bird_cv_quantity, na.rm = TRUE),
    SD_quantity = mean(bird_sd_quantity, na.rm = TRUE),
    PV_quantity = mean(bird_pv_quantity, na.rm = TRUE),
    SE_quantity = mean(bird_se_quantity),
    
    lower_ci_pv_quantity = mean(pv_lower_ci_quantity),
    upper_ci_pv_quantity = mean(pv_upper_ci_quantity),
    lower_ci_cv_quantity = mean(cv_lower_ci_quantity),
    upper_ci_cv_quantity = mean(cv_upper_ci_quantity),
    
    # Summary statistics for lmax_quantity
    CV_lmax_quantity = mean(bird_cv_lmax_quantity, na.rm = TRUE),
    SD_lmax_quantity = mean(bird_sd_lmax_quantity, na.rm = TRUE),
    PV_lmax_quantity = mean(bird_pv_lmax_quantity, na.rm = TRUE),
    SE_lmax_quantity = mean(bird_se_lmax_quantity),
    
    lower_ci_pv_lmax_quantity = mean(pv_lower_ci_lmax_quantity),
    upper_ci_pv_lmax_quantity = mean(pv_upper_ci_lmax_quantity),
    lower_ci_cv_lmax_quantity = mean(cv_lower_ci_lmax_quantity),
    upper_ci_cv_lmax_quantity = mean(cv_upper_ci_lmax_quantity),
    
    # Number of individuals
    n_individuals = n_distinct(band_number),
    mean_tes = mean(mean_tes),
    mean_quantity = mean(mean_quantity),
    mean_lmax_quantity = mean(mean_lmax_quantity)-0.001
  ) %>%
  as_tibble()

#Variability metrics
#Table S3 (infected only)
s.var.i <- summary_tibble.s.v.i %>%
  dplyr::select(primary_treatment_names, CV_tes, PV_tes, mean_tes, SE_tes, 
                CV_quantity, PV_quantity, mean_quantity, SE_quantity, 
                CV_lmax_quantity, PV_lmax_quantity, mean_lmax_quantity, SE_lmax_quantity, 
                n_individuals)

#write.csv(s.var.i, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Supplementary/TableS3.i.csv", row.names = TRUE)

#Remove high priming data points and error bars as variability from 2 points is not statistically valid
var_plot_df <- summary_tibble.s.v.i %>%
  filter(primary_treatment_names != "High Priming")

#Visualize difference in variability metrics Max Eye Score Secondary Infected Only
max_sec_tes.v.i<- ggplot() +
  
  geom_point(
    data = var_plot_df,
    aes(x = primary_treatment_names, y = PV_tes, color = primary_treatment_names),
    shape = 17, size = 3
  ) +
  geom_errorbar(
    data = var_plot_df,
    aes(x = primary_treatment_names, ymin = lower_ci_pv_tes, ymax = upper_ci_pv_tes,
        color = primary_treatment_names),
    width = 0, alpha = 0.75
  ) +
  scale_color_manual(values = pri_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    y = "Variability [PV] \n Maximum Eyescore",
    x = NULL,
    color = "Primary Treatment"
  )

max_sec_tes.v.i

#Visualize difference in variability metrics Max Pathogen Load Secondary
max_sec_path.v.i <- ggplot() +
  
  geom_point(
    data = var_plot_df,
    aes(x = primary_treatment_names, y = PV_lmax_quantity, color = primary_treatment_names),
    shape = 17, size = 3
  ) +
  geom_errorbar(
    data = var_plot_df,
    aes(x = primary_treatment_names, ymin = lower_ci_pv_lmax_quantity, ymax = upper_ci_pv_lmax_quantity,
        color = primary_treatment_names),
    width = 0, alpha = 0.75
  ) +
  scale_color_manual(values = pri_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    y = "Variability [PV] \n Log10(Maximum Pathogen Load)",
    x = NULL,
    color = "Primary Treatment"
  )

max_sec_path.v.i


#Visualize max TES Secondary
tes.max.s.i <- ggplot(max.sec %>% filter(inf_sec == 1), aes(x=primary_treatment_names, y=max_tes, color=primary_treatment_names))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_y_continuous(limits=c(0,6), breaks=c(0, 2, 4, 6))+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  labs(x=NULL, y="Maximum Eyescore", color="Primary Treatment", fill= "Primary Treatment")+
  theme(legend.position = "none")
tes.max.s.i

#Max log10(quantity) secondary
quant.max.s.log10.i <- ggplot(max.sec %>% filter(inf_sec == 1), aes(x=primary_treatment_names, y=max_quantity1, color=primary_treatment_names))+
  geom_hline(yintercept = 50, color="black", alpha=0.5, lty="dashed")+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.25, height=0, alpha=0.5)+
  scale_color_manual(values=pri_colors)+
  scale_fill_manual(values=pri_colors)+
  scale_y_log10(limits = c(1, 1e6))+
  labs(x=NULL, y="Log10(Maximum Pathogen Load)", color="Primary Treatment", fill= "Primary Treatment")
quant.max.s.log10.i


row_title <- function(label, size = 14) {
  wrap_elements(
    full = textGrob(label, gp = gpar(fontface = "bold", fontsize = size))
  )
}

row_title <- function(label, size = 16,
                      fill = "grey90", col = "black", lwd = 0.75,
                      height = 1) {
  
  wrap_elements(
    full = grobTree(
      rectGrob(
        height = unit(height, "npc"),
        gp = gpar(fill = fill, col = col, lwd = lwd)
      ),
      textGrob(
        label,
        gp = gpar(fontface = "bold", fontsize = size)
      )
    )
  )
}

tag_plot <- function(p, tag) {
  p +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(size = 20, face = "bold"),
      plot.tag.position = c(0, 0.98)  # top-left inside panel
    )
}


fig5.path <- (
  row_title("All Birds") /
    (tag_plot(quant.max.s.log10, "a") + tag_plot(max_sec_path.v, "b")) /
    row_title("Successfully Infected Birds Only") /
    (tag_plot(quant.max.s.log10.i, "c") + tag_plot(max_sec_path.v.i, "d"))
) +
  plot_layout(heights = c(0.06, .9, 0.06, .9)) &
  theme(legend.position = "none")

fig5.path

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/Fig5.png",
#         plot=fig5.path, width=12, height=10, dpi = 300)


fig6.tes <- (
  row_title("All Birds") /
    (tag_plot(tes.max.s, "a") + tag_plot(max_sec_tes.v, "b")) /
    row_title("Successfully Infected Birds Only") /
    (tag_plot(tes.max.s.i, "c") + tag_plot(max_sec_tes.v.i, "d"))
) +
  plot_layout(heights = c(0.06, .9, 0.06, .9)) &
  theme(legend.position = "none")

fig6.tes

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/Fig6.png",
#         plot=fig6.tes, width=12, height=10, dpi = 300)

#Compare Log10 Quantity vs Raw Quantity Variability
#Table S5
tableS5 <- bind_rows(
  summary_tibble.s.v %>%
    transmute(
      primary_treatment_names,
      CV = CV_quantity,
      PV = PV_quantity,
      CV_lower = lower_ci_cv_quantity,
      CV_upper = upper_ci_cv_quantity,
      PV_lower = lower_ci_pv_quantity,
      PV_upper = upper_ci_pv_quantity,
      type = "Raw Max Pathogen Load"
    ),
  summary_tibble.s.v %>%
    transmute(
      primary_treatment_names,
      CV = CV_lmax_quantity,
      PV = PV_lmax_quantity,
      CV_lower = lower_ci_cv_lmax_quantity,
      CV_upper = upper_ci_cv_lmax_quantity,
      PV_lower = lower_ci_pv_lmax_quantity,
      PV_upper = upper_ci_pv_lmax_quantity,
      type = "log10(Max Pathogen Load)"
    )
) %>%
  pivot_longer(cols = c(CV, PV), names_to = "metric", values_to = "value") %>%
  mutate(
    lower = if_else(metric == "CV", CV_lower, PV_lower),
    upper = if_else(metric == "CV", CV_upper, PV_upper),
    type = factor(type, levels = c("Raw Max Pathogen Load", "log10(Max Pathogen Load)"))
  )

#write.csv(tableS5, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Supplementary/TableS5.csv", row.names = TRUE)

# Figure S2
figS2 <- ggplot(tableS5, aes(x = primary_treatment_names, y = value, color = primary_treatment_names, shape = metric)) +
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, lty= metric), position = position_dodge(width = 0.5), width = 0) +
  facet_wrap(~type) +
  labs(
    #title = "Comparison of CV and PV with 95% Confidence Intervals",
    x = "Primary Treatment",
    y = "Variability Metric",
    color = "Primary Treatment",
    shape = "Variability Metric",
    lty = "Variability Metric"
  ) +
  scale_y_continuous(limits = c(0, 3)) +
  scale_color_manual(values = pri_colors)+
  scale_shape_manual(values = c(1, 16))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_bw(base_size = 14)

# ggsave(filename ="/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Figures/FigS2.png",
#         plot=figS2, width=9, height=6, dpi = 300)


#Infected birds only
tableS5.i <- bind_rows(
  summary_tibble.s.v.i %>%
    transmute(
      primary_treatment_names,
      CV = CV_quantity,
      PV = PV_quantity,
      CV_lower = lower_ci_cv_quantity,
      CV_upper = upper_ci_cv_quantity,
      PV_lower = lower_ci_pv_quantity,
      PV_upper = upper_ci_pv_quantity,
      type = "Raw Max Pathogen Load"
    ),
  summary_tibble.s.v.i %>%
    transmute(
      primary_treatment_names,
      CV = CV_lmax_quantity,
      PV = PV_lmax_quantity,
      CV_lower = lower_ci_cv_lmax_quantity,
      CV_upper = upper_ci_cv_lmax_quantity,
      PV_lower = lower_ci_pv_lmax_quantity,
      PV_upper = upper_ci_pv_lmax_quantity,
      type = "log10(Max Pathogen Load)"
    )
) %>%
  pivot_longer(cols = c(CV, PV), names_to = "metric", values_to = "value") %>%
  mutate(
    lower = if_else(metric == "CV", CV_lower, PV_lower),
    upper = if_else(metric == "CV", CV_upper, PV_upper),
    type = factor(type, levels = c("Raw Max Pathogen Load", "log10(Max Pathogen Load)"))
  )

#write.csv(tableS5.i, "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Supplementary//TableS5.i.csv", row.names = TRUE)

tableS5.i.plot <- tableS5.i %>%
  filter(primary_treatment_names != "High Priming")

# Figure S2 infected birds only
figS2i <- ggplot(tableS5.i.plot, aes(x = primary_treatment_names, y = value, color = primary_treatment_names, shape = metric)) +
  geom_hline(yintercept = 1, alpha=0.5)+
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, lty= metric), position = position_dodge(width = 0.5), width = 0) +
  facet_wrap(~type) +
  scale_x_discrete(drop = FALSE) +
  labs(
    #title = "Comparison of CV and PV with 95% Confidence Intervals",
    x = "Primary Treatment",
    y = "Variability Metric",
    color = "Primary Treatment",
    shape = "Variability Metric",
    lty = "Variability Metric"
  ) +
  scale_y_continuous(limits = c(0, 3)) +
  scale_color_manual(values = pri_colors)+
  scale_shape_manual(values = c(1, 16))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_bw(base_size = 14)

# ggsave(filename ="/Users/jesse/Documents/Virginia Tech/Research/EEID Heterogeneity Grant/Expt1A (VA94 Dose)/Manuscript/Final Scientific Reports Submission/Review/Figures/FigS2.png",
#         plot=figS2, width=9, height=6, dpi = 300)

#Are max eye scores equally variable across priming doses upon secondary challenge?
#Brown-Forsythe
leveneTest(max_tes ~ primary_treatment, data = max.tes.sec %>% filter(secondary_dose == "7000"), center = median)

dat <- max.tes.sec %>% filter(secondary_dose == "7000")

# Get all pairwise combinations of max tes
pairwise_results_max_tes <- combn(unique(dat$primary_treatment), 2, simplify = FALSE) %>%
  map_df(~{
    subdat <- dat %>% filter(primary_treatment %in% .x)
    test <- leveneTest(max_tes ~ primary_treatment, data = subdat, center = median)
    tibble(
      group1 = .x[1],
      group2 = .x[2],
      df = test[1, "Df"],
      df_denom = test[2, "Df"],
      F_value = test[1, "F value"],
      p_value = test[1, "Pr(>F)"]
    )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Adjust for multiple comparisons

# write_csv(pairwise_results_max_tes,
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/BF Results/BF_pairwise_results_max_tes.csv")

#Are max eye scores equally variable across priming doses upon secondary challenge?
#Brown-Forsythe Infected birds
leveneTest(max_tes ~ primary_treatment, data = max.tes.sec %>%
             filter(secondary_dose == "7000" & inf_sec == 1 & primary_treatment != "High"), center = median)

#Are max log10 pathogen loads equally variable across priming doses upon secondary challenge?
max.quant.sec$lmax_quantity <- log10(max.quant.sec$max_quantity+1)

dat <- max.quant.sec %>% filter(secondary_dose == "7000")
#Brown-Forsythe
leveneTest(lmax_quantity ~ primary_treatment, data = max.quant.sec %>% filter(secondary_dose == "7000"), center = median)

# Get all pairwise combinations of max log10 pathogen loads
pairwise_results_lmax_quantity <- combn(unique(dat$primary_treatment), 2, simplify = FALSE) %>%
  map_df(~{
    subdat <- dat %>% filter(primary_treatment %in% .x)
    test <- leveneTest(lmax_quantity ~ primary_treatment, data = subdat, center = median)
    tibble(
      group1 = .x[1],
      group2 = .x[2],
      df = test[1, "Df"],
      df_denom = test[2, "Df"],
      F_value = test[1, "F value"],
      p_value = test[1, "Pr(>F)"]
    )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Adjust for multiple comparisons

# write_csv(pairwise_results_lmax_quantity,
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/BF Results/BF_pairwise_results_lmax_quant.csv")

#Brown-Forsythe Infected only
leveneTest(lmax_quantity ~ primary_treatment, data = max.quant.sec %>%
             filter(secondary_dose == "7000" & inf_sec == 1 & primary_treatment != "High"), center = median)

# Get all pairwise combinations of max log10 pathogen loads
pairwise_results_lmax_quantity.i <- combn(unique(dat$primary_treatment), 2, simplify = FALSE) %>%
  map_df(~{
    subdat <- dat %>% filter(primary_treatment %in% .x & inf_sec == 1)
    test <- leveneTest(lmax_quantity ~ primary_treatment, data = subdat, center = median)
    tibble(
      group1 = .x[1],
      group2 = .x[2],
      F_value = test[1, "F value"],
      p_value = test[1, "Pr(>F)"]
    )
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Adjust for multiple comparisons


####Does Primary Treatment affect Maximum Pathology or log10 Pathogen load upon secondary challenge?####

##Pathology
#All birds
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

#Infected birds only
max.sec.i <- max.sec %>%
  filter(inf_sec==1)

kruskal_tes.i<- kruskal.test(max_tes ~ primary_treatment, data=max.sec.i)
kruskal_tes.i
#X2 = 6.26, df = 2, p = 0.044


dunn_tes.i <- dunn.test(max.sec.i$max_tes, max.sec.i$primary_treatment, method = "bonferroni")
dunn_tes.i

dunn_tes_df.i <- tibble(
  Comparison = dunn_tes.i$comparisons,
  Z_value = dunn_tes.i$Z,
  P_unadjusted = dunn_tes.i$P,
  P_adjusted = dunn_tes.i$P.adjusted
)

##Pathogen Load
#All birds
#Kruskal Wallis Tests
kruskal_quant<- kruskal.test(lmax_quantity ~ primary_treatment_names, data=max.sec)
kruskal_quant
#X2 = 13.604, df = 2, p = 0.00111

dunn_quant <- dunn.test(max.sec$lmax_quantity, max.sec$primary_treatment_names, method = "bonferroni")
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

#Infected Only
#Kruskal Wallis Tests
kruskal_quant.i<- kruskal.test(lmax_quantity ~ primary_treatment_names, data=max.sec.i)
kruskal_quant.i
#X2 = 4.59, df = 2, p = 0.10

dunn_quant.i <- dunn.test(max.sec.i$lmax_quantity, max.sec.i$primary_treatment_names, method = "bonferroni")
dunn_quant.i

dunn_quant_df <- tibble(
  Comparison = dunn_quant.i$comparisons,
  Z_value = dunn_quant.i$Z,
  P_unadjusted = dunn_quant.i$P,
  P_adjusted = dunn_quant.i$P.adjusted
)

# Export statistics
# write_xlsx(
#   list(KruskalWallis_TES = kruskal_tes_summary, DunnTest_TES = dunn_tes_df,
#      KruskalWallis_lQuant = kruskal_quant_summary, DunnTest_lQuant = dunn_quant_df),
#  "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Reinfection_Augments_Heterogeneity/Public/Results/Mean Pathogen Load and Pathology/combined_tes_lquant_sec_kruskal_dunn.xlsx"
# )
