#Prior pathogen exposure augments inter-individual heterogeneity in antibody levels and reinfection loads in a songbird-pathogen system

####95% CIs for CoV from gamma distributed susceptibility from Hawley et al., 2024####
  #Hawley, D. M., Pérez-Umphrey, A. A., Adelman, J. S., Fleming-Davies, A. E., Garrett-Larsen, J., Geary, S. J.,
    #Childs, L. M., & Langwig, K. E. (2024). Prior exposure to pathogens augments host heterogeneity in susceptibility
    #and has key epidemiological consequences. PLOS Pathogens, 20(9), e1012092. https://doi.org/10.1371/journal.ppat.1012092

#MIT License
# 
# Copyright (c) [2026] [Garrett-Larsen et al.]
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# REPRODUCIBILITY
#   - No setwd(); uses {here} for project-relative paths
#   - Assumes you open/run from the project root


# EXPECTED STRUCTURE
#   project_root/
#       Reinfection_Augments_Heterogeneity.R
#       Variability_Functions.R
#       CIs_Sus_2024.R
#       reinfection_response.rds
#       all_pseudosets.csv
#       deviance_params_noday41positives.csv


############################################################################

#Requires all_pseudosets.csv, deviance_params_noday41positives.csv
library(dplyr)

ps <- read.csv("all_pseudosets.csv")
dp <- read.csv("deviance_params_noday41positives.csv")

ps$shape <- ps$par1.gamma
ps$scale <- ps$par2.gamma

dp$shape <- dp$par1.gamma
dp$scale <- dp$par2.gamma

#Calculate CoV from gamma distribution 
dp <- dp %>%
  mutate(CV_gamma = 1 / sqrt(shape))

# Compute the gamma mean and standard deviation
dp <- dp %>%
  mutate(
    gamma_mean = par1.gamma * par2.gamma,  # Mean = shape * scale
    gamma_sd = sqrt(par1.gamma * par2.gamma^2),  # SD = sqrt(shape * scale^2)
    CV_gamma = gamma_sd / gamma_mean  # CV = std/mean
  )
dp
CoV_calc <- dp$CV_gamma
#Calculate CoV from pseudoreplicates 
ps <- ps %>%
  mutate(CV_gamma = 1 / sqrt(shape))

#Calculate mean CoV from pseudoreplicates
cv_summary <- ps %>%
  filter(bird.groups %in% c("0 va", "750 va", "30000 va")) %>%  # Filter for relevant groups
  group_by(bird.groups) %>%
  summarise(
    mean_CV = mean(CV_gamma, na.rm = TRUE),
    median_CV = median(CV_gamma, na.rm = TRUE),  # Also compute median to compare
    CI_lower_95 = quantile(CV_gamma, probs = 0.025, na.rm = TRUE),
    CI_upper_95 = quantile(CV_gamma, probs = 0.975, na.rm = TRUE),
    CI_lower_66 = quantile(CV_gamma, probs = 0.17, na.rm = TRUE),
    CI_upper_66 = quantile(CV_gamma, probs = 0.83, na.rm = TRUE),
  ) %>%
  mutate(bird.groups = factor(bird.groups, levels = c("0 va", "750 va", "30000 va"))) %>%
  arrange(bird.groups)

cv_summary

library(ggplot2)

ggplot(cv_summary, aes(x=bird.groups, y=mean_CV))+
  geom_point(data=dp, aes(x=group, y=CV_gamma))+
  geom_errorbar(aes(x=bird.groups, y=mean_CV, ymin = CI_lower_95, ymax= CI_upper_95))

ci_upper_boot <- cv_summary$CI_upper_95
ci_lower_boot <- cv_summary$CI_lower_95
