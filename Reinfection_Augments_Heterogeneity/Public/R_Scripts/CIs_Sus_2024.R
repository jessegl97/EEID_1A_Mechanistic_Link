#Reinfection Augemts Heterogeneity
####95% CIs for CoV from gamma distributed susceptibility from Hawley et al., 2024####
#Requires all_pseudosets.csv, deviance_params_noday41positives.csv


ps <- read.csv("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Susceptibility Paper Code/all_pseudosets.csv")
dp <- read.csv("/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/Susceptibility Paper Code/deviance_params_noday41positives.csv")

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
    CI_lower = quantile(CV_gamma, probs = 0.025, na.rm = TRUE),
    CI_upper = quantile(CV_gamma, probs = 0.975, na.rm = TRUE)
  ) %>%
  mutate(bird.groups = factor(bird.groups, levels = c("0 va", "750 va", "30000 va"))) %>%
  arrange(bird.groups)

cv_summary

ggplot(cv_summary, aes(x=bird.groups, y=mean_CV))+
  geom_point(data=dp, aes(x=group, y=CV_gamma))+
  geom_errorbar(aes(x=bird.groups, y=mean_CV, ymin = CI_lower, ymax= CI_upper))

ci_upper_boot <- cv_summary$CI_upper
ci_lower_boot <- cv_summary$CI_lower
