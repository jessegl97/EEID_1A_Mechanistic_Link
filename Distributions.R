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

m.ab <- m.ab %>%
  filter(band_number != 2505)

m.ab <- m.ab %>%
  filter(!(band_number %in% c(2274, 2514, 2469, 2520, 2494)))

#Fit gamma distribution

m.ab14 <- m.ab %>%
  filter(dpi == 14) %>%
  dplyr::select(dpi, elisa_od, band_number, primary_treatment)%>%
  na.omit()

fit <- fitdistr(m.ab14$elisa_od, "gamma")
# Extract the shape and rate parameters
shape <- fit$estimate["shape"]
rate <- fit$estimate["rate"]

# Define a function to fit the gamma distribution and extract parameters
fit_gamma <- function(data) {
  cleaned_data <- na.omit(data$elisa_od)
  fit <- fitdistr(cleaned_data, "gamma")
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  return(c(shape = shape, rate = rate))
}

results <- m.ab14 %>%
  group_by(primary_treatment) %>%
  reframe(fit_gamma(cur_data()))

ggplot(m.ab14, aes(x=elisa_od, fill=primary_treatment))+
  geom_histogram()+
  #geom_density()+
  facet_wrap(~primary_treatment)

####Chat GPT Shit
# Define a function to fit the gamma distribution and extract parameters
fit_gamma <- function(data) {
  cleaned_data <- na.omit(data$elisa_od)
  fit <- fitdistr(cleaned_data, "gamma")
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  scale <- 1 / rate
  return(list(shape = shape, rate = rate, scale = scale))
}

# Function to perform bootstrapping
bootstrap_gamma <- function(data, n_boot = 1000) {
  boot_results <- replicate(n_boot, {
    sample_data <- sample(data$elisa_od, replace = TRUE)
    fit <- fitdistr(sample_data, "gamma")
    return(fit$estimate)
  })
  ci_shape <- quantile(boot_results["shape", ], c(0.025, 0.975))
  ci_rate <- quantile(boot_results["rate", ], c(0.025, 0.975))
  ci_scale <- 1 / ci_rate
  return(list(ci_shape = ci_shape, ci_rate = ci_rate, ci_scale = ci_scale))
}

# Apply the functions to each group of primary_treatment
results <- m.ab14 %>%
  group_by(primary_treatment) %>%
  summarise(fit = list(fit_gamma(cur_data())),
            bootstrap = list(bootstrap_gamma(cur_data()))) %>%
  unnest_wider(fit) %>%
  unnest_wider(bootstrap)

# Print the results
print(results)

# Generate data for plotting
generate_gamma_data <- function(shape, scale, ci_shape, ci_scale, max_x) {
  x <- seq(0, max_x, length.out = 100)
  y <- dgamma(x, shape = shape, scale = scale)
  y_lower <- dgamma(x, shape = ci_shape[1], scale = ci_scale[1])
  y_upper <- dgamma(x, shape = ci_shape[2], scale = ci_scale[2])
  tibble(x = x, y = y, y_lower = y_lower, y_upper = y_upper)
}

max_x <- max(m.ab14$elisa_od, na.rm = TRUE)

plot_data <- results %>%
  rowwise() %>%
  mutate(data = list(generate_gamma_data(shape, scale, ci_shape, ci_scale, max_x))) %>%
  unnest(data)

mean_susceptibility <- m.ab14 %>%
  group_by(primary_treatment) %>%
  summarise(mean_susceptibility = mean(elisa_od, na.rm = TRUE))

# Combine the results with plot_data
plot_data <- results %>%
  rowwise() %>%
  mutate(data = list(generate_gamma_data(shape, scale, ci_shape, ci_scale, max_x))) %>%
  unnest(data) %>%
  left_join(mean_susceptibility, by = "primary_treatment")

# Plot the gamma distributions with confidence intervals, mean susceptibility, and facets
ggplot(plot_data, aes(x = x, y = y, color = primary_treatment, fill = primary_treatment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2) +
  geom_vline(aes(xintercept = mean_susceptibility), linetype = "dashed", size = 1) +
  facet_wrap(~ primary_treatment) + #, scales="free_y"
  geom_text(data = results, aes(x = max_x * 0.7, y = Inf, 
                                label = paste("Shape:", round(shape, 2), "\nScale:", round(scale, 2))), 
            vjust = 2, hjust = 0.5, color = "black") +
  labs(title = "Gamma Distributions with 95% Confidence Intervals and Mean IgG",
       x = "ELISA OD",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")

#Susceptibility distributions a la chatGPT
m.ab_sec <- m.ab %>%
  filter(dpi == 56) %>%  # Focus on data after dpi 42 for secondary infection
  group_by(band_number) %>%
  mutate(inf_sec = ifelse(any(coalesce(inf, 0) == 1), 1, 0)) %>%
  ungroup()

# Define a function to fit the gamma distribution and extract parameters
fit_gamma <- function(data) {
  cleaned_data <- na.omit(data$elisa_od)
  fit <- fitdistr(cleaned_data, "gamma")
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  scale <- 1 / rate
  return(list(shape = shape, rate = rate, scale = scale))
}

# Function to perform bootstrapping
bootstrap_gamma <- function(data, n_boot = 1000) {
  boot_results <- replicate(n_boot, {
    sample_data <- sample(data$elisa_od, replace = TRUE)
    fit <- fitdistr(sample_data, "gamma")
    return(fit$estimate)
  })
  ci_shape <- quantile(boot_results["shape", ], c(0.025, 0.975))
  ci_rate <- quantile(boot_results["rate", ], c(0.025, 0.975))
  ci_scale <- 1 / ci_rate
  return(list(ci_shape = ci_shape, ci_rate = ci_rate, ci_scale = ci_scale))
}

# Apply the functions to each secondary infection status group
results_sec <- m.ab_sec %>%
  group_by(inf_sec) %>%
  summarise(fit = list(fit_gamma(cur_data())),
            bootstrap = list(bootstrap_gamma(cur_data()))) %>%
  unnest_wider(fit) %>%
  unnest_wider(bootstrap)

# Print the results
print(results_sec)

# Generate data for plotting
generate_gamma_data <- function(shape, scale, ci_shape, ci_scale, max_x) {
  x <- seq(0, max_x, length.out = 100)
  y <- dgamma(x, shape = shape, scale = scale)
  y_lower <- dgamma(x, shape = ci_shape[1], scale = ci_scale[1])
  y_upper <- dgamma(x, shape = ci_shape[2], scale = ci_scale[2])
  tibble(x = x, y = y, y_lower = y_lower, y_upper = y_upper)
}

max_x <- max(m.ab_sec$elisa_od, na.rm = TRUE)

# Calculate mean susceptibility for each secondary infection status
mean_susceptibility_sec <- m.ab_sec %>%
  group_by(inf_sec) %>%
  summarise(mean_susceptibility = mean(elisa_od, na.rm = TRUE))

plot_data_sec <- results_sec %>%
  rowwise() %>%
  mutate(data = list(generate_gamma_data(shape, scale, ci_shape, ci_scale, max_x))) %>%
  unnest(data) %>%
  left_join(mean_susceptibility_sec, by = "inf_sec")

# Plot the gamma distributions with confidence intervals, mean susceptibility, and facets
ggplot(plot_data_sec, aes(x = x, y = y, color = as.factor(inf_sec), fill = as.factor(inf_sec))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper), alpha = 0.2) +
  geom_vline(aes(xintercept = mean_susceptibility, color = as.factor(inf_sec)), linetype = "dashed", size = 1) +
  facet_wrap(~ inf_sec) +
  geom_text(data = results_sec, aes(x = max_x * 0.7, y = Inf, 
                                    label = paste("Shape:", round(shape, 2), "\nScale:", round(scale, 2))), 
            vjust = 2, hjust = 0.5, color = "black") +
  labs(title = "Gamma Distributions with 95% Confidence Intervals and Mean Susceptibility (Secondary Infection)",
       x = "ELISA OD",
       y = "Density",
       color = "Secondary Infection Status",
       fill = "Secondary Infection Status") +
  theme_minimal() +
  theme(legend.position = "none")

#### Susceptibility Distributions
sus <- m.ab %>%
  filter(dpi == 56)%>%
  group_by(primary_treatment, secondary_dose)%>%
  summarise(prop_inf = (sum(inf_sec)/n_distinct(band_number)))
sus

ggplot(sus, aes(x=secondary_dose, y=prop_inf, color=fct_rev(primary_treatment)))+
         geom_point()+
  scale_x_continuous(trans = 'log10', limits = c(0.1, 10000))
#+
  facet_wrap(~primary_treatment)

sus_dist <- m.ab %>%
  dplyr::select(dpi, primary_treatment, secondary_dose, band_number, elisa_od, inf, inf_sec, quantity1)%>%
  filter(dpi == 56)%>%
  na.omit()


summary(sus_dist)
head(sus_dist)

####More ChatGPT shit: this time beta regression of susceptibility distributions
library(betareg)

dose_response_model <- glm(inf_sec ~ log10(secondary_dose) + primary_treatment, data=sus_dist, family=binomial)

summary(dose_response_model)

#predict probabilities from logistic model
#maybe for antibodies do the same thing but with antibodies instead of inf_sec?
sus_dist <- sus_dist %>%
  mutate(predicted_prob = predict(dose_response_model, type= "response"))

beta_model <- betareg(predicted_prob ~ 1, data=sus_dist)

summary(beta_model)

#extract alpha and beta parameters
alpha <- beta_model$coefficients$mean
beta <- beta_model$coefficients$precision

# Define the function to calculate the variance
beta_variance <- function(alpha, beta) {
  return((alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
}

# Calculate the variance
variance <- beta_variance(alpha, beta)
print(paste("Variance (sigma):", variance))

#plot observed data and fitted dose-response curve
ggplot(sus_dist, aes(x=log10(secondary_dose), y=inf_sec))+
  geom_jitter(height=0.01, width=0.1, alpha=0.5)+
  geom_point(data=sus, aes(x=log10(secondary_dose), y=prop_inf))+
  stat_smooth(method="glm", method.args = list(family="binomial"), se=FALSE, aes(color=fct_rev(primary_treatment)))+
  labs(title="Dose-Response Curve", x="MG Dose", y= "Infection Probability")+
  facet_wrap(~primary_treatment)

#plot predicted probabilities and beta distribution
dose_beta <-ggplot(sus_dist, aes(x= predicted_prob))+
  geom_histogram(aes(y=..density.., fill=fct_rev(primary_treatment)), bins=70, alpha=0.5) + 
  #..density.. makes the histogram a density instead by making area under histogram sum to 1
  stat_function(fun=dbeta, args = list(shape1 = -alpha, shape2 = beta), color="red", size=1)+
  labs(title = "Beta Distribution Fit", x = "Predicted Probability of Infection", y = "Density")+
  facet_wrap(~primary_treatment)

dose_beta

##ANTIBODIES
#predict probabilities from logistic model
#maybe for antibodies do the same thing but with antibodies instead of inf_sec?
m.ab.w <- m.ab %>% 
  dplyr::select(dpi, band_number, primary_treatment, secondary_dose,
                sex, elisa_od, inf_sec)%>%
  filter(dpi ==14)

wibird <- m.ab.w %>% pivot_wider(names_from = dpi,
                                 names_glue = "{.value}_{dpi}",
                                 values_from = elisa_od)



wibird <- wibird %>% mutate(totalTES = sum())
wibird<- wibird %>% na.omit(elisa_od_14)


glm.ab14<- glm(inf_sec ~ elisa_od_14 + secondary_dose, data=wibird, family=binomial())

summary(glm.ab14)

wibird <- wibird %>%
  mutate(predicted_prob = predict(glm.ab14, type= "response"))

beta_model <- betareg(predicted_prob ~ 1, data=wibird)

summary(beta_model)

#extract alpha and beta parameters
alpha <- beta_model$coefficients$mean
beta <- beta_model$coefficients$precision

# Define the function to calculate the variance
beta_variance <- function(alpha, beta) {
  return((alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)))
}

# Calculate the variance
variance <- beta_variance(alpha, beta)
print(paste("Variance (sigma):", variance))

#plot observed data and fitted dose-response curve
ggplot(sus_dist, aes(x=log10(secondary_dose), y=inf_sec))+
  geom_jitter(height=0.01, width=0.1, alpha=0.5)+
  geom_point(data=sus, aes(x=log10(secondary_dose), y=prop_inf))+
  stat_smooth(method="glm", method.args = list(family="binomial"), se=FALSE, aes(color=fct_rev(primary_treatment)))+
  labs(title="Dose-Response Curve", x="MG Dose", y= "Infection Probability")+
  facet_wrap(~primary_treatment)

#plot predicted probabilities and beta distribution from antibodies
ab_beta <- ggplot(wibird, aes(x= predicted_prob))+
  geom_histogram(aes(y=..density.., fill=fct_rev(primary_treatment)), bins=70, alpha=0.5) +
  stat_function(fun=dbeta, args = list(shape1 = alpha, shape2 = beta), color="red", size=1)+
  labs(title = "Beta Distribution Fit", x = "Predicted Probability of Infection", y = "Density")+
  facet_grid(~primary_treatment)

dose_beta
ab_beta



###### New Plan ####
##Fit dose response curves and antibody curve distributions. Compare the susceptibility distributions from these
##Heterogeneous Infection:


