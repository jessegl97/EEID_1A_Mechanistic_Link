#Functions for Reinfection_Final_Code.R
#### Define functions for calculating CV, PV, and V2 ####
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

SE <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))  # Standard deviation / sqrt(n)
}