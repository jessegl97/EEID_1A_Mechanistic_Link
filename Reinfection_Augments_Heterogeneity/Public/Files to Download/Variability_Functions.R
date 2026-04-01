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
#Functions for Reinfection_Augments_Heterogeneity.R
#### Define functions for calculating CV, PV, V2, and SE ####
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