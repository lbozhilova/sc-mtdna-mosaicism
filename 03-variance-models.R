####################################
#### Studying variance effects #####
####################################

# Last edited: 06/04/23 by LVB

# Description: Model the normalised heteroplasmy variance in two ways:
# - firstly, using a linear model to determine whether mutation, tissue, time  
#   are explanatory variables, and
# - secondly, simulating the Johnston & Jones (2016) model.

#----- Packages
library("tidyverse")
library("MASS")

#----- Load data
load("data/parsed/02-kimura-fits.RData")

#----- Linear model
df <- sc_kimura_fits %>% 
  filter(day >= "P6") %>% 
  mutate(daynum = as.numeric(str_remove(day, "P")))
  
model_null <- lm(nvar ~ 1, data = df)
model_full <- lm(nvar ~  mutation * daynum * tissue, data = df)
model_star <- lm(nvar ~  mutation + daynum + tissue, data = df)

summary(model_full)
# Call:
#   lm(formula = nvar ~ mutation * daynum * tissue, data = df)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.166977 -0.037684 -0.005153  0.040323  0.139544 
#
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         1.182e-01  1.915e-02   6.175 5.07e-08 ***
# mutationm.5019                     -5.105e-02  2.708e-02  -1.885   0.0640 .  
# daynum                              7.115e-04  8.762e-05   8.120 1.98e-11 ***
# tissuespleen                       -6.875e-03  2.708e-02  -0.254   0.8004    
# mutationm.5019:daynum               2.894e-04  1.239e-04   2.335   0.0227 *  
# mutationm.5019:tissuespleen         1.356e-02  3.830e-02   0.354   0.7244    
# daynum:tissuespleen                -1.789e-04  1.239e-04  -1.444   0.1537    
# mutationm.5019:daynum:tissuespleen  1.726e-04  1.753e-04   0.985   0.3283    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05651 on 64 degrees of freedom
# Multiple R-squared:  0.8515,	Adjusted R-squared:  0.8352 
# F-statistic: 52.42 on 7 and 64 DF,  p-value: < 2.2e-16

# Perform stepwise AIC model selection starting with model_star, and bounded by
# model_null and model_full.
model_best <- stepAIC(model_star, 
                      scope = list(lower = model_null, upper = model_full))
summary(model_best)
# Call:
#   lm(formula = nvar ~ mutation + daynum + mutation:daynum, data = df)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.169179 -0.034338 -0.006481  0.035006  0.137342 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            1.148e-01  1.370e-02   8.380 4.47e-12 ***
# mutationm.5019        -4.427e-02  1.938e-02  -2.285   0.0255 *  
# daynum                 6.220e-04  6.269e-05   9.921 7.53e-15 ***
# mutationm.5019:daynum  3.757e-04  8.866e-05   4.237 6.97e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05718 on 68 degrees of freedom
# Multiple R-squared:  0.8384,	Adjusted R-squared:  0.8313 
# F-statistic: 117.6 on 3 and 68 DF,  p-value: < 2.2e-16

#----- Johnston and Jones (2016) model
# Model parameters:
#   - n: avergae copy number immediately after cell divison
#           note: average copy number is therefore  n' = 3 * n / 2
#   - g: number of cell generations
#   - t: time (days)
#   - tau: rate of mtDNA turnover
#           note: half-life is tau / log(2)
seg_var <- function(n, g) { (1 - 1 / (2 * n))^g }
rep_var <- function(n, t, tau) { 4 * t / (3 * n * tau) }
h_var <- function(n, g, t, tau) { 1 - seg_var(n, g) + rep_var(n, t, tau) }

jjsim_df <- expand_grid(nc = seq(100, 2000, 100),
                        t = seq(10, 500, 10)) %>% 
  mutate(g = t / 3) %>% 
  mutate(quiescent = h_var(nc, 0, t, tau = 21 / log(2)),
         quiescent_tau = h_var(nc, 0, t, tau = 7 / log(2)),
         proliferating = h_var(nc, g, t, tau = 7 / log(2))) %>% 
  pivot_longer(4:6, names_to = "cells", values_to = "nvar")
jjsim_df$cells <- factor(jjsim_df$cells, levels = c("proliferating", "quiescent", "quiescent_tau"))

#----- Save
save(jjsim_df, model_full, model_best, file = "data/parsed/03-models.RData")
