# EDSD 2023
# Course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 3: Introduction to statistical analysis of age-period-cohort effects

rm(list=ls())
source("code/lab00_prep_session.R")

hmd <- read_rds("data_input/hmd_dts_pop_v2.rds")

# let's have a look at the mortality experience of the US females
cd <- "USA"
sx <- "male"

amin <- 0
amax <- 80
pmin <- 1950
pmax <- 2021

dt <- 
  hmd %>%
  filter(code == cd,
         sex == sx,
         age %in% amin:amax,
         year %in% pmin:pmax)

unique(dt$age)
unique(dt$year)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Descriptive APC analysis ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Lexis surface of mortality change 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# offer a great view of the dynamics over time (APC)
plot_change("USA", "male", amin, amax, pmin, pmax)

# to simplify the analysis we reduce the categories by
# grouping ages and periods in 5-year groups 
# (cohorts (P-A) will be automatically grouped too)
dt2 <- 
  dt %>% 
  rename(A = age,
         P = year) %>% 
  # creating new age categories grouped in 5 years 
  mutate(A = A - A%%5,
         P = P - P%%5,
         # creating a variable for birth cohort
         C = P - A) %>% 
  # aggregating deaths and exposures by age and period intervals
  group_by(A, P, C) %>% 
  summarise(dts = sum(dts) %>% round(),
            pop = sum(pop) %>% round()) %>% 
  ungroup() %>% 
  # estimating death rates (/100K)
  mutate(mx = 1e5*dts/pop) %>% 
  # removing weird rates
  filter(!is.na(mx) & mx!= Inf & mx!= 0)

dt2

# 4 classical plots for descriptive APC 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In APC literature it is often suggested to build exploratory plots of the 
# linking the three dimensions Age-Period-Cohort

# a. Rates versus age at death for different periods:
# rates in the same age-group connected.
p1 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(A, mx, col = factor(P), group = P))+
  scale_y_log10()+
  theme_bw()
p1
# b. Rates versus age at death for different birth cohorts:
# rates in the same birth-cohort connected.
p2 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(A, mx, col = factor(C), group = C))+
  scale_y_log10()+
  theme_bw()
p2

p1+p2
# c. Rates versus date of death:
# rates in the same age-group connected by period.
p3 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(P, mx, col = factor(A), group = A))+
  scale_y_log10()+
  theme_bw()
p3
# d. Rates versus date of date of birth:
# rates in the same age-group connected by cohort.
p4 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(C, mx, col = factor(A), group = A))+
  scale_y_log10()+
  theme_bw()
p4

(p1+p2)/(p3+p4)

# now, APC models...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Statistical APC analysis ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unique(dt2$P) %>% length()
# age- and period-specific death rates
dt2 %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  scale_y_log10()+
  theme_bw()


# fitting a fully linear APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# age
m_a <- glm(dts ~ A, offset = log(pop), family = poisson, data = dt2)
summary(m_a)
aic_a <- m_a$aic

dt2 %>% 
  mutate(pred_a = predict(m_a) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_a), col = "red")+
  scale_y_log10()+
  theme_bw()

# age-period
m_ap <- glm(dts ~ A + P, offset = log(pop), family = poisson, data = dt2)
summary(m_ap)
aic_ap <- m_ap$aic

dt2 %>% 
  mutate(pred_a = predict(m_a) %>% exp() * 1e5/ pop) %>% 
  mutate(pred_ap = predict(m_ap) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_a), col = "red")+
  geom_line(aes(A, pred_ap, group = P), col = "blue")+
  scale_y_log10()+
  theme_bw()

# age-period-cohort
m_ac <- glm(dts ~ A + P + C, offset = log(pop), family = poisson, data = dt2)
summary(m_ac)
aic_ac <- m_ac$aic

# the model didn't like it!!: it removes either P or C when fitting
# because of perfect linear dependence between the three variables (A = P - C)
# ... it is the identification problem!!
# 𝜽_𝟎+𝒂(𝜷_𝒂−𝜷_𝒄)+𝒑(𝜷_𝒑+𝜷_𝒄)


# Decomposing the APC model into linear and nonlinear effects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# cross-sectional age death rates 
h_apc <- glm(dts ~ factor(A) - 1  + I(P-1905) + factor(P) + factor(C), 
             offset = log(pop), family = poisson, data = dt2)
summary(h_apc)


# Age model ====
# ~~~~~~~~~~~~~~
# fitting a model with age as a categorical variable (nonlinear)
m_a <- glm(dts ~ factor(A), offset = log(pop), family = poisson, data = dt2)
summary(m_a)

# plotting predicted values
dt2 %>% 
  mutate(pred_a = predict(m_a) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_a), col = "blue")+
  scale_y_log10()+
  theme_bw()

# mortality age structure is changing over time!!
# we can assume they move proportionally for all ages over time 
# (same change in all ages) let's add the period into the model to see 

# Age-Period model ====
# ~~~~~~~~~~~~~~~~~~~~~
# identical change each period (assumed as a continuous change)
m_ap_lnr <- glm(dts ~ factor(A) + P, offset=log(pop), family = poisson, data = dt2)
m_ap_lnr

# plot with age-specific death rates changing constantly over time
p_ap_lnr <- 
  dt2 %>% 
  mutate(pred_ap_lnr = predict(m_ap_lnr) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_ap_lnr, group = P), col = "red")+
  scale_y_log10()+
  theme_bw()+
  labs(title = "Age-Period (linear)")
p_ap_lnr

# constant change over periods (linear effect)
m_ap_lnr$coefficients
(exp(-0.01303669)-1)*100

# now let's allow for different changes each period 
# (period assumed as a factor variable also)
m_ap_nlr <- glm(dts ~ factor(A) + factor(P), offset=log(pop), family = poisson, data = dt2)
m_ap_nlr

# plot with age-specific death rates changing non-constantly over time
p_ap_nlr <- 
  dt2 %>% 
  mutate(pred_ap_nlr = predict(m_ap_nlr) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_ap_nlr, group = P), col = "blue")+
  scale_y_log10()+
  theme_bw()+
  labs(title = "Age-Period (nonlinear)")
p_ap_nlr

# comparison between both approaches
p_ap_lnr+p_ap_nlr


# but we could also assume that changes are not occurring proportionally 
# over period but over birth cohorts
# let's add the cohort to see 


# Age-Cohort model ====
# ~~~~~~~~~~~~~~~~~~~~~

# identical change each cohort (assumed as a continuous change)
m_ac_lnr <- glm(dts ~ factor(A) + C, offset=log(pop), family = poisson, data = dt2)
m_ac_lnr

# different change each cohort (assumed as a factor variable)
m_ac_nlr <- glm(dts ~ factor(A) + factor(C), offset=log(pop), family = poisson, data = dt2)
m_ac_nlr

# plot with age-specific death rates changing constantly over Cohorts
p_ac_lnr <- 
  dt2 %>% 
  mutate(pred_ac_lnr = predict(m_ac_lnr) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_ac_lnr, group = C), col = "red")+
  scale_y_log10()+
  theme_bw()+
  labs(title = "Age-Cohort (linear)")

# plot with age-specific death rates changing non-constantly over Cohorts
p_ac_nlr <- 
  dt2 %>% 
  mutate(pred_ac_nlr = predict(m_ac_nlr) %>% exp() * 1e5/ pop) %>% 
  ggplot()+
  geom_point(aes(A, mx, col = P))+
  geom_line(aes(A, pred_ac_nlr, group = C), col = "blue")+
  scale_y_log10()+
  theme_bw()+
  labs(title = "Age-Cohort (nonlinear)")

# comparison between both approaches
p_ac_lnr+p_ac_nlr

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# what do we have so far?

# how similar are the linear models?
p_ac_lnr+p_ap_lnr
m_ap_lnr
m_ac_lnr
# what is the difference?
# exactly the same fitting, same deviance, same AIC...
# ... same model, different parameterization

# both are the same because the slope of change over time (P* or C* coefficient) 
# captures the **drift**, i.e., the sum of the true period and cohort slopes 
# (𝜷_𝒑+𝜷_𝒄).

# this is known as the **Age-drift model**

# the drift coefficient is in log scale; we have to transform it to obtain 
# the change of mortality over time

# how fast is changing mortality over time?
get_drift(m_ap_lnr)

# but age coefficients are different!
# because they are interpreted differently
# AdP: exp(coeffs) of age correspond to the death rates in the period of 
# reference (cross-sectional rates)
# AdC: exp(coeffs) of age correspond to the death rates in the cohort of 
# reference (longitudinal rates)

# this is called the Age-drift model


# how similar are the nonlinear models?
p_ac_nlr+p_ap_nlr
m_ap_nlr
m_ac_nlr
# the models are quite different
# these are the Age-Period and Age-Cohort models


# APC model ====
# ~~~~~~~~~~~~~~

# what happen if we just include everything?
s_apc <- glm(dts ~ factor(A) + factor(P) + factor(C), 
             offset = log(pop), family = poisson, data = dt2 )
s_apc

# it seems to be working... what is the issue with this?
# look at all the APC categories..., only the categories of references should 
# be lacking of coefficients

# why does the last cohort have a NA as coefficient?
# that is the way how R (other statistical software in general) treats 
# perfect multicollinearity: removes one variable for fitting the model
# this is equivalent of constraining both the first and last cohort categories 
# to be equal and zero

# the first dimension (period or cohort) added to the equation will absorb the 
# linear effect (drift), and the last one will be "detrended".  


# Holford approach 
# ~~~~~~~~~~~~~~~~
# Separate linear (drift) and nonlinear effects (APC), 
# both period and cohort effects will be "detrended" 
# the P and C reference are the first and last categories of each

# a few clues for manipulating the formula in the model
# the fitting will be the same but the interpretation of coeffcients will be
# more intuitive:
# adding a "-1" in the formula removes the intercept, and in this case, age 
# coefficients can be interpreted directly as rates.
# "I()" allows us to manipulate the formula, in this case to change the period 
# of reference in the liner trend (drift)
# "relevel()" allows us to change the category of reference for nonlinear effects
# it is possible to select the category by location or value
# relevel(factor(P), 10) | relevel(factor(P), "1950")

# cross-sectional age death rates 
h_apc <- glm(dts ~ factor(A) - 1  + I(P-1970) + factor(P) + factor(C), 
              offset = log(pop), family = poisson, data = dt2)
summary(h_apc)

# longitudinal age death rates 
h_acp <- glm(dts ~ factor(A) - 1  + I(C-1950) + factor(C) + factor(P), 
              offset = log(pop), family = poisson, data = dt2)
summary(h_acp)

# *drift*
extract_coeffs(h_apc)$drift
extract_coeffs(h_acp)$drift

drift_h <- (round(extract_coeffs(h_apc)$drift, 3) - 1)*100

# extract coefficients for plot
coef_h_apc <- extract_coeffs(h_apc)$coeffs
coef_h_acp <- extract_coeffs(h_acp)$coeffs

# APC and ACP models are the same (same drift, fitting, etc.) but age-specific 
# death rates are located either on the period or the cohort of reference 
bind_rows(coef_h_apc %>% 
            mutate(model = "APC"),
          coef_h_acp %>% 
            mutate(model = "ACP")) %>% 
  ggplot()+
  geom_line(aes(value, effect, group = model, col = model))+
  facet_wrap(~tdim, scales = "free")+
  scale_y_log10()+
  theme_bw()+
  labs(title = paste0("ADRs at the period/cohort of reference, detrended period and cohort effects, with drift = ", drift_h, "%"))


# comparing the fitting of the models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# how justified is to saturate the model with the three temporal variables APC?
# how much improves the model when adding each temporal dimension?
aic_a <- tibble(model = "A", aic = m_a$aic, dev = m_a$deviance)
aic_ad <- tibble(model = "Ad", aic = m_ap_lnr$aic, dev = m_ap_lnr$deviance)
aic_ap <- tibble(model = "AP", aic = m_ap_nlr$aic, dev = m_ap_nlr$deviance)
aic_ac <- tibble(model = "AC", aic = m_ac_nlr$aic, dev = m_ac_nlr$deviance)
aic_apc <- tibble(model = "APC", aic = h_apc$aic, dev = h_apc$deviance)

bind_rows(aic_a,
          aic_ad,
          aic_ap,
          aic_ac,
          aic_apc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Carstensen approach ====
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Separate linear (drift) and nonlinear effects (APC), 
# both period and cohort effects will be detrended, with sum/average = 0; 
# the reference is the linear trend;
# period and cohort curves are interpreted as the relative risks (RR)
# compared to the overall average / liner trend 

# advantages: 
# more flexibility in the modeling, allowing for 
# - semi-parametric terms (splines instead of factors)
# - different ways to extract the trend
# - better confidence intervals
# - better to interpret

# disadvantage: 
# - much more complex to model

# ... but the nice package "Epi" (made by Carstensen too) makes life much easier
# you can find the documentation here:
# https://cran.r-project.org/web/packages/Epi/Epi.pdf

library(Epi)

# In the Epi package we need to rename the variables as 
# A: age
# P: year
# D: deaths
# Y: exposures

dt_carst <- 
  dt2 %>% 
  rename(D = dts,
         Y = pop) %>% 
  select(-mx, -C)


apc_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdPC", 
          scale = 10^5)

apc_c$Drift
apc.plot(apc_c)
# I don't like it

# I made a function for plotting estimates from Carstensen models
plot_carst(apc_c)


# now longitudinal age-specific death rates 
acp_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_c)



# the function is very versatile and has many ways to parameterize
# for instance, not grouping ages and periods, 

# age in single-year of period and age
dt_carst_x1 <- 
  dt %>% 
  select(D = dts,
         Y = pop,
         A = age,
         P = year)

acp_factor <- 
  apc.fit(dt_carst_x1, 
          model = "factor", 
          npar = c(A=15, P=10, C=15),
          dr.extr = "Y", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_factor)

# we can fit splines for the nonlinear effects instead of categorical variables 
# smoothed effects, more convenient for analysis of the trend
acp_splines <- 
  apc.fit(dt_carst_x1, 
          model = "bs", 
          ref.c = 1950,
          npar = c(A=15, P=10, C=15),
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_splines)



# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assignment in class: ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Compare the age-period-cohort effects by sex in a country of choice
# what is the drift? how to interpret it?
# what is the cohort having the worst outcome of all? what is the reference 
# for this comparison?
# how is mortality changing? are secular mortality changes being driven 
# by period- or cohort-based changes?

