# EDSD 2023
# course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 2: Introduction to visual analysis of age-period-cohort effects

rm(list=ls())
source("code/lab00_prep_session.R")

# load data on deaths and exposures from the HMD
# (already downloaded all available information, see lab00 script)
hmd <- read_rds("data_input/hmd_dts_pop.rds")
unique(hmd$code)

# Lets select a group
# Spain, all sexes, ages 0-100, period 1910-2010 
cd <- "ESP"
sx <- "total"

amin <- 0
amax <- 100
pmin <- 1910
pmax <- 2010

dt <- 
  hmd %>%
  filter(code == cd,
         sex == sx,
         age %in% amin:amax,
         year %in% pmin:pmax) %>% 
  mutate(mx = 1e5*dts/pop,
         log_mx = log(mx))

pmort_pref <- 
  list(scale_y_log10(),
         theme_bw())

# how is mortality at age 30y in Spain, between 1910 and 2010?
ag <- 30

dt %>% 
  filter(age == ag) %>% 
  ggplot()+
  geom_point(aes(year, mx), alpha = 0.2, size = 1)+
  geom_line(aes(year, mx))+
  pmort_pref

# how was the age structure of mortality during the Spanish flu in Spain?
dt %>% 
  filter(year == 1918) %>% 
  ggplot()+
  geom_point(aes(age, mx), alpha = 0.2, size = 1)+
  geom_line(aes(age, mx))+
  pmort_pref

# What about tyhe civil war?
dt %>% 
  filter(year == 1938) %>% 
  ggplot()+
  geom_point(aes(age, mx), alpha = 0.2, size = 1)+
  geom_line(aes(age, mx))+
  pmort_pref

# It is tricky to analyze trends in stochastic processes due to the random 
# variations. To improve the analysis it is suggested to smooth the death rates.
# There are multiple methods for smoothing data. Let's see a method conceived to
# smooth mortality rates

# 1. smoothing rates ====
# ~~~~~~~~~~~~~~~~~~~~~~~
# using p-splines is a very good practice for smoothing mortality rates, and 
# it has the option of evaluating the fitting using statistical tests for 
# selecting the best fitting, e.g., BIC (1), AIC (2) (less penalization), or 
# an arbitrary penalization parameter

# Giancarlo Camarda has a great package for this: "MortalitySmooth" (2012)
# https://www.jstatsoft.org/article/view/v050i01
# however the package is not in CRAN anymore, but the good and kind Tim Riffe 
# has a copy in his GitHub. We already installed from him (see lab00 script).

# 1.1. smoothing in one dimension ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We need to enter the information in vector forms (deaths and exposures)
dts <- 
  dt %>% 
  filter(year == 1938) %>% 
  pull(dts)

pop <- 
  dt %>% 
  filter(year == 1938) %>% 
  pull(pop)

fit_1d_bic <-
  Mort1Dsmooth(x = 1:length(dts), 
               y = dts, 
               offset = log(pop), 
               # accounting for overdispersion (in Poisson context, this is the 
               # case when the variance is larger than the average)
               overdispersion = TRUE, 
               # using AIC
               method = 1)

# extracting the smoothed log rates 
fit_1d_bic$logmortality[,1]

dts2 <- 
  dt %>% 
  filter(year == 1938) %>% 
  mutate(mx_bic = 1e5*exp(fit_1d_bic$logmortality[,1]))

dts2 %>% 
  ggplot()+
  geom_point(aes(age, mx))+
  geom_line(aes(age, mx_bic))+
  pmort_pref

# could we play with the penalization of the spline?
fit_1d_aic <-
  Mort1Dsmooth(x = 1:length(dts), 
               y = dts, 
               offset = log(pop), 
               # accounting for overdispersion
               overdispersion = TRUE, 
               # using BIC, that penalizes a bit less 
               method = 2)

fit_1d_ary <-
  Mort1Dsmooth(x = 1:length(dts), 
               y = dts, 
               offset = log(pop), 
               # accounting for overdispersion
               overdispersion = TRUE, 
               # using arbitrary penalization lambda, choosing a high value
               # (more penalty)
               method = 3, lambda = 1e5)

dts3 <- 
  dts2 %>% 
  mutate(mx_sm_aic = 1e5*exp(fit_1d_aic$logmortality[,1]),
         mx_sm_ary = 1e5*exp(fit_1d_ary$logmortality[,1]))

sms <- 
  dts3 %>% 
  select(age, starts_with("mx")) %>% 
  gather(-age, key = type, value = mx)

sms %>% 
  ggplot()+
  geom_line(aes(age, mx, col = type))+
  pmort_pref

# Interpolation ====
# ~~~~~~~~~~~~~~~~~~
# Can we interpolate?
# could be useful to generate counter factual scenarios
# For interpolating with Mortality Smooth we just have to add weights, 
# 1's for data to consider, 0's where we want to interpolate

# For instance:
# In Spain mortality was very high during the civil war (1936-1939) and 
# the neighboring years
# How was mortality experienced by those aged 25?
ag <- 25

cw <- 
  dt %>% 
  filter(age == ag)

cw %>% 
  ggplot()+
  geom_point(aes(year, mx), alpha = 0.2, size = 0.5)+
  geom_line(aes(year, mx))+
  labs(title = paste0("death rates at age ", ag))+
  pmort_pref

# now, we can smooth this to see a more clear pattern
dts_cw <- 
  cw %>% 
  pull(dts)

pop_cw <- 
  cw %>% 
  pull(pop)

fit_1d <-
  Mort1Dsmooth(x = cw$year, 
               y = dts_cw, 
               offset = log(pop_cw), 
               # accounting for overdispersion
               overdispersion = TRUE, 
               # using AIC
               method = 1)

cw2 <-
  cw %>% 
  mutate(mx_sm = 1e5*exp(fit_1d$logmortality[,1]))

cw2 %>% 
  ggplot()+
  geom_point(aes(year, mx), alpha = 0.2, size = 0.5)+
  geom_line(aes(year, mx))+
  geom_line(aes(year, mx_sm), col = "red")+
  scale_y_log10()+
  labs(title = paste0("death rates at age ", ag))+
  theme_bw()


# could we construct a counterfactual scenario of mortality without the civil war? 
# removing the period 1935-1944 for interpolation

ws <- 
  cw %>% 
  select(year) %>% 
  mutate(w = ifelse(year %in% 1935:1944, 0, 1)) %>% 
  pull(w)

fit_1d_i <-
  Mort1Dsmooth(x = cw$year, 
               y = dts_cw, 
               w = ws,
               offset = log(pop_cw), 
               overdispersion = TRUE, 
               method = 1)

cw3 <- 
  cw2 %>% 
  mutate(mx_sm_i = 1e5*exp(fit_1d_i$logmortality[,1])) %>% 
  select(year, starts_with("mx")) %>% 
  gather(-year, key = type, value = mx)

cw3 %>% 
  ggplot()+
  geom_line(aes(year, mx, col = type))+
  scale_y_log10()+
  theme_bw()+
  geom_vline(xintercept = c(1935, 1944), linetype = "dashed")



# 1.2. smoothing in two dimensions ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's look at the whole data in a Lexis surface
brks <- quantile(c(min(dt$log_mx), max(dt$log_mx)), probs = seq(0, 1, 0.25))
lbls <- exp(brks) %>% round()
  
redo_lexis_shape(pmin,pmax,amin,amax)

p1 <- 
  dt %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  geom_contour(bins = 30, col="black", size=.15, alpha=0.8)+
  lexis_shape
  
p1

# we can see  a lot of random variations in death rates
# what if we need smoothness not only on one temporal dimension but in two of them?
# We can do it with the Mortality Smooth package too!!

# We need to give the information in matrix forms (age x period)
ylist <- unique(dt$year) %>% sort() # all periods in data, ordered
alist <- unique(dt$age) %>% sort() # all ages in data, ordered

# mortality matrix
deaths <- 
  matrix(dt$dts, 
         nrow = length(alist), 
         ncol = length(ylist), 
         byrow = F)
colnames(deaths) <- ylist
rownames(deaths) <- alist

# lets have a look of the deaths matrix 
deaths

# population at risk (population/exposure) matrix
exposure <- 
  matrix(dt$pop, 
         nrow = length(alist),
         ncol = length(ylist), 
         byrow = F)
colnames(exposure) <- ylist
rownames(exposure) <- alist


# smoothing mortality in 2D using p-splines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_2d <- 
  Mort2Dsmooth(x = alist, 
               y = ylist, 
               Z = deaths, 
               offset = log(exposure),
               overdispersion = TRUE, 
               method = 2)

# extracting the smoothed log rates 
fit_2d$logmortality
# transforming them to smoothed deaths
mx_smooth <- 
  (exp(fit_2d$logmortality) * 1e5)

# from matrix to tidy (long) form
smt <- 
  mx_smooth %>% 
  # keeping the row names as age values
  as_tibble(rownames = "age") %>%
  mutate(age = age %>% as.double()) %>% 
  # reshaping to long form (tidy)
  gather(-age, key = year, value = mx) %>% 
  mutate(type = "m_smoothed",
         year = as.integer(year)) %>% 
  # replacing missing values and estimating log_rates (/100k)
  replace_na(list(mx = 0)) %>% 
  mutate(log_mx = log(mx))

# rapid plot of smoothed rates

# defining brackets again
brks <- quantile(c(min(smt$log_mx), max(smt$log_mx)), probs = seq(0, 1, 0.25))

# labels in rates, so we estimate exposures of log_rates
lbls <- 
  round(exp(brks))

# full plot of smoothed rates
p_sm <- 
  smt %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k", 
                     breaks = brks, labels = lbls) +
  geom_contour(bins = 30, col="black", size=.15, alpha=0.8)+
  lexis_shape
p_sm

# Comparison actual vs smoothed rates 
# in the same plot using patchwork package
p1+p_sm

# inspection of the smoothing
comp <- 
  dt %>% 
  left_join(smt %>% select(-type) %>% select(year, age, log_mx_sm = log_mx))

# in different years
yrs <- c(1960, 1980, 2000, 2010)
comp %>% 
  filter(year %in% yrs) %>% 
  ggplot()+
  # scatter plot of observed mortality 
  geom_point(aes(age, log_mx, col = factor(year), group = year),
             size = 1, alpha = 0.7)+
  # line of smoothed mortality
  geom_line(aes(age, log_mx_sm, col = factor(year), group = year))+
  scale_x_continuous(expand = c(0,0), breaks = seq(amin, amax, 10)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Age", y = "Death Rates (Log)")+
  theme_bw()

# at different ages
ags <- c(0, 10, 50, 60, 80)
comp %>% 
  filter(age %in% ags) %>% 
  ggplot()+
  # scatter plot of observed mortality 
  geom_point(aes(year, log_mx, col = factor(age), group = age),
             size = 1, alpha = 0.7)+
  # line of smoothed mortality
  geom_line(aes(year, log_mx_sm, col = factor(age), group = age))+
  scale_x_continuous(expand = c(0,0), breaks = seq(pmin, pmax, 10)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Period", y = "Death Rates (Log)")+
  theme_bw()



# we can have a quick estimation of the excess mortality from the civil war
# same principle of interpolation, but the weights are also a matrix
cw_dt <- 
  dt %>% 
  filter(year %in% 1920:1960)

# Again, we need to give the information in matrix forms (age x period)
ylist_cw <- unique(cw_dt$year) %>% sort() # all periods in data, ordered
alist_cw <- unique(cw_dt$age) %>% sort() # all ages in data, ordered

# mortality matrix
deaths <- 
  matrix(cw_dt$dts, 
         nrow = length(alist_cw), 
         ncol = length(ylist_cw), 
         byrow = F)
colnames(deaths) <- ylist_cw
rownames(deaths) <- alist_cw

# lets have a look of the deaths matrix 
deaths

# population at risk (population/exposure) matrix
exposure <- 
  matrix(cw_dt$pop, 
         nrow = length(alist_cw),
         ncol = length(ylist_cw), 
         byrow = F)
colnames(exposure) <- ylist_cw
rownames(exposure) <- alist_cw

ws_m <- 
  cw_dt %>% 
  mutate(w = ifelse(year %in% 1935:1944, 0, 1))

ws <- 
  matrix(ws_m$w, 
         nrow = length(alist_cw),
         ncol = length(ylist_cw), 
         byrow = F)
colnames(ws) <- ylist_cw
rownames(ws) <- alist_cw

fit_2d_int <- 
  Mort2Dsmooth(x = alist_cw, 
               y = ylist_cw, 
               Z = deaths, 
               W = ws,
               offset = log(exposure),
               overdispersion = TRUE, 
               method = 2)

# extracting and transforming estimates
fit_2d_int$logmortality
# what units are these?

# from matrix to tidy (long) form
smt_int <- 
  (exp(fit_2d_int$logmortality) * 1e5) %>% 
  # keeping the row names as age values
  as_tibble(rownames = "age") %>%
  mutate(age = age %>% as.double()) %>% 
  # reshaping to long form (tidy)
  gather(-age, key = year, value = mx) %>% 
  mutate(type = "m_smoothed",
         year = as.integer(year)) %>% 
  # replacing missing values and estimating log_rates (/100k)
  replace_na(list(mx = 0)) %>% 
  mutate(log_mx = log(mx))

# quick look
redo_lexis_shape(pmin = 1920,
                 pmax = 1960,
                 amin = 0,
                 amax = 100)
p_cw <- 
  cw_dt %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k", 
                     breaks = brks, labels = lbls) +
  geom_contour(bins = 30, col="black", size=.15, alpha=0.8)+
  lexis_shape

p_sm_i <- 
  smt_int %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k", 
                     breaks = brks, labels = lbls) +
  geom_contour(bins = 30, col="black", size=.15, alpha=0.8)+
  lexis_shape

p_cw + p_sm_i

exc_cw <- 
  cw_dt %>% 
  left_join(smt_int %>% 
              select(age, year, mx_int = mx, log_mx_int = log_mx)) %>% 
  mutate(exc = mx/mx_int)

exc_cw %>%
  ggplot(aes(x = year, y = age, z = exc))+
  geom_tile(aes(fill = exc)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k") +
  lexis_shape

# at different ages
ags <- c(0, 10, 50, 60, 80)
exc_cw %>% 
  filter(age %in% ags) %>% 
  ggplot()+
  # scatter plot of observed mortality 
  geom_point(aes(year, log_mx, col = factor(age), group = age),
             size = 1, alpha = 0.7)+
  # line of smoothed mortality
  geom_line(aes(year, log_mx_int, col = factor(age), group = age))+
  scale_x_continuous(expand = c(0,0), breaks = seq(pmin, pmax, 10)) +
  scale_y_continuous(expand = c(0,0))+
  labs(x = "Period", y = "Death Rates (Log)")+
  theme_bw()


smt_int_fit <- 
  predict(fit_2d_int, se.fit = TRUE)$fit %>% 
  as_tibble(rownames = "age") %>%
  # reshaping to long form (tidy)
  gather(-age, key = year, value = log_mx) %>% 
  mutate(type = "m_smoothed",
         year = as.integer(year)) %>% 
  # replacing missing values and estimating log_rates (/100k)
  replace_na(list(log_mx = 0))

smt_int_se <- 
  predict(fit_2d_int, se.fit = TRUE)$se.fit %>% 
  as_tibble(rownames = "age") %>%
  # reshaping to long form (tidy)
  gather(-age, key = year, value = se) %>% 
  mutate(type = "m_smoothed",
         year = as.integer(year)) %>% 
  # replacing missing values and estimating log_rates (/100k)
  replace_na(list(se = 0))

smt_int_ci <- 
  smt_int_fit %>% 
  left_join(smt_int_se) %>% 
  mutate(mx_int = exp(log_mx)*1e5,
         ll = exp(log_mx-1.96*se)*1e5,
         ul = exp(log_mx+1.96*se)*1e5,
         age = age %>% as.integer()) %>% 
  select(age, year, mx_int, ll, ul)

fit_2d_int

exc_cw <- 
  cw_dt %>% 
  left_join(smt_int_ci) %>% 
  mutate(exc = ifelse(dts > ul, mx/mx_int, 0))

exc_cw %>%
  ggplot(aes(x = year, y = age, z = exc))+
  geom_tile(aes(fill = exc)) +
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k") +
  lexis_shape

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rates of mortality change ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# mortality change over periods
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Percentage of mortality change within the same age ((year t)/(year t-1) - 1) * 100
db_per <- 
  smt2 %>% 
  group_by(age) %>% 
  mutate(ch = ((mx / lag(mx)) - 1) * 100,
         ch2 = mx/lag(mx)) %>% 
  ungroup() %>% 
  drop_na()

######

# Lexis surface of mortality change
# ggplot has an option for divergent scales (scale_fill_gradient2)
redo_lexis_shape(pmin, pmax, amin, amax)
db_per %>% 
  ggplot(aes(x = year, y = age, z = ch)) +
  geom_tile(aes(fill = ch))+
  geom_contour(breaks = 0, col = "black", size = .15, alpha = 0.8)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, space = "Lab", na.value = "grey50", 
                       guide = "colourbar", aesthetics = "fill",
                       name = "Mortality\nchange %")+
  lexis_shape
  
# but we do not have control of the brackets, scale, etc.
# e.g., ggplot does not allow to create a divergent scale with four colors.
# we can make a personal color palette and legend bar with the function ColorRampalette



# Adding more control to Lexis surface of mortality change
##########################################################

# Color palette with two color scales (25 levels each), 
qt <- 25
# one for decrease of mortality (from green to blue), 
# another for mortality increases (from yellow to red)
col_scale <- 
  c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
    "black",
    colorRampPalette(c("yellow", "red"), space = "Lab")(qt))

# Definition of brackets for the scale of change
val <- unique(db_per$ch)
val <- log(unique(db_per$ch2))
# separate negative and positive values
pval <- val[val>0] # all positive values of change (mortality deterioration)
nval <- val[val<0] # all negative values of change (mortality improvement)

# identification of brackets for positive values (minimum + 23 quantiles + maximum)
# pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
pcop <- quantile(pval, prob = seq(0, 1, 1/qt))
pcop[length(pcop)] <- pcop[length(pcop)]*1.01
# the same as above but for negative values (minimum + 23 quantiles + maximum)
# ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
ncop <- quantile(nval, prob = seq(0, 1, 1/qt))
ncop[1] <- ncop[1]*1.01
# chain of brackets 25 ranges for negative changes, central value of no change (0), 
# and 25 ranges for positive changes
breaks_mc <- c(ncop, pcop) 
as.numeric(breaks_mc)
# transforming them to percentage values
(exp(breaks_mc)-1)*100

# adding to each value of change (continuous) the corresponding 
# bracket (a discrete interval)
db_per2 <- 
  db_per %>% 
  mutate(ch_cut = cut(ch, breaks = (exp(breaks_mc)-1)*100)) 

db_per2 <- 
  db_per %>% 
  mutate(ch_cut = cut(ch, breaks = (exp(breaks_mc)-1)*100)) 

cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
col_brks <- cuts[c(1, 10, 20, 25, 30, 40, 50)]

# assigning a color from our scale to each bracket
col_values <- setNames(col_scale, cuts)

# Plot of mortality change over periods
db_per2 %>%
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  # adding the color palette constructed above
  scale_fill_manual(name = "Mortality\nchange %",
                    breaks = col_brks,
                    values = col_values)+
  #adding contour lines when the slope is 0
  geom_contour(aes(z = ch), breaks = 0, col="black", alpha=.7, size=.3)+ 
  lexis_shape

# mortality change over age
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Percentage of mortality change within the same period ((age x)/(age x-1) - 1) * 100
db_age <- 
  smt2 %>% 
  group_by(year) %>% 
  mutate(ch = ((mx / lag(mx)) - 1) * 100) %>% 
  ungroup() %>% 
  drop_na()

# Definition of brackets for the scale of change
val <- unique(db_age$ch)
# separate negative and positive values
pval <- val[val>0] # all positive values of change (mortality deterioration)
nval <- val[val<0] # all negative values of change (mortality improvement)

# identification of brackets for positive values (minimum + 23 quantiles + maximum)
pcop <- quantile(pval, prob = seq(0, 1, 1/qt))
pcop[length(pcop)] <- pcop[length(pcop)]*1.01
# the same as above but for negative values (minimum + 23 quantiles + maximum)
ncop <- quantile(nval, prob = seq(0, 1, 1/qt))
ncop[1] <- ncop[1]*1.01
# chain of brackets 25 ranges for negative changes, central value of no change (0), 
# and 25 ranges for positive changes
breaks_mc <- c(ncop, 0, pcop) 

as.numeric(breaks_mc)

# adding to each value of change (continuous) the corresponding 
# bracket (a discrete interval)
db_age2 <- 
  db_age %>% 
  mutate(ch_cut = cut(ch, breaks = breaks_mc)) 

cuts <- db_age2 %>% pull(ch_cut) %>% unique() %>% sort()
col_brks <- cuts[c(1, 10, 20, 25, 30, 40, 50)]

# assigning a color from our scale to each bracket
col_values <- setNames(col_scale, cuts)

# Plot of mortality change over periods
db_age2 %>%
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  # adding the color palette constructed above
  scale_fill_manual(name = "Mortality\nchange %",
                    breaks = col_brks,
                    values = col_values)+
  #adding contour lines when the slope is 0
  geom_contour(aes(z = ch), breaks = 0, col="black", alpha=.7, size=.3)+ 
  lexis_shape


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assignment in class: ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Let's have a look at the mortality trends in different populations within the 
# same country (HMD countries).
# Group exercise:
# 1) Build a Lexis surface of mortality change for each of the following groups 
# 2) Identify the presence of period or cohort effects in each group
# 3) Discuss which possible mechanisms could have caused the observed disturbances 
# in mortality, and the differences between both populations. You can have a quick 
# look on internet to find some clues
# 4) Choose a member of the group to present in 2 minutes the Lexis surfaces, 
# the identified effects, and the mechanisms that could be driving those 
# disturbances in mortality

# Group 1
# Russia (females vs males)

# Group 2
# New Zealand (Maories vs non-Maories)

# Group 3
# German females (East vs West) (hint: Wir Kinder vom Bahnhof Zoo)

# Group 4
# USA (females vs males)

# Group 5
# Spain (females vs males)
