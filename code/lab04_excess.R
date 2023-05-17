# EDSD 2023
# Course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 4: Introduction to estimation of excess mortality

rm(list=ls())
library(tidyverse)
library(lubridate)
library(ISOweek)

# loading mortality data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# from STMF
# ~~~~~~~~~
# downloading the last version of STMF Mortality input data zip 
# hint: input STMF is much better than the output because it comes in 
# finer age groups resolution (usually 5-y)! 

# download.file("https://www.mortality.org/File/GetDocument/Public/STMF/Inputs/STMFinput.zip", 
#               "data_input/STMFinput.zip")

# list of country codes in STMF
zip_files <- unzip("data_input/STMFinput.zip", list = TRUE)

# lets look at Spanish data
cd <- "ESP"
spain_file <- zip_files %>% filter(str_detect(Name, cd)) %>% pull(Name)

esp_dt <- 
  read_csv(unz("data_input/STMFinput.zip", spain_file)) %>% 
  mutate(Week = as.double(Week))

esp_dt
unique(esp_dt$Type)
unique(esp_dt$Access)
unique(esp_dt$Age)

# adding date to each ISO week, using the package ISOweek
# we need the week in this format: 2000-W01-7
esp_dt2 <- 
  esp_dt %>% 
  rename_with(str_to_lower) %>% 
  mutate(isoweek = paste0(year, "-W", sprintf("%02d", week), "-7"),
         date = ISOweek2date(isoweek)) %>% 
  filter(age == "TOT",
         sex == "b",
         year >= 2010) %>% 
  select(popcode, year, week, dts = deaths, date)

esp_dt2 %>% 
  ggplot()+
  geom_line(aes(date, dts))+
  theme_bw()

ggsave("figures/weekly_mortality_spain.png",
       w = 8, h = 2)


# ######################################
# estimating the baseline mortality ====
# ######################################


# Yearly average-week ==== 
# ~~~~~~~~~~~~~~~~~~~~~~~~
# the simplest approach (maybe the worst)

# estimating the average weekly deaths for the whole training period 
# as the baseline
yw_av <- 
  esp_dt2 %>% 
  filter(year <= 2019) %>%
  summarise(bsn = mean(dts)) %>% 
  pull(bsn)

esp_yw_av <- 
  esp_dt2 %>% 
  mutate(bsn = yw_av)
  
esp_yw_av %>% 
  ggplot()+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "#8ac926")+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()

ggsave("figures/weekly_mortality_spain_year_av_week_bsn.png",
       w = 8, h = 2)


# Week-specific average ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
wk_av <- 
  esp_dt2 %>% 
  filter(year <= 2019) %>%
  group_by(week) %>% 
  summarise(bsn = mean(dts)) %>% 
  ungroup()

esp_wk_av <- 
  esp_dt2 %>% 
  left_join(wk_av)

esp_wk_av %>% 
  ggplot()+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "#1982c4")+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()

ggsave("figures/weekly_mortality_spain_av_week_bsn.png",
       w = 8, h = 2)

# it seems not bad at all!!

# # let's estimate the excess
# 
# # weekly excess deaths
# exc_av <- 
#   esp_wk_av %>% 
#   filter(date >= "2020-03-01") %>% 
#   mutate(exc = dts-bsn)
# 
# exc_av %>% 
#   ggplot()+
#   geom_line(aes(date, exc))+
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   theme_bw()
# 
# # estimating total excess by year
# tot_exc_av <- 
#   exc_av %>% 
#   group_by(year) %>% 
#   summarise(dts = sum(dts),
#             bsn = sum(bsn),
#             exc = sum(exc)) %>% 
#   ungroup() %>% 
#   mutate(psc = dts / bsn - 1)


# Excess mortality using a Poisson model ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# exposures
# ~~~~~~~~~
# we need: weekly deaths and exposures
# issue: population counts only available annually
# solution: interpolation

# loading total population counts from WPP 
pop <- read_rds("data_input/wpp2022_pop.rds")

# selecting the Spanish population
esp_pop <- 
  pop %>% 
  filter(name == "Spain")
  
# from annual to weekly exposures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wee need weekly exposures between 2010 and 2023
# then, we can interpolate between 2009 and 2024

# first, create a grid with empty values
pop_empty <- 
  # 52 weeks per year
  expand_grid(year = 2009:2024, week = 1:52) %>% 
  # adding extra weeks for leap years 2009, 2015, and 2020
  bind_rows(tibble(year = c(2009, 2015, 2020), week = 53)) %>% 
  # arrange it chronologically
  arrange(year, week)

pop_inter <- 
  pop_empty %>%  
  left_join(esp_pop %>%
              select(year, pop) %>% 
              mutate(week = 26)) %>% 
  mutate(t = 1:n(),
         isoweek = paste0(year, "-W", sprintf("%02d", week), "-7"),
         date = ISOweek2date(isoweek))

# extracting values for interpolation
xs <- pop_inter %>% drop_na() %>% pull(t) # weeks with data
ys <- pop_inter %>% drop_na() %>% pull(pop) # available pop data
ts <- pop_inter %>% pull(t) # all weeks in which we need estimates

# smoothing using cubic splines
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the "spline()" function allows to interpolate by constraining the curve 
# to match the original data. In other words, it is strictly interpolation 
# rather than smoothing  

pop_inter2 <- 
  pop_inter %>% 
  mutate(pop2 = spline(xs, ys, xout = ts)$y)

# visualizing annual values and weekly estimates
pop_inter2 %>% 
  ggplot()+
  geom_line(aes(date, pop2))+
  geom_point(aes(date, pop), col = "red")+
  theme_bw()

# weekly exposures to use
pop2 <- 
  pop_inter2 %>% 
  select(-pop, -t, -isoweek) %>% 
  rename(pop = pop2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# merging deaths and exposures
esp_dt3 <- 
  esp_dt2 %>% 
  left_join(pop2) %>% 
  mutate(exposure = pop / 52, # exposure in person-weeks
         t = 1:n(), # a variable for the secular trend 
         w = ifelse(date <= "2020-03-15", 1, 0)) # a variable for the weights 


# GAM models allow us to include parametric and semiparametric terms together 
# the "mgcv" package is very convenient for working with GAM models
# https://cran.r-project.org/web/packages/mgcv/mgcv.pdf

library(mgcv)

# fitting the model
gam_model <- 
  gam(dts ~
        # linear term (exponential outside Poisson) for the secular trend
        t +
        # cyclic spline term for the seasonal trend
        s(week, bs = 'cp') + 
        # controlling for population changes over time
        offset(log(exposure)), 
      # to avoid including the pandemic period in the baseline estimation
      weights = w, 
      data = esp_dt3, 
      # using a quasipoisson distribution to account for overdispersion
      family = "quasipoisson") 

# how does it look?
summary(gam_model)

# weekly slope (t) = 0.00018984
# how much is that in 10 years?
((1 + 0.00018984)^52)^10

# example for predicting estimates
predict(gam_model, newdata = esp_dt3, type = "response", se.fit = TRUE)
predict(gam_model, newdata = esp_dt3)

# obtaining estimates for the three models
esp_bsn <- 
  esp_dt3 %>% 
  mutate(bsn = predict(gam_model, 
                       type = "response"),
         se = predict(gam_model, 
                      newdata = esp_dt3, 
                      type = "response", 
                      se.fit = TRUE)$se.fit,
         ll = bsn - 1.96*se,
         ul = bsn + 1.96*se)

esp_bsn %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = ll, ymax = ul), alpha = 0.5, fill = "red")+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn), col = "red")+
  theme_bw()

# excess estimation
esp_exc <- 
  esp_bsn %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31") %>% 
  mutate(exc_all = dts - bsn, # all excess
         exc_pos = ifelse(exc_all > 0, exc_all, 0), # only positive excess
         exc_sig = ifelse(dts > ul, exc_all, 0)) # only positive ans significant excess
  

# comparison
esp_exc %>%
  select(date, starts_with("exc")) %>% 
  gather(-date, key = type, value = excess) %>% 
  ggplot()+
  geom_line(aes(date, excess, col = type))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw()

# annual difference between alternative Poisson excess estimates
esp_exc_yr <- 
  esp_exc %>% 
  group_by(year) %>% 
  summarise(exc_all = sum(exc_all),
            exc_sig = sum(exc_sig),
            exc_pos = sum(exc_pos)) %>% 
  ungroup() %>% 
  mutate(diff_abs_pos = exc_pos- exc_all,
         diff_abs_sig = exc_sig- exc_all,
         diff_rel_pos = diff_abs_pos / exc_all,
         diff_rel_sig = diff_abs_sig / exc_all)

esp_exc_yr
# when excluding negative excess we add between 5% and 80% excess deaths by year

esp_exc_yr %>% 
  group_by() %>% 
  summarise(exc_all = sum(exc_all),
            exc_sig = sum(exc_sig),
            exc_pos = sum(exc_pos)) %>% 
  ungroup() %>% 
  mutate(diff_abs_pos = exc_pos- exc_all,
         diff_abs_sig = exc_sig- exc_all,
         diff_rel_pos = diff_abs_pos / exc_all,
         diff_rel_sig = diff_abs_sig / exc_all)

# 18% higher excess between 2020 and 2022 when excluding negative excess 
# estimates!! Not negligible at all

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# comparing excess estimation approaches ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comp <- 
  esp_yw_av %>% 
  select(year, week, date, dts, bsn) %>% 
  mutate(type = "weekly_average") %>% 
  bind_rows(esp_wk_av %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "weekly_spc_average")) %>% 
  bind_rows(esp_bsn %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "poisson_model"))

# plotting the three baselines
comp %>% 
  ggplot()+
  geom_line(aes(date, dts))+
  geom_line(aes(date, bsn, col = type))+
  theme_bw()

ggsave("figures/weekly_mortality_spain_baselines.png",
       w = 8, h = 2)

wk_comp <- 
  comp %>% 
  mutate(exc = dts - bsn,
         psc = dts / bsn) %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31")

# visualizing weekly excess deaths
wk_comp %>% 
  ggplot()+
  # geom_line(aes(date, exc))+
  geom_line(aes(date, exc, col = type))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw()

# visualizing cumulative excess
wk_comp %>% 
  group_by(type) %>% 
  mutate(exc_cum = cumsum(exc)) %>% 
  ggplot()+
  # geom_line(aes(date, exc))+
  geom_line(aes(date, exc_cum, col = type))+
  theme_bw()

ggsave("figures/spain_weekly_cum_excess.png",
       w = 8, h = 2)

# obtaining annual excess
yr_comp_exc <- 
  comp %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31") %>% 
  mutate(exc = dts - bsn) %>% 
  group_by(year, type) %>% 
  summarise(exc = sum(exc)) %>% 
  ungroup() %>%
  spread(type, exc)

yr_comp_diff <- 
  yr_comp_exc %>% 
  mutate(diff_abs_wa = weekly_average - poisson_model,
         diff_rel_wa = weekly_average / poisson_model,
         diff_abs_wsa = weekly_spc_average - poisson_model,
         diff_rel_wsa = weekly_spc_average / poisson_model) %>% 
  select(year, starts_with("diff"))

tt_comp_diff <- 
  yr_comp_exc %>% 
  summarise(diff_abs_wa = sum(weekly_average) - sum(poisson_model),
            diff_rel_wa = sum(weekly_average) / sum(poisson_model),
            diff_abs_wsa = sum(weekly_spc_average) - sum(poisson_model),
            diff_rel_wsa = sum(weekly_spc_average) / sum(poisson_model)) 

tt_comp_diff

# Looking at the three years 2020-2022, estimates using averages double 
# those obtained from the Poisson model!!!
