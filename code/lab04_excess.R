# EDSD 2023
# Course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 4: Introduction to estimation of excess mortality

rm(list=ls())
source("code/lab00_prep_session.R")

# 3 baselines: 3 colors
cols <- brewer.pal(3, "Dark2")


# loading mortality data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# from STMF
# ~~~~~~~~~
# downloading the last version of STMF Mortality input data zip 
# hint: input STMF is much better than the output because it comes in 
# finer age groups resolution (usually 5-y)! 

# list of country codes in STMF
zip_files <- unzip("data_input/STMFinput.zip", list = TRUE)

zip_files

# lets look at Spanish data, total sex, all ages, since 2015
cd <- "ESP"
sx <- "b"
ag <- "TOT"
ymin <- 2015

file_name <- zip_files %>% filter(str_detect(Name, cd)) %>% pull(Name)

dt <- 
  read_csv(unz("data_input/STMFinput.zip", file_name)) %>% 
  mutate(Week = as.double(Week))

dt

# adding date to each ISO week, using the package ISOweek
# we need the week in this format: 2000-W01-7
dt2 <- 
  dt %>% 
  # renaming all variablees to lower case
  rename_with(str_to_lower) %>% 
  # adding the date
  mutate(isoweek = paste0(year, "-W", sprintf("%02d", week), "-7"),
         date = ISOweek2date(isoweek)) %>% 
  filter(age == ag,
         sex == sx,
         year >= ymin) %>% 
  select(code = popcode, year, week, dts = deaths, date)

dt2 %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  theme_bw()



# ######################################
# estimating the baseline mortality ====
# ######################################


# average-week ====
# ~~~~~~~~~~~~~~~~~~~~~~~~
# the simplest approach (maybe the worst)

# estimating the average weekly deaths for the whole training period
# as the baseline
w_av <- 
  dt2 %>% 
  filter(year <= 2019) %>%
  summarise(bsn = mean(dts)) %>% 
  pull(bsn)

dt_w_av <- 
  dt2 %>% 
  mutate(bsn = w_av)

dt_w_av %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[1])+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()

# Week-specific average ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
ws_av <- 
  dt2 %>% 
  filter(year <= 2019) %>%
  group_by(week) %>% 
  summarise(bsn = mean(dts)) %>% 
  ungroup()

dt_ws_av <- 
  dt2 %>% 
  left_join(ws_av)

dt_ws_av %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[2])+
  geom_vline(xintercept = ymd("2020-03-15"), linetype = "dashed")+
  theme_bw()

# it seems not bad at all!!


# Excess mortality using a Poisson model ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# we need: weekly deaths and exposures

# exposures
# ~~~~~~~~~
# issue: population counts only available annually
# solution: interpolation

# loading total population counts from WPP 
pop <- read_rds("data_input/wpp2022_pop.rds")

# selecting the Spanish population
pop2 <- 
  pop %>% 
  filter(code == cd) %>%
  select(year, pop) %>% 
  mutate(week = 26)

# from annual to weekly exposures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# wee need weekly exposures between 2010 and 2023
# then, we can interpolate between 2009 and 2024

# first, create a grid with all weeks that we need to populate with 
# weekly exposures
pop_empty <- 
  # 52 weeks per year
  expand_grid(year = 2009:2024, week = 1:52) %>% 
  # adding extra weeks for leap years 2009, 2015, and 2020
  bind_rows(tibble(year = c(2009, 2015, 2020), week = 53)) %>% 
  # arrange it chronologically
  arrange(year, week)

# preparing data for interpolation
pop_inter <- 
  pop_empty %>%  
  # adding annual population in midyear to week 26 (the midyear week!)
  left_join(pop2) %>% 
  mutate(t = 1:n(), # creating a continuous variable for week sequence 
         # transforming year-week to date, using ISOweek2date() function
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

# extracting predictions from the model
interpolated_pop <- spline(xs, ys, xout = ts)

pop_inter2 <- 
  pop_inter %>% 
  mutate(pop2 = interpolated_pop$y)

# visualizing annual values and weekly estimates
pop_inter2 %>% 
  ggplot()+
  geom_line(aes(date, pop2), linewidth = 1)+
  geom_point(aes(date, pop), col = "red", size = 3)+
  theme_bw()

# weekly exposures to use
pop3 <- 
  pop_inter2 %>% 
  select(-pop, -t, -isoweek) %>% 
  rename(pop = pop2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# merging deaths and exposures
dt3 <- 
  dt2 %>% 
  left_join(pop3) %>% # merging with weekly population
  mutate(exposure = pop / 52, # exposure in person-weeks
         t = 1:n(), # a variable for the secular trend 
         w = ifelse(date <= "2020-03-01", 1, 0)) # a variable for the weights 


# GAM models allow us to include parametric and semiparametric terms together 
# the "mgcv" package is very convenient for working with GAM models
# https://cran.r-project.org/web/packages/mgcv/mgcv.pdf

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
      data = dt3, 
      # using a quasipoisson distribution to account for overdispersion
      family = "quasipoisson") 

# how does it look?
summary(gam_model)

# weekly slope (t) = 0.00018984
# how much is that in 10 years?
((1 + 0.00018984)^52)^10

# example for predicting estimates
bsn_poi <- predict(gam_model, newdata = dt3, type = "response", se.fit = TRUE)

# obtaining estimates for the three models
bsn <- 
  dt3 %>% 
  mutate(bsn = bsn_poi$fit,
         se = bsn_poi$se.fit,
         ll = bsn - 1.96*se,
         ul = bsn + 1.96*se)

bsn %>% 
  ggplot()+
  geom_ribbon(aes(date, ymin = ll, ymax = ul), alpha = 0.5, fill = cols[3])+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn), linewidth = 1, col = cols[3])+
  theme_bw()

# excess estimation
exc <- 
  bsn %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31") %>% 
  mutate(exc = dts - bsn) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# comparing excess estimation approaches ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# binding together the three estimates

bsns <- 
  bind_rows(dt_w_av %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "weekly_average"),
            dt_ws_av %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "weekly_spc_average"),
            bsn %>% 
              select(year, week, date, dts, bsn) %>% 
              mutate(type = "poisson_model"))

# plotting the three baselines
p_bsns <- 
  bsns %>% 
  ggplot()+
  geom_line(aes(date, dts), linewidth = 1)+
  geom_line(aes(date, bsn, col = type), linewidth = 1)+
  scale_color_manual(values = cols)+
  theme_bw()
p_bsns

excs <- 
  bsns %>% 
  mutate(exc = dts - bsn,
         psc = dts / bsn) %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31")

# visualizing weekly excess deaths
excs %>% 
  ggplot()+
  geom_line(aes(date, exc, col = type), linewidth = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_manual(values = cols)+
  theme_bw()

# visualizing cumulative excess
excs_cum <- 
  excs %>% 
  group_by(type) %>% 
  mutate(exc_cum = cumsum(exc)) %>% 
  ggplot()+
  # geom_line(aes(date, exc))+
  geom_line(aes(date, exc_cum, col = type))+
  theme_bw()

# obtaining annual excess
# yr_exc <- 
#   excs %>% 
#   filter(date >= "2020-03-15",
#          date <= "2022-12-31") %>% 
#   group_by(year, type) %>% 
#   summarise(exc = sum(exc)) %>% 
#   ungroup() %>%
#   spread(year, exc) %>% 
#   mutate(total = `2020` + `2021` + `2022`)
# 
# yr_exc

yr_exc <- 
  excs %>% 
  filter(date >= "2020-03-15",
         date <= "2022-12-31") %>% 
  group_by(year, type) %>% 
  summarise(exc = sum(exc)) %>% 
  ungroup()

yr_exc %>% 
  spread(year, exc)


tots <- 
  yr_exc %>% 
  spread(type, exc) %>% 
  summarise(poisson_model = sum(poisson_model),
            weekly_average = sum(weekly_average),
            weekly_spc_average = sum(weekly_spc_average),
            year = "total")

tots

t1 <- round(tots$poisson_model)
t2 <- round(tots$weekly_average)
t3 <- round(tots$weekly_spc_average)

d2 <- paste0(round((t2 / t1 - 1)*100, 1), "%")
d3 <- paste0(round((t3 / t1 - 1)*100, 1), "%")

# plotting excess estimates by year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yr_exc %>% 
  ggplot(aes(type, exc, fill = factor(year))) + 
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(x = "poisson_model", y = 0), label = paste0(t1, "\nreference"), vjust = -1)+
  geom_text(aes(x = "weekly_average", y = 0), label = paste0(t2, "\n", d2), vjust = -1)+
  geom_text(aes(x = "weekly_spc_average", y = 0), label = paste0(t3, "\n", d3), vjust = -1)+
  labs(fill = "year")+
  coord_cartesian(expand = 0)+
  theme_bw()


# Looking at the three years 2020-2022, estimates using averages double 
# those obtained from the Poisson model!!!


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assignment in class: ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimate excess mortality in the US and France using the three methods and calculate 
# the potential bias when using the average approach

excess <- 
  obtain_excess(cd = "USA", sx = "b", ag = "TOT", ymin = 2015)

excess[[1]]
excess[[2]]
excess[[3]] %>% spread(type, exc)

excess <- 
  obtain_excess(cd = "FRA", sx = "b", ag = "TOT", ymin = 2015)

excess[[1]]
excess[[2]]
excess[[3]] %>% spread(type, exc)

# ok, that was quick...

# have fun in Euskal Erkidegoa!!
# thanks for all, you are a great group, it was great to give this course :)

