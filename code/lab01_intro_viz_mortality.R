# EDSD 2023
# course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 1: Introduction to mortality visualization

rm(list=ls())
source("code/lab00_prep_session.R")

# 1. Colors ====
# ~~~~~~~~~~~~~~

# standard palettes
# ~~~~~~~~~~~~~~~~~
# a few packages make life easier, for instance "RColorBrewer"

# display all available palettes
# (sequential, qualitative, and divergent palettes)
display.brewer.all()

# what about color-blind friendly palettes
display.brewer.all(colorblindFriendly = TRUE)

# choosing a different number of colors
display.brewer.pal(5, "Set3")

# show me the hexadecimal color codes 
brewer.pal(5, "Set3")

# other ideas
# https://colorbrewer2.org/
# https://coolors.co/
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html


# build our own palette
# ~~~~~~~~~~~~~~~~~~~~~
col1 <- "#0019FF"
col2 <- "#FF0000"

qt <- 5

scale1 <- colorRampPalette(c(col1, col2), space = "Lab")(qt)
scale1

# 2. unit scales for plotting mortality ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# It is a good practice to use a log scale when plotting mortality rates 
# at different ages, proportional changes, and ratios, among others...
hmd <- read_rds("data_input/hmd_dts_pop.rds")
head(hmd)
unique(hmd$code)

# Selecting Spain
ct <- "ESP"

dt <- 
  hmd %>% filter(code == "ESP") %>% 
  # estimating death rates and their log
  mutate(mx = 1e5*dts/pop,
         log_mx = log(mx)) %>% 
  filter(age <= 100)

# mortality rates
# ~~~~~~~~~~~~~~~
# plotting age-specific death rates between 2015 and 2020
mxs <- 
  dt %>% 
  filter(year %in% 2015:2020) %>% 
  ggplot()+
  geom_line(aes(age, mx, col = factor(year)))+
  facet_grid(~sex)+
  scale_color_manual(values = c(scale1, "black"))+
  theme_bw()

mxs

# log scale in the y axis when plotting death rates by age
# as we know increase is exponential with age
dt %>% 
  filter(year %in% 2015:2020) %>% 
  ggplot()+
  geom_line(aes(age, log_mx, col = factor(year)))+
  facet_grid(~sex)+
  scale_color_manual(values = c(scale1, "black"))+
  theme_bw()

# but how to interpret the values?

# better to use the log scale keeping the original values
# super easy in ggplot
mxs+scale_y_log10()



# 3. plotting estimates of multiplicative order ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ratios
# ~~~~~~
sexr <- 
  dt %>% 
  filter(year == 2020) %>% 
  select(age, sex, mx) %>% 
  spread(sex, mx) %>% 
  mutate(sexr = male/female,
         log_sexr = log(sexr))

p_sexr <- 
  sexr %>% 
  ggplot()+
  geom_point(aes(age, sexr))+
  # adding a line to identify the equality (ratio = 1)
  geom_hline(yintercept = 1, linetype  = "dashed")+
  theme_bw()

p_sexr
# any issues with this scale?
# by principle, opposite values should be equidistant, and the two scales your 
# are measuring should have the same magnitude
# it is different the meaning of scale when our thinking is 
# additive or multiplicative
# e.g., in the case of ratios it is a multiplicative relation

# in additive terms, the opposite of 2 is -2, the center is 0
# in multiplicative terms? what is the opposite of double? what is the center?
# what are the possible values in each direction?

# a solution for this is to use log scales
# the distance between 0.5 and 1 is the same as between 1 and 2
# the distance between 0 and 1 is the same as between 1 and infinite
# lets apply log scales
p_log_sexr <- 
  sexr %>% 
  ggplot()+
  geom_point(aes(age, log_sexr))+
  # adding a line to identify the equality (ratio = 1)
  geom_hline(yintercept = 1, linetype  = "dashed")+
  theme_bw()

p_log_sexr

# what happened? how to interpret the values? even the reference is not in the 
# proper place
# again, use the ggplot option to use the log scales while keeping the actual 
# labels
p_sexr+scale_y_log10(breaks = c(0.5, 1, 2, 4))

# the same principle works for relative risk (RR) measures!!




# 4. Lexis diagrams ====
# ~~~~~~~~~~~~~~~~~~~~~~
# the main visualization tool for demographers!!

amin <- 0
amax <- 100
pmin <- 1910
pmax <- 2020

p_lexis <- 
  ggplot()+
  # adding lines for periods
  geom_vline(xintercept = seq(pmin, pmax, 10), 
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  # adding lines for ages
  geom_hline(yintercept = seq(amin, amax, 10), 
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  theme_bw()
  
p_lexis

# lets plot the life of someone who was born in 1960
bt_coh <- 1960

p_lexis+
  geom_segment(aes(x = bt_coh, y = 0, xend = pmax, yend = pmax - bt_coh), 
               colour = "red")

# the age and period scales are different, and there is unnecessary extra space
p_lexis+
  geom_segment(aes(x = bt_coh, y = 0, xend = pmax, yend = pmax - bt_coh), 
               colour = "red")+
  # cutting extra space
  coord_equal(expand = 0)


# a few additional things
p_lexis+
  geom_segment(aes(x = bt_coh, y = 0, xend = pmax, yend = pmax - bt_coh), 
               colour = "red")+
  coord_equal(expand = 0)+
  # adding proper labels to both axis
  scale_x_continuous(breaks = seq(pmin, pmax, 10))+
  scale_y_continuous(breaks = seq(amin, amax, 10))+
  # adding axis titles
  labs(y = "Age", x = "Period")+
  # adding cohorts
  geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
              linetype = "dashed", color = "grey30", linewidth = .2, alpha = 0.8)

# saving preferences in an object
# tip: we can create a list of customizations and then apply that to the plot
lexis_shape <- 
  list(
    geom_vline(xintercept = seq(pmin, pmax, 10), 
               linewidth = 0.2, linetype = "dashed", 
               alpha = 0.8, color = "grey30"),
    geom_hline(yintercept = seq(amin, amax, 10), 
               linewidth = 0.2, linetype = "dashed", 
               alpha = 0.8, color = "grey30"),
    # adding cohorts
    geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8),
    coord_equal(expand = 0),
    # adding proper labels to both axis
    scale_x_continuous(breaks = seq(pmin, pmax, 10)),
    scale_y_continuous(breaks = seq(amin, amax, 10)),
    # adding axis titles
    labs(y = "Age", x = "Period"),
    theme_bw()
  )

ggplot()+
  lexis_shape

redo_lexis_shape(pmin = 1970,
                 pmax = 2010,
                 amin = 20,
                 amax = 80)

ggplot()+
  lexis_shape

redo_lexis_shape(pmin,pmax,amin,amax)

  
# 5. Lexis surfaces ====
# ~~~~~~~~~~~~~~~~~~~~~~
# lets plot the mortality experience in Spain between 1950 and 2020, 
# between ages 0 and 100
dt2 <- 
  dt %>% 
  filter(year %in% pmin:pmax,
         age %in% amin:amax)
  
# quick plot of log_rates
dt2 %>%
  ggplot(aes(x = year, y = age, z = mx))+
  geom_tile(aes(fill = mx))

# we can already apply what we used for the diagram
dt2 %>%
  ggplot(aes(x = year, y = age, z = mx))+
  geom_tile(aes(fill = mx))+
  lexis_shape

# lets add a nice scale color, package viridis has some nice options
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
dt2 %>%
  ggplot(aes(x = year, y = age, z = mx))+
  geom_tile(aes(fill = mx))+
  # adding a viridis scale
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     # breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  # contour lines are useful for looking at level changes over time
  geom_contour(bins = 20, col = "black", linewidth = .15, alpha = 0.8)+
  lexis_shape

# you can modify the option parameter for other nice palettes
dt2 %>%
  ggplot(aes(x = year, y = age, z = mx))+
  geom_tile(aes(fill = mx))+
  # adding a viridis scale
  scale_fill_viridis(option = "F", discrete = F,  direction = -1, 
                     # breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  # contour lines are useful for looking at level changes over time
  geom_contour(bins = 20, col = "black", linewidth = .15, alpha = 0.8)+
  lexis_shape

# what about the scales?
# can we select the breaks? can we select the labels to show?
brks <- quantile(c(min(dt2$mx), max(dt2$mx)), probs = seq(0, 1, 0.25))
lbls <- round(brks)

dt2 %>%
  ggplot(aes(x = year, y = age, z = mx))+
  geom_tile(aes(fill = mx))+
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  geom_contour(bins = 20, col = "black", size = .15, alpha = 0.8)+
  lexis_shape

# nice, but... what about the log rule for death rates?
# lets use the log then

dt2 %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx))+
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Mortality\nrate /100k") +
  geom_contour(bins = 20, col = "black", size = .15, alpha = 0.8)+
  lexis_shape
  
  
# are the scales better? how to interpret them?
s_lexis_log+
  scale_z_log10()

# ups!!!!!

# select a sample and convert it
brks <- quantile(c(min(dt2$log_mx), max(dt2$log_mx)), probs = seq(0, 1, 0.25))
lbls <- round(exp(brks))

dt2 %>%
  ggplot(aes(x = year, y = age, z = log_mx))+
  geom_tile(aes(fill = log_mx))+
  # adding axis titles
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     breaks = brks, labels = lbls,
                     name = "Mortality\nrate /100k") +
  lexis_shape

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Assignment in class: ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~

# Option 1.
# ~~~~~~~~~
# Construct a lexis surface of the sex ratio of mortality in a country of 
# your choice between 1900 and 2000
dt

# Option 2.
# ~~~~~~~~~
# Construct a lexis surface of age specific fertility rates ratios 
# between two countries of your choice with available data for all available 
# years
# fertility data, from the HFD, is here:
asfr <- 
  read_rds("data_input/hfd_asfr.rds") 
asfr
unique(asfr$code)


# in both cases:
# ~~~~~~~~~~~~~~
# use the appropriate scales, labels, and colors that allows an easy and proper
# interpretation of the ratios

# add contour lines to distinguish the direction of the ratios
# what kind of color scale is the most appropriate? 
# build a color scale if there is no a good one available

# add cohort labels to the Lexis surface

# Is there anything to interpret?
# Any additional analysis that can be done?

