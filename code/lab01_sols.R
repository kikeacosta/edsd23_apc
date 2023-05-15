# EDSD 2023
# course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Lab 1: Introduction to mortality visualization
# Possible solutions to the proposed assignment 

rm(list=ls())
source("code/lab00_prep_session.R")

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



# Option 1. ====
# ~~~~~~~~~~~~~~

hmd <- read_rds("data_input/hmd_dts_pop.rds")
head(hmd)

ct <- "ESP"

dt <- 
  hmd %>% filter(code == "ESP") %>% 
  # estimating death rates and their log
  mutate(mx = 1e5*dts/pop,
         log_mx = log(mx)) %>% 
  filter(age <= 100)

# ratios
sexr2 <- 
  dt %>% 
  select(year, age, sex, mx) %>% 
  spread(sex, mx) %>% 
  mutate(sexr = male/female,
         log_sexr = log(sexr)) %>% 
  # removing strange/extreme values
  filter(!is.na(log_sexr) & log_sexr != Inf & log_sexr != -Inf)

# the easiest solution for the divergent scale is offered by ggplot2:
# scale_fill_gradient2()


# but first, define the Lexis shape
redo_lexis_shape(min(sexr2$year),
                 max(sexr2$year),
                 min(sexr2$age),
                 max(sexr2$age))

# define breaks
brks <- quantile(c(min(sexr2$log_sexr), max(sexr2$log_sexr)), 
                 probs = seq(0, 1, 0.25))
lbls <- round(exp(brks), 1)
lbls
# if the labelos are not ideal we can redefine them in the opposite direction 
# (first defining the labels, then the breaks) so they have more sense for 
# interpretation. E.g., 
lbls <- c(0.5, 1, 2, 4, 6)
brks <- log(lbs)

sexr2 %>% 
  ggplot(aes(x = year, y = age, z = log_sexr))+
  geom_tile(aes(fill = log_sexr))+
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0,
                       space = "Lab",
                       na.value = "grey50",
                       breaks = brks,
                       labels = lbls,
                       name = "Ratio")+
  geom_contour(aes(z = log_sexr), breaks = 0, col="black", alpha=.7, size=.3)+ 
  lexis_shape+
  labs(title = "Mortality sex ratios in Spain (males/females)")


# But what if we want to have full control of the color scale?
# like for instance decide to have fixed colors for the same values?

# then we can build our own divergent scale using the function colorRampPalette()

# how many scales in each direction?
qt <- 15
# 15 degrees for decrease of mortality (from green to blue), 
# another 15 for mortality increases (from yellow to red)
col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
               "white",
               colorRampPalette(c("yellow", "red"), space = "Lab")(qt))

length(col_scale)

# Definition of brackets for the scale of change
val <- unique(sexr2$log_sexr)
# separate negative and positive values
pval <- val[val>0] # all positive values of log sex ratio (higher male mortality)
nval <- val[val<0] # all negative values of log sex ratio (higher female mortality)

# identification of brackets for positive values
# same quantiles as for colors
seq(0, 1, 1/qt)
pcop <-quantile(pval, prob = seq(0, 1, 1/qt))
# the same as above but for negative values 
ncop <-quantile(nval, prob = seq(0, 1, 1/qt))

# adjusting brakets so all values can be included
pcop[length(pcop)] <- pcop[length(pcop)]*1.1
ncop[1] <- ncop[1]*1.01


# chain of n brackets  for negative log sex ratios, central value of same mortality (0), 
# and n ranges for positive log sex ratios
breaks_mc <- c(ncop, pcop) 
length(breaks_mc)

# adding to each value of change (continuous) the corresponding bracket (a discrete interval)
# trick here: exponenciate the breaks and work directly with the actual values 
# rather than with the logs, so labels are clearer and interpretation better  
sexr3 <- 
  sexr2 %>% 
  mutate(ch_cut = cut(sexr, breaks = exp(breaks_mc)))

# all categories
cuts <- sexr3 %>% pull(ch_cut) %>% unique() %>% sort()
# assigning color to each category
col_values <- setNames(col_scale, cuts)

sexr3 %>% 
  ggplot(aes(year, age, z = ch_cut)) +
  geom_tile(aes(fill = ch_cut))+
  # adding the color palette constructed above
  scale_fill_manual(values = col_values, 
                    breaks = cuts[c(1, 5, 10, 15, 20, 25, 30)],
                    name = "Ratio")+ 
  #adding contour lines when the slope is 0
  geom_contour(aes(z = log_sexr), breaks = 0, col="black", alpha=.7, size=.3)+ 
  lexis_shape+
  labs(title = "Lexis surface of mortality sex ratios in Spain",
       x="Period", y="Age")+
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 12),
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    plot.background = element_rect(fill = "white", colour = "transparent")
  )


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Option 2. ====
# ~~~~~~~~~~~~~~

# plotting fertility rates ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in addition to what we have already done, this time I am adding the labels as 
# percentage difference between both countries, instead of just ratios

asfr <- 
  read_rds("data_input/hfd_asfr.rds") 

unique(asfr$code)
cd <- c("ESP", "SWE", "JPN")

asfr2 <- 
  asfr %>% 
  filter(code %in% cd) 

pmin <- 1940
redo_lexis_shape(pmin,
                 max(asfr2$year),
                 min(asfr2$age),
                 max(asfr2$age))

asfr2 %>% 
  filter(year >= pmin) %>% 
  ggplot(aes(year, age, z = asfr))+
  geom_tile(aes(fill = asfr))+
  scale_fill_viridis(option = "C", discrete = F,  direction = -1, 
                     name = "Fertility\nrate /100k")+
  facet_wrap(~code)+
  lexis_shape


frs <- 
  asfr2 %>% 
  filter(year >= pmin,
         code %in% c("ESP", "SWE")) %>% 
  spread(code, asfr) %>% 
  mutate(fr = SWE / ESP,
         log_fr = log(fr))

redo_lexis_shape(1940, 2020, 15, 45)

bks <- c(-2,-1,0,1,2)
# labels converted to percentage values
lbs <- paste0(round((exp(bks) -1)*100), "%")

frs %>% 
  filter(age %in% 15:45,
         year %in% 1940:2020) %>% 
  ggplot(aes(year, age, z = log_fr))+
  geom_tile(aes(fill = log_fr))+
  scale_fill_gradient2(breaks = bks, labels = lbs, 
                       name = "% difference")+
  lexis_shape+
  labs(title = "Fertility rates ratios Sweden / Spain")


