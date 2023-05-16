# EDSD 2023
# Course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Preparing environment

# installing missing packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
libs <- c("devtools",
          "tidyverse", 
          "remotes",
          "readr",
          "haven", 
          "viridisLite", 
          "viridis", 
          "mgcv",
          "HMDHFDplus",
          "ISOweek",
          "patchwork",
          "RColorBrewer")

for (i in libs){
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}
options(timeout = 600)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("timriffe/MortalitySmooth")
# if (!("wpp2022" %in% rownames(installed.packages()))) remotes::install_github("PPgp/wpp2022", force = TRUE)

# Loading required packages 
lapply(libs, require, character.only = T)
library("MortalitySmooth")
# library("wpp2022")

# avoiding scientific notation
options(scipen=999)
# let's keep the same seed for reproducibility of results
set.seed(2019) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# population data from WPP 2022 ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading population data from WPP 2022 estimates

if (!file.exists("data_input/wpp2022_pop.rds")) {
  
  data(pop1dt)
  data(popproj1dt)
  
  pop_hist <- 
    pop1dt %>% 
    as_tibble() %>% 
    # filter(name =='Spain') %>% 
    select(name, year, pop, popM, popF)
  
  pop_proj <- 
    popproj1dt %>% 
    as_tibble() %>% 
    # filter(name =='Spain') %>% 
    select(name, year, pop, popM, popF)
  
  pop <- 
    bind_rows(pop_hist, pop_proj) %>% 
    arrange(name, year)
  
  # write_rds(pop, "data_input/total_annual_population_all_countries.rds")
  write_rds(pop, "data_input/wpp2022_pop.rds")
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deaths and exposures from the HMD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/hmd_dts_pop.rds")) {
 
  # HMD user and password
  # usethis::edit_r_environ()
  # hmd_us="acosta@demogr.mpg.de"
  # hmd_pw="Secreto_1"
  
  # getting HMD username and password from the R environment
  hmd_us <- Sys.getenv("hmd_us")
  hmd_pw <- Sys.getenv("hmd_pw")
  
  
  cat("Downloading deaths and exposures from HMD\n")
  cds_hmd <- getHMDcountries() %>% pull(CNTRY)
  
  hmd <- tibble()
  for(ct in cds_hmd){
    cat(paste0(ct, "\n"))
    chunk_d <- 
      readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
      as_tibble() %>%
      mutate(Code = ct)
    
    hmd <- 
      hmd %>%
      bind_rows(chunk_d)
  }
  
  hmd
  
  hmd_e <- tibble()
  for(ct in cds_hmd){
    cat(paste0(ct, "\n"))
    chunk_e <- 
      readHMDweb(ct, "Exposures_1x1", hmd_us, hmd_pw) %>%
      as_tibble() %>%
      mutate(Code = ct)
    
    hmd_e <- 
      hmd_e %>%
      bind_rows(chunk_e)
  }
  
  hmd2 <- 
    hmd %>% 
    rename_all(tolower) %>% 
    select(-openinterval) %>% 
    gather(female, male, total, key = sex, value = dts)
  
  hmd_e2 <- 
    hmd_e %>% 
    rename_all(tolower) %>% 
    select(-openinterval) %>% 
    gather(female, male, total, key = sex, value = pop)
  
  hmd_all <- 
    hmd2 %>% 
    left_join(hmd_e2) 
  
  write_rds(hmd_all, "data_input/hmd_dts_pop.rds")
  
}

# other HMD data?
# New Zealand Maories and non Maories 
# Germany East and West


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-specific fertility rates from the HFD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/hfd_asfr.rds")) {
  
  hmd_us <- Sys.getenv("hmd_us")
  hmd_pw <- Sys.getenv("hmd_pw")
  
  # some countries having issues when trying to download... 
  # excluding them for now
  cts_issues <- c("BGR", "CAN", "POL", "KOR", "RUS", "CHE")
  
  cat("Downloading ASFR from HFD\n")
  
  cds_hfd <- 
    getHFDcountries() %>% 
    filter(!CNTRY %in% cts_issues) %>% 
    pull(CNTRY)
  
  hfd <- tibble()
  for(ct in cds_hfd[c(1:3, 6:43)]){
    cat(paste0(ct, "\n"))
    
    chunk_f <- 
        readHFDweb(ct, "asfrRR", hmd_us, hmd_pw) %>%
        as_tibble() %>%
        mutate(Code = ct)
    
    hfd <- 
      hfd %>%
      bind_rows(chunk_f) %>% 
      unique()
  }
  hfd
  
  hfd2 <- 
    hfd %>% 
    rename_all(tolower) %>% 
    select(-openinterval)
  
  write_rds(hfd2, "data_input/hfd_asfr.rds")
  
}

# # parity data
# # ~~~~~~~~~~~
# if (!file.exists("data_input/hfd_parity.rds")) {
#   
#   cts_par <- c("ESP", "SWE", "JPN", "ITA", "DEUTE", "DEUTW")
#   
#   hfd_par <- tibble()
#   for(ct in cts_par){
#     cat(paste0(ct, "\n"))
#     
#     chunk_p <- 
#       readHFDweb(ct, "asfrRRbo", hmd_us, hmd_pw) %>%
#       as_tibble() %>%
#       mutate(Code = ct)
#     
#     hfd_par <- 
#       hfd_par %>%
#       bind_rows(chunk_p) %>% 
#       unique()
#   }
#   hfd_par
#   
#   hfd_par2 <- 
#     hfd_par %>% 
#     select(-OpenInterval)
#   
#   write_rds(hfd_par2, "data_input/hfd_parity.rds")
#   
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# redefining the Lexis shape specs ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
redo_lexis_shape <- 
  function(pmin, pmax, amin, amax){
    lexis_shape <<-  
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
  }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot a Lexis surface of mortality change over time ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_change <- function(c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  db2 <- db %>% 
    mutate(deaths = deaths + 1,
           Mx = 100000 * deaths / exposure, 
           log_m = log(Mx)) %>% 
    filter(country == c,
           sex == s,
           age >= amin & age <= amax, 
           year >= ymin & year <= ymax)
  
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$deaths, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$exposure, nrow = length(alist), ncol = length(ylist), byrow=F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  
  # smoothing mortality with best AIC (that is method=2, also possible best BIC with method=1)
  fit <- 
    Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                 overdispersion = TRUE, method = 2)
  
  # transforming them to smoothed deaths
  mx_smooth <- (exp(fit$logmortality) * 100000)
  
  # from matrix to tidy form (adjusting ages to start in 0)
  smt <- 
    mx_smooth %>% 
    as_tibble() %>% 
    rownames_to_column() %>% 
    rename(age = rowname) %>% 
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age) - 1 + amin2,
           year = as.integer(year))
  
  # replacing missing values and estimating log_rates (/100k)
  smt2 <- 
    smt %>% 
    replace_na(list(Mx = 0)) %>% 
    mutate(log_m = log(Mx))
  
  db_per <- 
    smt2 %>% 
    group_by(age) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  qt <- 25
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_per$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_per2 <- 
    db_per %>% 
    mutate(ch_cut = cut(db_per$ch, breaks = breaks_mc)) 
  
  # all categories
  cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  p_change <- 
    db_per2 %>% 
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values, 
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)], 
                      name = "Mortality\nchange %")+ 
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black", alpha=.7, size=.3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", size = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", size = .20, alpha = 0.8) +
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", size = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
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
  
  return(p_change)
}





