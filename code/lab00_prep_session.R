# EDSD 2023
# Course in analysis of mortality disturbances 
# Instructor: Enrique Acosta (CED)
# Preparing environment

# installing missing packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
libs <- c("devtools",
          "tidyverse", 
          "lubridate",
          "remotes",
          "readr",
          "viridisLite",
          "viridis", 
          "mgcv",
          "HMDHFDplus",
          "ISOweek",
          "countrycode",
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
    arrange(name, year) %>% 
    mutate(code = countrycode(name, origin = "country.name",
                              destination = "iso3c")) %>% 
    drop_na(code)
  
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deaths and exposures from the full HMD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!file.exists("data_input/hmd_dts_pop_vcomplement.rds")) {
  
  dts_files <- 
    list.files("data_input/Deaths_1x1") %>% 
    as_tibble() %>% 
    filter(str_detect(value, "GBR_SCO|DEUTE|DEUTW|NZL_MA|NZL_NM")) %>% 
    pull(value)

  pop_files <- 
    list.files("data_input/Exposures_1x1") %>% 
    as_tibble() %>% 
    filter(str_detect(value, "GBR_SCO|DEUTE|DEUTW|NZL_MA|NZL_NM")) %>% 
    pull(value)
  
  # loading deaths from all countries in HMD
  db_d <- tibble()
  for(i in 1:length(dts_files)){
    
    txt_file <- dts_files[i]
    
    print(txt_file)
    
    temp <- 
      read_tsv(paste0("data_input/Deaths_1x1/", txt_file),
               skip = 1) %>%
      rename(var =1) %>% 
      mutate(var = str_replace_all(var, "\\.", "")) %>% 
      separate(1, c("year", "age", "female", "male", "total")) %>% 
      mutate(code = str_replace(txt_file, ".Deaths_1x1.txt", ""))
    
    db_d <- db_d %>% 
      bind_rows(temp)
  }
  
  # loading exposures from all countries in HMD
  db_p <- tibble()
  for(i in 1:length(pop_files)){
    
    txt_file <- pop_files[i]
    
    print(txt_file)
    
    temp <- 
      read_tsv(paste0("data_input/Exposures_1x1/", txt_file),
               skip = 1) %>%
      rename(var =1) %>% 
      mutate(var = str_replace_all(var, "\\.", "")) %>% 
      separate(1, c("year", "age", "female", "male", "total")) %>% 
      mutate(code = str_replace(txt_file, ".Exposures_1x1.txt", ""))
    
    db_p <- db_p %>% 
      bind_rows(temp)
  }
  
  db_d2 <- 
    db_d %>% 
    gather(male, female, total, key = sex, value = dts) %>% 
    left_join(db_p %>% 
                gather(male, female, total, key = sex, value = pop)) %>% 
    mutate(dts = dts %>% as.double(),
           pop = pop %>% as.double(),
           age = age %>% as.integer(),
           year = year %>% as.integer())
  
  write_rds(db_d2, "data_input/hmd_dts_pop_vcomplement.rds")
}

if (!file.exists("data_input/hmd_dts_pop_v2.rds")) {
  
  hmd1 <- read_rds("data_input/hmd_dts_pop.rds")
  hmd2 <- read_rds("data_input/hmd_dts_pop_vcomplement.rds")
    
  out <-  
    bind_rows(hmd1, hmd2)
  
  write_rds(out, "data_input/hmd_dts_pop_v2.rds")
  
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# weekly mortality from the STMF ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/STMFinput.zip")) {
  download.file("https://www.mortality.org/File/GetDocument/Public/STMF/Inputs/STMFinput.zip",
                "data_input/STMFinput.zip")
}


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
  db2 <- 
   hmd %>% 
    mutate(deaths = dts + 1,
           Mx = 100000 * dts / pop, 
           log_m = log(Mx)) %>% 
    filter(code == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax)
  
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$dts, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$pop, nrow = length(alist), ncol = length(ylist), byrow=F)
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
    as_tibble(rownames = "age") %>%
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age),
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
    geom_contour(aes(z = ch), breaks = 0, col="black", 
                 alpha=.7, linewidth = .3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8)+
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


# extracting coefficients from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_coeffs <- function(mod){
  tp1 <- 
    coef(summary(mod)) %>% 
    as_tibble(rownames = "coeff") %>% 
    mutate(tdim = case_when(str_detect(coeff, "\\(A\\)") ~ "Age",
                            str_detect(coeff, "\\(P\\)") ~ "Period",
                            str_detect(coeff, "\\(C\\)") ~ "Cohort",
                            str_detect(coeff, "I\\(C") ~ "Drift",
                            str_detect(coeff, "I\\(P") ~ "Drift"),
           effect = exp(Estimate),
           ll = exp(Estimate - 1.96*`Std. Error`),
           ul = exp(Estimate + 1.96*`Std. Error`)) %>% 
    separate(coeff, c("trash", "value"), sep = "\\)") %>% 
    mutate(value = value %>% as.double()) %>% 
    select(tdim, value, effect, ll, ul)
  
  ps_model <- tp1 %>% filter(tdim == "Period") %>% pull(value)
  ps_data <- unique(dt2$P)
  ps_miss <- ps_data[!(ps_data %in% ps_model)]
  
  cs_model <- tp1 %>% filter(tdim == "Cohort") %>% pull(value)
  cs_data <- unique(dt2$C)
  cs_miss <- cs_data[!(cs_data %in% cs_model)]
  
  tp2 <- 
    tp1 %>% 
    bind_rows(tibble(tdim = "Period", value = ps_miss, effect = 1, ll = 1, ul = 1),
              tibble(tdim = "Cohort", value = cs_miss, effect = 1, ll = 1, ul = 1)) %>% 
    mutate(tdim = factor(tdim, levels = c("Age", "Period", "Cohort"))) %>% 
    arrange(tdim, value) %>% 
    filter(tdim != "Drift")
  
  drift <- 
    tp1 %>% 
    filter(tdim == "Drift") %>% pull(effect)
  
  out <- 
    list("coeffs" = tp2, "drift" = drift)
  return(out)
}

# extracting the drift from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_drift <- 
  function(mdl){
  
  pcoeff <- 
    coef(summary(mdl))[,1] %>%
    as_tibble(rownames = "coeff") %>%
    filter(coeff %in% c("P", "C")) %>%
    pull(value)
  
  # this is in log scale, let's transform it in percentage value
  # change of mortality between periods (%)
  drift <- (exp(pcoeff) - 1)*100
  drift
  
  # interpretation
  cat(paste0(round(drift, 3),
             "% change in mortality each period/cohort category"))
}


# plotting Carstensen APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_carst <- 
  function(mod){
    
    age_eff <- 
      mod$Age %>% 
      as_tibble() %>% 
      mutate(dim = "Age")
    
    per_eff <- 
      mod$Per %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Period")
    
    coh_eff <- 
      mod$Coh %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Cohort")
    
    bind_rows(coh_eff, per_eff)
    
    p_a_eff <- 
      age_eff %>% 
      ggplot()+
      geom_line(aes(Age, Rate))+
      scale_y_log10()+
      facet_wrap(~dim, scales = "free_x")+
      theme_bw()
    
    p_pc_eff <- 
      bind_rows(coh_eff, per_eff) %>% 
      mutate(dim = factor(dim, levels = c("Period", "Cohort"))) %>% 
      ggplot()+
      geom_ribbon(aes(Year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3)+
      geom_line(aes(Year, RR))+
      scale_y_log10(breaks = c(0.2, 0.5, 0.7, 0.8, 1, 1.2, 1.5, 2, 4, 5))+
      scale_x_continuous(breaks = seq(1800, 2020, 10))+
      theme_bw()+
      facet_grid(~dim, scales = "free_x", space = "free_x")+
      geom_hline(yintercept = 1, linetype = "dashed")
    
    p_a_eff+p_pc_eff+ 
      plot_layout(widths = c(1, 2.5))
  }
