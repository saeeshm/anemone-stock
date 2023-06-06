# Author: Saeesh Mangwani
# Date: 2023-05-12

# Description: This is the initialization script for the 2017 anemone GDM. It
# implements the following operations:
# - Creates the weekly timestep dataset used for fitting the CatDyn GDMs
# - Fits an exploratory CatDyn null model, used for tuning initial parameters
# - Writes the finalized initialization parameters to a CSV

# ==== Libraries ====
library(dplyr)
library(readr)
library(lubridate)
# library(CatDyn)
source("scripts/CDMN.functions.r")

# ==== Global variables ====

# The year for which the stock is being modelled
modyear <- 2017

# ==== File paths ====

# Path to catch and effort data (weekly summarized)
ceff_path <- 'data/daily_ceff_sale_ts.csv'

# Path to mean weight table
indvwt_path <- 'data/indiv_wt_yd.csv'

# Path to output folder where results from this model fitting will be stored
out_dir <- file.path('output/gdm', modyear)
dir.create(out_dir)
dir.create(file.path(out_dir, 'plots'))

# Path to an RData file that can store the results from this and all subsequent
# models
rdpath <- file.path(out_dir, 'gdm_base.Rdata')

# ==== Reading data ====
ceff <- read_csv(ceff_path)
indvwt <- read_csv(indvwt_path)

# ==== Preparing combined dataset for modelling ====

# Filtering only data for the relevant year
yrdf_daily <- ceff %>% 
  mutate(year = year(date), yd = yday(date), .after=date) %>% 
  filter(year == modyear) %>% 
  # Joining individual weight data
  left_join(indvwt, by = 'yd') %>% 
  # Selecting only required columns
  select(date, year, 
         'catch' = catch_kg_tot, 
         'eff_hrs' = eff_h_tot,
         'eff_fsh' = eff_num_fisher_tot,
         weight_g)

# Summarizing by isoweek
yrdf <- yrdf_daily %>% 
  mutate(year = year(date), .after='date') %>% 
  mutate(month = month(date), .after='date') %>% 
  mutate(week = isoweek(date), .after='date') %>% 
  # If any isoweek assignment in January is 52, it refers to the isoweek from
  # the previous year. For now, just reclassifying this to be part of the first
  # week
  # mutate(week = ifelse((month < 2) & (week == 52), 1, week)) %>%
  # Alternatively, these can just be removed
  filter(!(month == 1 & week > 50)) %>%
  group_by(week) %>% 
  summarize(
    'catch' = sum(catch),
    'eff_hrs' = sum(eff_hrs),
    'eff_fsh' = sum(eff_fsh),
    'weight_kg' = mean(weight_g)/1000
  ) %>% 
  mutate(timestep = 50+(1:nrow(.)), .after='week')

# ==== Creating a CatDyn object ====

# Setting season dates for the year
season_range <- c(paste0(modyear, '-01-02'), paste0(modyear, '-12-30'))

# Creating a catdyn object
anm_cdyn <- as.CatDynData(
  x=yrdf,
  step="week",
  fleet.name="mano",
  coleff=4,
  colcat=3,
  colmbw=6,
  unitseff="nhours",
  unitscat="kg",
  unitsmbw="kg",
  nmult="thou",
  season.dates=season_range
)

# Plotting diagnostics from the base CatDyn object
png(file=file.path(out_dir, 'plots', 'base_diagnostics.png'))
plot(x=anm_cdyn, mark=TRUE, offset=c(9,10), hem="N")
dev.off()

# ==== Defining initial parameters ====

# Timestep range:
ts_range <- range(anm_cdyn$Data[[1]]$time.step)

# Initialization parameters - iteratively tune these based on results from the
# exploratory model below
params <- c(
  'M' = 0.03, 
  'N0' = 4200, 
  'k' = 0.00005, 
  'alpha' = 1.2, 
  'beta' = 1.0
)
# All parameters must be in the log scale for use in the model
log_params <- log(params)

# ==== Tuning initial parameters - Pure depletion model ====

# I.e this shows exploratory diagnostics based on directly using the above
# defined parameters, not trying to optimize the fit. Useful for an exploratory
# analysis to select good starting conditions, as this will improve results when
# we actually fit the model.

# Examining the ludic model using just initial parameters (No fit)
gdm_pure_expl <- catdynexp(
  x=anm_cdyn,
  p=0,
  par=log_params,
  dates=ts_range,
  distr="apnormal"
)

# Plotting results and diagnostics - and saving to disk
png(file=file.path(out_dir, 'plots', 'null_model_nofit.png'))
plot(x=gdm_pure_expl,
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=1.9,
     Biom.ypos=0.8,
     Cat.tstep=dim(yrdf)[1],
     Cat.xpos=1.9,
     Cat.ypos=0.7)
dev.off()

# ==== Tuning initial parameters - 1P-1P model ====

# Setting in and out pulse dates - inferred roughly from the exploratory model
pin <- 7
pout <- 37

# Timestep range:
ts_range <- c(
  head(anm_cdyn$Data[[1]]$time.step, 1),
  pin,
  pout,
  tail(anm_cdyn$Data[[1]]$time.step, 1)
)

# Initialization parameters - iteratively tune these based on results from the
# exploratory model below
params <- c(
  'M' = 0.03, 
  'N0' = 4200, 
  'p_in1' = 200,
  'p_out1' = 50,
  'k' = 0.00005, 
  'alpha' = 1.2, 
  'beta' = 1.0
)

# All parameters must be in the log scale for use in the model
log_params <- log(params)

# Examining the ludic model using just initial parameters (No fit)
gdm_1p_expl <- catdynexp(
  x=anm_cdyn,
  p=-1,
  par=log_params,
  dates=ts_range,
  distr="apnormal"
)

# Plotting results and diagnostics - and saving to disk
png(file=file.path(out_dir, 'plots', 'null_model_1p_nofit.png'))
plot(x=gdm_1p_expl,
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=1.9,
     Biom.ypos=0.8,
     Cat.tstep=dim(yrdf)[1],
     Cat.xpos=1.9,
     Cat.ypos=0.7)
dev.off()

# ==== Tuning initial parameters - 2P-2P model ====

# Setting in and out pulse dates - inferred roughly from the exploratory model
pin1 <- 14
pout1 <- 24
pin2 <- 33
pout2 <- 37

# Timestep range:
ts_range <- c(
  head(anm_cdyn$Data[[1]]$time.step, 1),
  pin1,
  pout1,
  pin2,
  pout2,
  tail(anm_cdyn$Data[[1]]$time.step, 1)
)

# Initialization parameters - iteratively tune these based on results from the
# exploratory model below
params <- c(
  'M' = 0.03, 
  'N0' = 4200, 
  'p_in1' = 200,
  'p_out1' = 50,
  'p_in2' = 200,
  'p_out2' = 50,
  'k' = 0.00005, 
  'alpha' = 1.2, 
  'beta' = 1.0
)

# All parameters must be in the log scale for use in the model
log_params <- log(params)

# Examining the ludic model using just initial parameters (No fit)
gdm_2p_expl <- catdynexp(
  x=anm_cdyn,
  p=-2,
  par=log_params,
  dates=ts_range,
  distr="apnormal"
)

# Plotting results and diagnostics - and saving to disk
png(file=file.path(out_dir, 'plots', 'null_model_2p_nofit.png'))
plot(x=gdm_2p_expl,
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=1.9,
     Biom.ypos=0.8,
     Cat.tstep=dim(yrdf)[1],
     Cat.xpos=1.9,
     Cat.ypos=0.7)
dev.off()

# ==== Writing outputs to disk ====

# Prepared weekly timestep CatDyn dataset for 2017
yrdf %>% 
  write_csv(file.path(out_dir, paste0('gdm_input_yrdf_', modyear, '.csv')))

# TUNED initial parameters (this will write the `params` object as-is. Ensure it
# contains the correct parameters before writing)
as_tibble(params) %>% 
  mutate(param = names(params), .before=value) %>% 
  write_csv(file.path(out_dir, paste0('gdm_params_', modyear, '.csv')))

# Null model object saved to the Rdata file
mod_objs <- ls()[ls() %in% c('anm_cdyn', 'gdm_pure_expl','gdm_1p_expl', 'gdm_2p_expl')]
save(list=mod_objs, file=rdpath)
