# Author: Saeesh Mangwani
# Date: 2023-05-25

# Description: Fitting GD models for all pulse hypotheses constructed from
# script 2

# ==== Libraries ====
library(dplyr)
library(readr)
library(readxl)
library(stringr)
library(purrr)
library(optimx)
# library(CatDyn)
source("scripts/_help_funcs.R")
source("scripts/CDMN.functions.R")

# ==== Global variables ====

# Current year being modelled
modyear <- 2017

# Distributions to test
dists <- c('apnormal', 'aplnormal')
dists <- setNames(dists, dists)

# Optimization methods to test
optms <- c("spg", "CG", "Nelder-Mead")

# ==== File paths ====

# Path to base output folder for this year's results
out_dir <- file.path('output/gdm', modyear)

# Path to prepared CatDyn time series dataset (from script 011_*)
yrdf_path <- list.files(out_dir, pattern = 'yrdf', full.names = T)

# Path to dataframe containing initialization parameters
param_path <- list.files(out_dir, pattern = 'param', full.names = T)

# Path to dataframe containing hypotheses identified from script 2
hypdf_path <- list.files(out_dir, pattern = 'hypotheses', full.names = T)

# Path to Rdata file where initialization conditions from the model are stored
base_rdpath <- file.path(out_dir, 'gdm_base.Rdata')

# Path to folder where model results will be saved as Rdata files
mod_dir <- file.path(out_dir, 'models')

# Creating a sub-folder for GDM plots from pulse-hypothesis models
gdm_plot_dir <- file.path(out_dir, 'plots', 'gdm')
dir.create(gdm_plot_dir)

# ==== Reading data ====

# Catch/effort dataframes
yrdf <- read_csv(yrdf_path)
paramdf <- read_csv(param_path)

# Hypothesis dataframe
hypdf <- read_excel(hypdf_path) %>% 
  group_by(type) %>% 
  mutate(grpid = 1:n()) %>% 
  ungroup() %>% 
  # Creating unique names for each hypothesis based on type
  mutate(modname = case_when(
    type == 0 ~ 'gdm_pure',
    type < 0 ~ paste0('gdm_', abs(type), 'p', abs(type), 'p_', grpid),
    type > 0 ~ paste0('gdm_', abs(type), 'p_', grpid)
  ))

# Initial model objects - loads the base CatDyn object dataset (anm_cdyn) and
# the exploratory model fit/null model (null_model_exp)
load(base_rdpath)

# ==== Defining initial parameters ====

# Setting season dates for the year
season_range <- c(paste0(modyear, '-01-02'), paste0(modyear, '-12-30'))

# Timestep range:
tsteps <- anm_cdyn$Data[[1]]$time.step

# Getting tuned initial parameters from the dataframe
params <- setNames(paramdf$value, paramdf$param)

# ==== Fitting models from pulse hypotheses ====

# Iterating over each hypothesis using its unique id
walk(hypdf$id, ~{
  # Getting the created name for this model
  modname <- hypdf %>% filter(id == .x) %>% pull(modname)
  print(paste('Hypothesis:', modname))
  
  # Creating an output Rdata path from this model name - this is where results
  # will be stored
  out_rdpath <- file.path(mod_dir, paste0(modname, '.Rdata'))
  
  # If the model object already exists, skipping
  if(file.exists(out_rdpath)) {
    print(paste0('Results already exist at: ', out_rdpath, '.Skipping...'))
    return()
  }
  
  # Getting number of pulses in this hypothesis
  pulse <- hypdf %>% filter(id == .x) %>% pull(type)
  
  # Getting pulses weeks for this hypothesis
  in_weeks <- hypdf %>% filter(id == .x) %>% pull(in_week)
  out_weeks <- hypdf %>% filter(id == .x) %>% pull(out_week)
  # Converting to numeric vectors
  in_weeks <- as.numeric(str_split(in_weeks, ',')[[1]])
  out_weeks <- as.numeric(str_split(out_weeks, ',')[[1]])
  
  # Defining the ts range
  ts_range <- get_ts_range(tsteps, pulse, in_weeks, out_weeks)
  # Getting relevant initialization parameters
  pars <- get_relv_params(params, pulse)
  
  print('--- Fitting models (all distributions and optimizers):')
  # Fitting models for this hypothesis, with all distributions and optimizations
  # (this takes time)
  fits_all <- imap(dists, ~{
    # Printing the name of the distribution
    print(paste('------', .y))
    # Fitting the model - this can take a WHILE
    tryCatch({
      model <- CatDynFit(
        # CatDyn dataset object - read from the RDS file, defined in script 1
        x=anm_cdyn,
        p=pulse,
        par=log(pars),
        dates=ts_range,
        distr=.x,
        method=optms,
        itnmax=10000
      )
      return(model)
    }, error=function(e){
      print(paste('Model with distribution', .y, 'failed to fit.'))
      print(e)
      return(NULL)
    })
  })
  
  # Removing NULL models, i.e models that failed to fit
  fits <- discard(fits_all, is.null)
  
  print('--- Getting predictions')
  # Getting predictions
  preds <- imap(fits, ~{
    distname <- .y
    model <- .x
    # Getting predictions from each optimization method
    outpreds <- map(optms, ~{
      print(paste0(' ------ ', distname, '-', .x))
      out <- tryCatch(CatDynPred(x=model, method=.x), error=function(e) return(NULL))
      return(get0('out'))
    })
    # Setting names
    outpreds <- setNames(outpreds, optms)
    # Removing null results - i.e where the model failed to converge
    outpreds <- discard(outpreds, is.null)
    # Returning the converged model list
    return(outpreds)
  })
  
  print('--- Plotting predictions')
  # Plotting prediction results and saving to disk
  iwalk(preds, ~{
    dist_name <- .y
    fits <- .x
    iwalk(fits, ~{
      # Creating filename for this plot
      fname <- paste0(modname, '_', dist_name, '_', .y, '.png')
      # print(fname)
      # Creating plot and saving to disk
      png(file=file.path(gdm_plot_dir, fname))
      qplot_cdyn(.x, yrdf, Cat.tstep=dim(yrdf)[1])
      dev.off()
    })
  })
  
  print('--- Saving results')
  # Combining the models with their predictions
  models <- imap(fits, ~{
    .x$Predictions <- preds[[.y]]
    return(.x)
  })
  # Writing to Rdata file
  save(models, file=out_rdpath)
})



