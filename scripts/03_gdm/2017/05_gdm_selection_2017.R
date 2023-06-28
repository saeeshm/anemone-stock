# Author: Saeesh Mangwani
# Date: 2023-06-16

# Description: A script that selects the "best models" (two from each hypothesis
# type) from among models that survive minimum quality tests performed in script
# 4. Following this, a single best model is identified based on visual
# inspection. Outputs from the best set and the single best are summarized as R
# data files, plots and excel sheets

# ==== Libraries ====
library(dplyr)
# library(tidyr)
library(purrr)
# library(stringr)
library(ggplot2)
source("scripts/_help_funcs.R")
source("scripts/CDMN.functions.R")

# ==== Global/User variables ====

# Current year being modelled
modyear <- 2017

# ==== File paths ====

# Path to folder where model selection results for this year are stored (also
# where results frmo this script will be saved)
out_dir <- file.path('output/gdm', modyear, 'selection')

# Path to Rdata file containing the surviving model
surv_mods_path <- file.path(out_dir, 'surviving_models.Rdata')

# Folder where the results from the final selected model(s) will be stored
sel_path <- file.path(out_dir, 'best_selected')

# ==== Reading data ====

# Surviving model set
survmods <- get(load(file=surv_mods_path))

# ==== Selecting best models by AIC (per hypothesis type and distribution) ====
selmods <- survmods %>% 
  group_by(Model, Distribution) %>% 
  group_modify(~{
    .x %>% 
      arrange(AIC) %>% 
      head(10)
  }) %>% 
  ungroup()

# How many of each hypothesis selected?
table(selmods$Model)

# ==== Examining average parameter Coeff variation for each model ====

# Getting a list of all parameter values for each model
paramlist <- pmap(selmods, function(hypothesis, Method, modobj,...){
  pars <- CatDynPar(modobj, Method)
})
paramlist <- setNames(paramlist, (selmods$uid))

# Getting mean coefficient of variation across all estimates for each model
mean_cvptab <- imap_dfr(paramlist, ~{
  list(
    'uid' = as.integer(.y),
    'mean_cvp' = mean(.x$CVpCent, na.rm=T),
    'num_missing' = length(which(is.na(.x$CVpCent))),
    'prop_missing' = length(which(is.na(.x$CVpCent)))/nrow(.x),
    'which_missing' = paste(.x$Parameter[which(is.na(.x$CVpCent))], collapse=',')
  )
})

# Removing models with more than half missing a coefficient of variance estimate
selmods <- mean_cvptab %>% 
  filter(prop_missing < 0.5) %>% 
  select(uid, mean_cvp) %>% 
  # Then joining back to the selected model set to obtain the best models
  left_join(selmods, by = 'uid')

# How many models of each type remaining
table(selmods$Model)

# ==== Selecting and plotting the "Top" 2 models from each model type ====
topmods <- selmods %>% 
  group_by(Model) %>% 
  group_modify(~{
    # Top is defined as ranked based on:
    .x %>% 
      # AIC, lowest average CVP, lowest mean and median parameter correlations
      arrange(AIC, mean_cvp, abs(mn), abs(mdn), abs(sd)) %>%
      head(2)
  }) %>% 
  ungroup()

# ==== Manually examining best model fits to pick a single best ====

# Getting predicted results from each "best" model
toppreds <- pmap(topmods, function(hypothesis, Method, modobj, ...){
  modobj$Predictions[[Method]]
})
names(toppreds) <- topmods$hypothesis
names(toppreds)

# Selecting which model to examine
modnum <- 5
paste("Examining hypothesis:", names(toppreds)[modnum])

# Plotting model fit
yrdf <- topmods[modnum,]$modobj[[1]]$Data$Data$mano
qplot_cdyn(model = toppreds[[modnum]], yrdf = yrdf, Cat.tstep = dim(yrdf)[[1]])

# Checking parameters
CatDynPar(topmods[modnum,]$modobj[[1]], topmods[modnum,]$Method) %>% 
  mutate(across(c(Estimates, CVpCent), ~{round(.x, 5)}))

# Adding selection flag to best model dataframe - identifies the selected 1 or 2
# models I am choosing for this year
selindex <- c(5,6)
selids <- topmods$uid[selindex]
topmods <- topmods %>% 
  mutate(Sel.Model = ifelse(uid %in% selids, T, F))

# ==== Saving results to disk ====

# Saving top model set ----------

# Saving top model table to an Rdata file
saveobjs <- ls()[ls() %in% c('topmods')]
save(list = saveobjs, file=file.path(out_dir, 'top_models.Rdata'))

# Saving only summarizing information to a simpler excel sheet
topmods %>% 
  select(-c(mn, mdn, sd, perc_abv_08, perc_abv_06)) %>% 
  select(-c(uid, modobj, corvals)) %>% 
  xlsx::write.xlsx(x=., file=file.path(out_dir, 'top_models_summary.xlsx'),
                   sheetName = "Sheet1",
                   col.names = TRUE, append = FALSE)

# Saving selected model(s) results ----------

# Full Rdata object
topmods %>% 
  filter(Sel.Model) %>% 
  save(., file=file.path(sel_path, 'best_selected.Rdata'))

# Summary of the best model results
topmods %>% 
  filter(Sel.Model) %>% 
  select(-c(mn, mdn, sd, perc_abv_08, perc_abv_06)) %>% 
  select(-c(uid, modobj, corvals)) %>% 
  xlsx::write.xlsx(x=., file=file.path(sel_path, 'selected_best_summary.xlsx'),
                   sheetName = "Sheet1",
                   col.names = TRUE, append = FALSE)

# Plots of the best model(s)
topmods %>% 
  filter(Sel.Model) %>% 
  pwalk(function(modobj, Method, hypothesis, ...){
    fname <- paste0(hypothesis,'.png')
    yrdf <- modobj$Data$Data$mano
    png(file=file.path(sel_path, fname))
    qplot_cdyn(model = modobj$Predictions[[Method]], 
               yrdf = yrdf, 
               Cat.tstep = dim(yrdf)[[1]])
    dev.off()
  })

# Parameters from the best model(s) to a dataframe
topmods %>% 
  filter(Sel.Model) %>% 
  pmap_dfr(function(modobj, Method, hypothesis, uid, ...){
    CatDynPar(modobj, Method) %>% 
      mutate(across(c(Estimates, CVpCent), ~{round(.x, 5)})) %>% 
      mutate(hypothsis = hypothesis, uid = uid)
  }) %>% 
  xlsx::write.xlsx(x=., file=file.path(sel_path, 'selected_best_params.xlsx'),
                 sheetName = "Sheet1",
                 col.names = TRUE, append = FALSE)
