# Author: Saeesh Mangwani
# Date: 2023-05-26

# Description: Diagnosing results from hypothesis model fitting for selecting
# the best models

# ==== Libraries ====
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(moments)
source("scripts/_help_funcs.R")
source("scripts/CDMN.functions.R")

# ==== Global/User variables ====

# Current year being modelled
modyear <- 2017

# Distributions tested
dists <- c('apnormal', 'aplnormal')
dists <- setNames(dists, dists)

# Optimization methods tested
optms <- c("spg", "CG", "Nelder-Mead")

# ==== File paths ====

# Path to base output folder for this year's results
out_dir <- file.path('output/gdm', modyear)

# Path to folder where model results will be saved as Rdata files
mod_dir <- file.path(out_dir, 'models')

# Path to folder where selected model results (from this script will be stored)
result_dir <- file.path(out_dir, 'selection')

# ==== Reading data ====

# Listing all model Rdata objects, and reading their data
fnames <- str_remove_all(list.files(mod_dir), '\\.Rdata')
result_paths <- list.files(mod_dir, full.names=T) %>% setNames(fnames)
results_full <- map(result_paths, ~{
  print(.x)
  return(get(load(file=.x)))
})

# Filtering result lists that are empty (all models failed to converge)
empty_index <- (map_int(results_full, length) != 0)
results <- results_full[empty_index]

# ==== Helper functions ====

# Abbreviating distribution names
abb_dists <- function(dnames){
  case_when(dnames == 'apnormal' ~ 'apn')
}

# ==== Preparing the datalist for use with CatDynExp ====
modeldf <- imap_dfr(results, ~{
  # Getting the hypothesis name for this model set
  hyp <- .y
  dists <- .x
  # Renaming model names, with scenario and abbreviated distribution
  names(dists) <- paste0(hyp, '-', substr(names(dists), 1, 3))
  # Checking how many successful optimizations methods were present for each
  # distribution
  optm_names <- map(dists, ~{names(.x$Predictions)})
  # Removing the Predictions object - we don't need it anymore
  # dists <- map(dists, ~{.x$Predictions <- NULL; return(.x)})
  # Converting the list of optimization names per model to a dataframe -
  # helps ensure names and models are correctly aligned
  optmdf <- imap_dfr(optm_names, ~{
    tibble('hypothesis' = .y, 'optimizer' = .x)
  })
  
  # Joining the CatDyn model to this dataframe, using the model name (this
  # ensures the same model is repeated for every optimization method that
  # converged for it)
  optmdf <- optmdf %>% 
    mutate(model = dists[hypothesis])
})

# ==== Running model comparison via CatDyn ====

# Separating hypothesized and pure models since the function cannot diagnose
# them together
puremods <- modeldf %>% filter(str_detect(hypothesis, 'pure'))
imods <- modeldf %>% filter(str_detect(modeldf$hypothesis, '_\\d{1}p(?!\\d{1}p)_'))
iomods <- modeldf %>% filter(str_detect(modeldf$hypothesis, '_\\d{1}p\\d{1}p_'))

# Getting model diagnostics
purecomp <- CatDynSum(x=puremods$model, season="2017", method=puremods$optimizer) %>% 
  mutate(hypothesis = puremods$hypothesis, .before="Fleet")
icomp <- CatDynSum(x=imods$model, season="2017", method=imods$optimizer) %>% 
  mutate(hypothesis = imods$hypothesis, .before="Fleet")
iocomp <- CatDynSum(x=iomods$model, season="2017", method=iomods$optimizer) %>% 
  mutate(hypothesis = iomods$hypothesis, .before="Fleet")

# Cleaning names
purecomp <- purecomp %>% select(-contains('mano'))
names(icomp) <- str_remove_all(names(icomp), '\\.mano')

# joining all diagnostics results into a single dataframe
modcomp_main <- purecomp %>% 
  bind_rows(icomp) %>% 
  bind_rows(iocomp) %>% 
  # Adding a uid column
  mutate(uid = 1:nrow(.), .before='hypothesis') %>% 
  as_tibble() %>% 
  # Adding all the model results as columns, so we can access them
  left_join(modeldf %>% 
              rename('Method' = optimizer, 'modobj' = model), 
            by = c('hypothesis', "Method"))

# Saving the prepared data list and the model comparison object to an Rdata file
saveobjs <- ls()[ls() %in% c('modeldf', 'modcomp_main')]
save(list = saveobjs, file=file.path(result_dir, 'model_comparison.Rdata'))

# ==== Model selection based on comparison results ====

# Temporary objects that will get moved around as we filter
modcomp <- modcomp_main
print(paste('Starting with a total of', nrow(modcomp), 'models'))

# Initial filtering ----------
modcomp <- modcomp %>% 
  # Removing NA AIC values - models with no convergance, so no diagnostics
  filter(!is.na(AIC)) %>% 
  # Removing values where the maximum abs gradient is less than 1 (above this is
  # implausible values)
  filter(Max.Abs.Grads.<1)

print(paste(nrow(modcomp), 'models remaining'))

# Ranking models by correlation structure ----------

# Extracting data for obtaining correlation structure
modcors <- pmap(modcomp, function(hypothesis, Method, modobj, ...){
  cors <- modobj$Model[[Method]]$Cor
  return(list(
    'hyp' = hypothesis,
    'optm' = Method,
    'cors' = cors
  ))
})

# Calculating correlation summary statistics
modcordf_full <- map_dfr(modcors, ~{
  hypname <- .x$hyp
  optm <- .x$optm
  # Getting correlation values
  corvals <- round(.x$cors, 3)
  # Converting to a symmetrical matrix
  cormx <- matrix(corvals, nrow=sqrt(length(corvals)), 
                  ncol=sqrt(length(corvals)), 
                  byrow=T)
  # Extracting upper and lower triangular vals - this removes self-correlation
  # values, and only keeps the ones we actually want to visualize
  cvals_cln <- c(cormx[upper.tri(cormx)], cormx[lower.tri(cormx)])
  # Calculating statistics
  optmdf <- tibble(
    'hypothesis' = hypname,
    'optimizer' = optm,
    'mn' = mean(cvals_cln),
    'mdn' = ifelse(all(abs(cvals_cln) == 1), 1, median(cvals_cln)),
    'sd' = sd(cvals_cln),
    'skew' = skewness(cvals_cln),
    'kurt' = kurtosis(cvals_cln),
    'perc_abv_08' = mean(abs(cvals_cln) > 0.8),
    'perc_abv_06' = mean(abs(cvals_cln) > 0.6),
    'corvals' = list(cvals_cln)
  ) %>% 
    mutate(across(contains('perc'), ~{round(.x, 3)})) %>% 
    mutate(across(c('mn', 'mdn', 'sd', 'skew', 'kurt'), ~{round(.x, 3)})) %>% 
    mutate(mn_mdn_diff = mn-mdn)
  return(optmdf)
})

modcordf_full %>% filter(hypothesis == 'gdm_1p_1-apl')
test <- results[['gdm_1p_1']]['aplnormal']
CatDynCor(test, ttl=rep('test', 1), method='spg', arr=c(1,1))
hist(test$aplnormal$Model$spg$Cor, breaks=10)

# Initial correlation-structure filtering
modcordf <- modcordf_full %>% 
  # Removing results where median correlation is 1 - that means all parameters
  # were perfectly correlated with one another
  filter(mdn != 1) %>% 
  # Removing models where more than half of the correlations are above 60%
  filter(!perc_abv_06 > 0.5) 
  # Ordering based on the proportions of correlations higher than 0.8 (we want
  # this to be low), and then based on the absolutel mean, sd, kurtosis, and
  # skew (we want all to be low)
  # arrange(perc_abv_06, perc_abv_08, abs(mn), abs(sd), abs(kurt), abs(skew)) %>% 
  # Creating a rank variable
  # mutate(corrank = 1:nrow(.))

# Joining correlation ranks back to the original dataframe - this filters models
# removed during this step
modcomp <- modcordf %>% 
  rename('Method' = optimizer) %>% 
  left_join(modcomp, by=c('hypothesis', 'Method'))

# How many models left?
print(paste(nrow(modcomp), 'models remaining'))
# Which kinds of pulse hypotheses have persisted?
table(modcomp$Model)
  
# Ranking models by AIC (per distribution) ----------
modcomp <- modcomp %>% 
  group_by(Model, Distribution) %>% 
  group_modify(~{
    .x %>% 
      arrange(AIC) %>% 
      mutate(aicrank = 1:nrow(.))
  }) %>% 
  ungroup()

# ==== Selecting the "best" 5 models per method + distribution ====

# Selection currently based on an even weighting between AIC and correlation
# rank
selmods <- modcomp %>% 
  group_by(Model, Distribution) %>% 
  group_modify(~{
    .x %>% 
      arrange(aicrank) %>% 
      head(5)
  })
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

# ==== Selecting and plotting the "best" model from each model type ====
bestmods <- selmods %>% 
  group_by(Model) %>% 
  group_modify(~{
    .x %>% 
      # Best is defined as ranked based on lowest average CVP
      arrange(mean_cvp) %>% 
      # And lowest values on indicators for coorelation structure
      arrange(perc_abv_06, abs(mn), abs(mdn), abs(sd)) %>% 
      arrange(AIC) %>% 
      select(-c(mn, mdn, sd, skew, kurt, perc_abv_08, perc_abv_06, mn_mdn_diff)) %>% 
      head(2)
  })

bestpredlist <- pmap(bestmods, function(hypothesis, Method, modobj, ...){
  modobj$Predictions[[Method]]
})
names(bestpredlist) <- bestmods$Model
names(bestpredlist)

# Plotting - visually making a selection on the best model
qplot_cdyn(bestpredlist[[8]], 
           bestmods[1,]$modobj[[1]]$Data$Data$mano, 
           Cat.tstep=dim(bestmods[1,]$modobj[[1]]$Data$Data$mano)[[1]]
)

# Checking parameters
CatDynPar(bestmods[5,]$modobj[[1]], bestmods[5,]$Method)

# ==== General distribution of parameter values ====
# test <- modcomp
# modcomp <- test

# Getting a list of all parameter values for each model in the selected set
paramlist <- pmap(selmods, function(hypothesis, Method, modobj,...){
  CatDynPar(modobj, Method)
})
paramlist <- setNames(paramlist, (selmods$uid))

# Getting a long-form single dataframe giving only Estimates for each param
pardf <- imap_dfr(paramlist, ~{
  .x %>% 
    select(-Timing, -CVpCent) %>% 
    mutate(Parameter = str_replace_all(Parameter, '\\.1\\/', '1_per_')) %>% 
    mutate(Parameter = str_replace_all(Parameter, '\\.', '_')) %>% 
    pivot_wider(names_from = 'Parameter', values_from = 'Estimates') %>% 
    mutate(uid = .y, .before='M1_per_week')
})

# Joining hypothesis information
pardf <- modcomp %>% 
  select(uid, hypothesis, Method) %>% 
  left_join(pardf %>% 
              mutate(uid = as.integer(uid)), by='uid')

# Checking distributions
hist(pardf$M1_per_week)
hist(pardf$N0_thou)
hist(pardf$Recruitment_thou_Wave1)
hist(pardf$Recruitment_thou_Wave2)
hist(pardf$Spawning_thou_Wave1)
hist(pardf$Spawning_thou_Wave2)
hist(pardf$k1_per_nhours)
hist(pardf$alpha)
hist(pardf$beta)


