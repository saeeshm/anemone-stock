# Author: Saeesh Mangwani
# Date: 2023-05-26

# Description: Diagnosing results from the full set of models fitted in script
# 3. 
# This script performs all quality filters to only retain models that meet
# minimum standards defined prior to this analysis. This script also visualizes
# parameter distributions and identified candidates weeks for pulse hypothesis
# among the surviving model set.
# This analysis gives a sense of the breadth and distribution of parameter
# values, useful not only for informing the reasonable selection of a best model
# but also to understand the variability in parameter estimates

# ==== Libraries ====
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
source("scripts/_help_funcs.R")
source("scripts/CDMN.functions.R")

# ==== Global/User variables ====

# Current year being modelled
modyear <- 2017

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
table(empty_index)
results <- results_full[empty_index]

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
save(list = saveobjs, file=file.path(result_dir, 'full_model_comparison.Rdata'))

# ==== Filtering models by quality ====

# Temporary result dataframe that will get moved around as we filter
modcomp <- modcomp_main
print(paste('Starting with a total of', nrow(modcomp), 'models'))

# Filtering by maximum absolute gradient ----------

# Removing NA AIC values - models with no convergance, so no diagnostics
modcomp <- modcomp %>% 
  # Removing values where the maximum abs gradient is less than 1 (above this is
  # implausible values)
  filter(Max.Abs.Grads.<1)
print(paste(nrow(modcomp), 'models remaining'))

# Filtering by Convergence code  ----------

# Getting convergence codes for each model
convs <-  pmap_dfr(modcomp, function(uid, Method, modobj,...){
  convcode <- modobj$Model[[Method]]$converg
  tibble('uid' = uid, 'conv_code' = as.character(convcode))
})

# Joining convergance code back to the model dataframe
modcomp <- modcomp %>% 
  left_join(convs, by='uid') %>% 
  # Filtering for only convergance codes 0 or 2 (convergance, or pretty close)
  filter(conv_code %in% c('0', '2'))

print(paste(nrow(modcomp), 'models remaining'))

# Filtering by KKTS conditions ----------
kkts <- pmap_dfr(modcomp, function(uid, Method, modobj,...){
  kkvals <- modobj$Model[[Method]]$kkt
  kktab <- kkvals %>% mutate('uid' = uid)
  rownames(kktab) <- NULL
  return(kktab)
})

# Joining kkts data back to the model dataframe
modcomp <- modcomp %>% 
  left_join(kkts, by='uid') %>%
  # Removing missing kkt models, these we consider failure
  filter(!is.na(kkt1) & !is.na(kkt1)) %>% 
  # Keeping only models where at least 1 condition is true (not super
  # restrictive but oh well)
  filter(kkt1 | kkt2)

print(paste(nrow(modcomp), 'models remaining'))

# Filtering by parameter correlations ----------

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
    # 'skew' = skewness(cvals_cln),
    # 'kurt' = kurtosis(cvals_cln),
    'perc_abv_08' = mean(abs(cvals_cln) > 0.8),
    'perc_abv_06' = mean(abs(cvals_cln) > 0.6),
    'corvals' = list(cvals_cln)
  ) %>% 
    mutate(across(c('mn', 'mdn', 'sd', 'perc_abv_08', 'perc_abv_06'), ~{round(.x, 3)}))
    # mutate(mn_mdn_diff = mn-mdn)
  return(optmdf)
})

# Removing models where more than 50% of parameter correlations are higher than
# 0.6 
modcordf <- modcordf_full %>% 
  # Removing results where median correlation is 1 - that means all parameters
  # were perfectly correlated with one another
  filter(mdn != 1) %>% 
  # Removing models where more than half of the correlations are above 60%
  filter(!perc_abv_06 > 0.5) 

# Filtering models by joining the filtered by correlation table with the full
# table - retains only models that clear the filters
modcomp <- modcordf %>% 
  rename('Method' = optimizer) %>% 
  left_join(modcomp, by=c('hypothesis', 'Method'))

# How many models left?
print(paste(nrow(modcomp), 'models remaining'))

# Saving surviving models in an Rdata file
saveobjs <- ls()[ls() %in% c('modcomp')]
save(list = saveobjs, file=file.path(result_dir, 'surviving_models.Rdata'))

# Exporting surviving model information to an excel sheet
modcomp %>% 
  select(-c(mn, mdn, sd, perc_abv_08, perc_abv_06)) %>% 
  select(-c(uid, modobj, corvals, Sel.Model)) %>% 
  xlsx::write.xlsx(x=., file=file.path(result_dir, 'surviving_models_summary.xlsx'),
                   sheetName = "Sheet1",
                   col.names = TRUE, append = FALSE)

# ==== Diagnosing surviving models ====

# How many mods of each type
mcount <- setNames(as.vector(table(modcomp$Model)), 
                   names(table(modcomp$Model)))
mcount

# Proportional to how many were intially present for each type?
mcount_full <- setNames(as.vector(table(modcomp_main$Model)), 
                        names(table(modcomp_main$Model)))
round(mcount/mcount_full[names(mcount)], 3)

# Visualizing parameter variability ----------

# Getting all the converged parameter values for each model as a dataframe
paramtab <- pmap_dfr(modcomp, function(uid, Method, modobj,...){
  pars <- modobj$Predictions[[Method]]$Model$Parameters
  data.frame(as.list(pars)) %>% 
    mutate(uid = uid, .before=M)
}) %>% tibble()

# Joining paramater data with model attribute information
paramtab <- modcomp %>% 
  select(uid, hypothesis, Model, Distribution, Method) %>% 
  left_join(paramtab, by='uid')

# Function for quick plotting distributions (with an option to filter for
# outliers)
plotDist <- function(pvar, grpvar, paramtab, title, filter_outliers=F){
  if(filter_outliers){
    iqr <-  IQR(paramtab[[pvar]], na.rm=T)
    quants <- quantile(paramtab[[pvar]], probs=c(0.25, 0.75), na.rm=T)
    lb <- quants[1]-(iqr*3)
    ub <- quants[2]+(iqr*3)
    paramtab <- paramtab %>% 
      filter(!!sym(pvar) >= lb & !!sym(pvar) <= ub) 
  }
  paramtab %>%
    ggplot(aes(y=!!sym(pvar), x=!!sym(grpvar), colour=!!sym(grpvar))) +
    geom_boxplot(alpha = 0.5, show.legend=F) +
    geom_jitter(alpha = 0.3, show.legend=F) + 
    # geom_hline(yintercept = median(paramtab$M, na.rm=T), 
    #            linetype='dashed', colour='darkgrey', linewidth=0.8) +
    geom_hline(yintercept = mean(paramtab[[pvar]], na.rm=T), 
               linetype='dashed', colour='firebrick', linewidth=0.8) +
    theme_minimal(11, 'serif') +
    theme(plot.background=element_rect(colour=NA, fill='white')) +
    labs(x=NULL, 
         title=title, 
         subtitle = ifelse(filter_outliers, 'Outliers removed', ''))
}

# Diagnosis plot directory
dplot_dir <- file.path(result_dir, 'param_diagnosis_plots')
if (!dir.exists(dplot_dir)) dir.create(dplot_dir)

# All combinations of parameters, with names that will be used as their titles
params <- c('M', 'N0', 'k', 'alpha', 'beta', 'P1', 'Q1', 'P2', 'Q2') %>% 
  setNames(c('Mortality', 'Initial population', 
             'K-parameter', 'Alpha', 'Beta', 
             'Magnitude of in-pulse 1', 'Magnitude of out-pulse 1',
             'Magnitude of in-pulse 2', 'Magnitude of out-pulse 2'))

# Plotting parameter distributions for all parameters (full and outliers removed)
iwalk(params, ~{
  print(.y)
  # Plotting by Hypothesis type:
  pname <- paste0(.x, '_by_hyp_full.png')
  plotDist(.x, 'Model', paramtab, .y, F) %>% 
    ggsave(file.path(dplot_dir, pname), plot = .)
  
  pname <- paste0(.x, '_by_hyp_outrm.png')
  plotDist(.x, 'Model', paramtab, .y, T) %>% 
    ggsave(file.path(dplot_dir, pname), plot = .)
  
  # Plotting by Optimization method:
  pname <- paste0(.x, '_by_optm_full.png')
  plotDist(.x, 'Method', paramtab, .y, F) %>% 
    ggsave(file.path(dplot_dir, pname), plot = .)
  
  pname <- paste0(.x, '_by_optm_outrm.png')
  plotDist(.x, 'Method', paramtab, .y, T) %>% 
    ggsave(file.path(dplot_dir, pname), plot = .)
})

# Visualizing spread of identified pulse timings ----------

# Getting pulse dates from all models (not-just-surviving) as a dataframe
datestab <- pmap_dfr(modcomp, function(uid, Method, modobj,...){
  dts <- modobj$Predictions[[Method]]$Model$Dates
  names(dts) <- str_remove_all(names(dts), '\\.entry|\\.exit')
  data.frame(as.list(dts)) %>% 
    mutate(uid = uid, .before=ts.start)
})

# Getting attributes from model dataframe
datestab <- modcomp %>% 
  select(uid, hypothesis, Model, Distribution, Method) %>% 
  left_join(datestab %>% select(-ts.start, -ts.end), by='uid')

# Function for visualizing
plotPtime <- function(timevar, grpvar, title){
  plotdf <- datestab %>% 
    mutate(pvar = factor(!!sym(timevar))) %>% 
    filter(!is.na(pvar)) %>% 
    group_by(!!sym(grpvar)) %>% 
    # For each hypothesis type
    group_modify(~{
      # Total number of models in this hypothesis type
      num_mods <- nrow(.x)
      # Calculating the proportion of times each week was detected in this
      # hypothesis, relative to the total number of surviving models
      .x %>% 
        group_by(pvar) %>% 
        summarize('prop_detected' = n()/num_mods)
    }) %>% 
    ungroup()
  
  # Plotting
  plotdf %>% 
    ggplot(aes(x = pvar, y=prop_detected, fill=Model)) +
    geom_col(position='dodge', colour='black', linewidth=0.2, alpha=0.8) +
    theme_minimal(11, 'serif') +
    theme(plot.background=element_rect(colour=NA, fill='white')) + 
    labs(x = 'Candidate week', title=title, y = '% detected')
}

# Visualizing - In-pulse 1
plotPtime('ts.P1', 'Model', 'In-pulse 1') %>% 
  ggsave(file.path(dplot_dir, 'det_pulses_in_1.png'), plot = .)
plotPtime('ts.P2', 'Model', 'In-pulse 2') %>% 
  ggsave(file.path(dplot_dir, 'det_pulses_in_2.png'), plot = .)
plotPtime('ts.Q1', 'Model', 'Out-pulse 1') %>% 
  ggsave(file.path(dplot_dir, 'det_pulses_out_1.png'), plot = .)
plotPtime('ts.Q2', 'Model', 'Out-pulse 2') %>% 
  ggsave(file.path(dplot_dir, 'det_pulses_out_2.png'), plot = .)

# ==== Selecting best models ====
  
# Selecting the best 10 based on AIC (per hypothesis type and distribution) ----------
selmods <- modcomp %>% 
  group_by(Model, Distribution) %>% 
  group_modify(~{
    .x %>% 
      arrange(AIC) %>% 
      head(10)
  }) %>% 
  ungroup()
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

# ==== Selecting and plotting the "best" 2 models from each model type ====
bestmods <- selmods %>% 
  group_by(Model) %>% 
  group_modify(~{
    # Best is defined as ranked based on:
    .x %>% 
      # AIC, lowest average CVP, lowest mean and median parameter correlations
      arrange(AIC, mean_cvp, abs(mn), abs(mdn), abs(sd)) %>%
      # Removing unnnecessary columns
      # select(-c(mn, mdn, sd, skew, kurt, perc_abv_08, perc_abv_06, mn_mdn_diff)) %>% 
      head(2)
  }) %>% 
  ungroup()

# ==== Manually examining best model fits to pick a single best ====

# Getting predicted results from each "best" model
bestpreds <- pmap(bestmods, function(hypothesis, Method, modobj, ...){
  modobj$Predictions[[Method]]
})
names(bestpreds) <- bestmods$hypothesis
names(bestpreds)

# Selecting which model to examine
modnum <- 5
paste("Examining hypothesis:", names(bestpreds)[modnum])

# Plotting model fit
yrdf <- bestmods[modnum,]$modobj[[1]]$Data$Data$mano
qplot_cdyn(model = bestpreds[[modnum]], yrdf = yrdf, Cat.tstep = dim(yrdf)[[1]])

# Checking parameters
CatDynPar(bestmods[modnum,]$modobj[[1]], bestmods[modnum,]$Method) %>% 
  mutate(across(c(Estimates, CVpCent), ~{round(.x, 5)}))

# Adding selection flag to best model dataframe
selindex <- c(5,6)
selids <- bestmods$uid[selindex]
bestmods <- bestmods %>% 
  mutate(Sel.Model = ifelse(uid %in% selids, T, F))

# Saving best model table to an Rdata file
saveobjs <- ls()[ls() %in% c('bestmods')]
save(list = saveobjs, file=file.path(result_dir, 'best_models.Rdata'))

# Saving only summary parameters to a simplified excel sheet
bestmods %>% 
  select(-c(mn, mdn, sd, perc_abv_08, perc_abv_06, mn_mdn_diff)) %>% 
  select(-c(uid, modobj, corvals)) %>% 
  xlsx::write.xlsx(x=., file=file.path(result_dir, 'best_models_summary.xlsx'),
                   sheetName = "Sheet1",
                   col.names = TRUE, append = FALSE)
