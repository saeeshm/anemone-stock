# Author: Saeesh Mangwani
# Date: 2023-05-12

# Description: Fitting a pure depletion model to the 2017 dataset and using it
# to construct potential hypothesis for pulses

# ==== Libraries ====
library(dplyr)
library(readr)
library(xlsx)
library(tidyr)
library(lubridate)
library(purrr)
library(stringr)
library(optimx)
library(ggplot2)
# library(CatDyn)
source('scripts/_help_funcs.R')
source("scripts/CDMN.functions.r")

# ==== Global variables ====

# Current year being modelled
modyear <- 2017

# Distributions to test
dists <- c('apnormal', 'aplnormal')
dists <- setNames(dists, dists)

# Optimization methods to test
optms <- c("spg", "CG", "Nelder-Mead")

# ==== File paths ====

# Path to output folder for this year's results
out_dir <- file.path('output/gdm', modyear)
# Creating a folder for storing model result RData files
dir.create(file.path(out_dir, 'models'))

# Path to prepared CatDyn time series dataset (from script 011_*)
yrdf_path <- list.files(out_dir, pattern = 'yrdf', full.names = T)

# Path to dataframe containing initialization parameters
param_path <- list.files(out_dir, pattern = 'param', full.names = T)

# Path to Rdata file where the initialized model objects are stored (exploratory
# null model, and the base CatDyn dataframe)
base_rdpath <- file.path(out_dir, 'gdm_base.Rdata')

# Path to Rdata file where fitted pure GDM results will be saved
out_rdpath <- file.path(out_dir, 'models', 'gdm_pure.Rdata')

# Creating a sub-folder for pure_gdm plots
pure_gdm_plot_dir <- file.path(out_dir, 'plots', 'gdm')
dir.create(pure_gdm_plot_dir)

# Creating a folder for storing visualization results from the hypothesis
# detection process
hyps_out_dir <- file.path(out_dir, 'hypothesis_detection')
dir.create(hyps_out_dir)

# Path to the output hypothesis table
hyp_tab_out_path <- file.path(out_dir, paste0('gdm_hypotheses_', modyear, '.xlsx'))

# ==== Reading data ====

# Catch/effort dataframes
yrdf <- read_csv(yrdf_path)
paramdf <- read_csv(param_path)

# Initial model objects - loads the base CatDyn object dataset (anm_cdyn) and
# the exploratory model fit/null model (null_model_exp)
load(base_rdpath)

# ==== Defining initial parameters ====

# Setting season dates for the year
season_range <- c(paste0(modyear, '-01-02'), paste0(modyear, '-12-30'))

# Timestep range:
ts_range <- range(anm_cdyn$Data[[1]]$time.step)

# Getting tuned initial parameters from the dataframe
params <- setNames(paramdf$value, paramdf$param)

# All parameters must be in the log scale for use in the model
log_params <- log(params)
  
# ==== Fitting the ludic model (pure depletion - no perturbations) ====

# If the setup has already been run, reading the results
# file.remove(out_rdpath)
if(file.exists(out_rdpath)) load(out_rdpath)
# Running the models only if the results don't already exist
if(!exists('models')){
  # Fitting using all optimization methods and distributions defined above
  fits <- imap(dists, ~{
    # Printing the name of the distribution
    print(.y)
    # Fitting the model - this can take a WHILE
    base_model_fit <- CatDynFit(
      # CatDyn dataset object - read from the RDS file, defined in script 1
      x=anm_cdyn,
      p=0,
      par=log_params,
      dates=ts_range,
      distr=.x,
      method=optms,
      itnmax=50000
    )
  })
  
  # Getting predictions
  preds <- imap(fit_models, ~{
    print(.y)
    model <- .x
    # Getting predictions from each optimization method
    outpreds <- map(optms, ~{
      print(paste(' ---', .x))
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
  
  # Plotting prediction results and saving to disk
  iwalk(preds, ~{
    dist_name <- .y
    fits <- .x
    iwalk(fits, ~{
      # Creating filename for this plot
      fname <- paste0('gdm_pure_',dist_name, '_', .y, '.png')
      print(fname)
      # Creating plot and saving to disk
      png(file=file.path(pure_gdm_plot_dir, fname))
      qplot_cdyn(.x, yrdf, Cat.tstep=dim(yrdf)[1])
      dev.off()
    })
  })
  
  # Combining the models with their predictions
  models <- imap(fit_models, ~{
    .x$Predictions <- preds[[.y]]
    return(.x)
  })
  
  # Writing to Rdata file
  save(list='models', file=out_rdpath)
}else{
  print('Model results already exist! Reading them from the Rdata file.')
}

# ==== Identifying candidate hypotheses from the Catch-Spike ====

# Getting the CatDyn configured dataset, which contains the calculated
# catch-spike (we have only 1 fleet)
anmdat <- anm_cdyn$Data[[1]]

# Positive perturbation candidates ----------

# From the catchspike, identifying the 10 time-points with the highest values
val_index <- anmdat$spikecat %in% head(sort(anmdat$spikecat, decreasing = T), 10)

# Then the 10 values with the highest positive residuals based on a GAM
# m <- mgcv::gam(spikecat ~ s(time.step, bs = "cr"), data=anmdat)
# resid_index <- m$residuals %in% head(sort(m$residuals, decreasing = T), 10)

# Selecting the ones that meet both conditions are are positive (i.e > 0)
# pos_pertb_vals <- anmdat$spikecat[(anmdat$spikecat > 0) & (val_index & resid_index)]
pos_pertb_vals <- anmdat$spikecat[(anmdat$spikecat > 0) & (val_index)]

# Adding these selections to the dataframe
anmdat <- anmdat %>% 
  mutate(pos_pertb_cspike = spikecat %in% pos_pertb_vals)

# Plotting the catch-spike, with potential positive perturbations identified
pp_cspike <- anmdat %>% 
  ggplot(aes(x = time.step, y = spikecat)) + 
  geom_hline(aes(yintercept=0), linetype='dashed', colour='black') +
  geom_point(aes(colour=pos_pertb_cspike), show.legend = F) +
  scale_colour_manual(values=c('darkgrey', 'darkorange')) +
  # geom_smooth(method = "loess", span = 0.2, colour='firebrick') +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr"),
  colour='firebrick', alpha=0.2) +
  theme_minimal(11, 'serif') +
  theme(plot.background = element_rect(colour=NA, fill='white')) +
  labs(x = 'Week number', y = 'Catch Spike')
pp_cspike

# Saving to disk
ggsave(file.path(hyps_out_dir, 'pos_pertb_cspike.png'), 
       plot = pp_cspike, 
       width=7, height=5, dpi=300)

# Negative perturbation candidates ----------

# Examining residuals - here we are interested in the weeks with the biggest
# DROP relative to a previous week. So first computing a lag difference of the
# catch-spike
anmdat <- anmdat %>% 
  mutate(lagdiff = spikecat - lag(spikecat)) %>% 
  # Getting the lag difference from the previous time-step for each week as well
  mutate(lagdiff2 = lag(lagdiff)) %>% 
  rowwise() %>% 
  # Averaging the lag difference for the current and previous time-step - allows
  # us to prioritize weeks where there is a consistent drop over at least 2
  # weeks
  mutate(meanlag = mean(c(lagdiff, lagdiff2), na.rm=F)) %>% 
  ungroup() %>% 
  # Identifying the 10 points with the largest negative mean lag
  mutate(neg_pertb_cspike = meanlag %in% head(sort(.$meanlag), 10))

# If two identified weeks are in succession, selecting only the latter (this is
# the one which is more negative)
nnwks <- filter(anmdat, neg_pertb_cspike)$time.step
nnwks[which((nnwks - lag(nnwks)) == 1) - 1] <- NA
nnwks_cln <- na.omit(nnwks)

# Editing the negative pertb boolean variable to only indicate these filtered
# weeks
anmdat <- anmdat %>% 
  mutate(neg_pertb_cspike = ifelse(time.step %in% nnwks_cln, T, F))

# Plotting these lowest lag points (these are canditates for emigration spikes)
np_cspike <- anmdat %>% 
  ggplot() +
  geom_hline(aes(yintercept=0), linetype='dashed', colour='black') +
  geom_line(aes(x = time.step, y=meanlag), colour = 'firebrick', 
            alpha = 0.5) +
  geom_line(aes(x = time.step, y=spikecat), colour = 'darkgrey') +
  geom_point(aes(x = time.step, y=spikecat, colour=neg_pertb_cspike), show.legend=F) +
  scale_colour_manual(values=c('NA', 'darkorange')) +
  # geom_point(aes(x = time.step, y=meanlag, colour=neg_pertb_cspike), 
  #            show.legend = F) +
  # geom_point(aes(x = time.step, y=lagdiff), colour = 'darkgreen') +
  # geom_line(aes(x = time.step, y=lagdiff), colour = 'darkgreen') +
  # geom_point(aes(x = time.step, y=lagdiff2), colour = 'orange') +
  # geom_line(aes(x = time.step, y=lagdiff2), colour = 'orange') +
  # geom_point(aes(x = time.step, y=lagdiff3), colour = 'darkgreen') +
  # geom_line(aes(x = time.step, y=lagdiff3), colour = 'darkgreen') +
  theme_minimal(11, 'serif') +
  theme(plot.background = element_rect(colour=NA, fill='white')) +
  labs(x = 'Week number', y = 'Catch Spike')
np_cspike

# Saving to disk
ggsave(file.path(hyps_out_dir, 'neg_pertb_cspike.png'), 
       plot = np_cspike, 
       width=7, height=5, dpi=300)

# ==== Identifying candidate hypotheses from model residuals ====

# Getting deviance residuals from each model as a dataframe
drs <- imap_dfr(models, ~{
  distname <- .y
  optms <- .x$Predictions
  resid_df <- imap_dfr(optms, ~{
    resids <- .x$Model$Results$Deviance.Residuals
    weeks <- .x$Model$Results$Period.week
    iqr <- IQR(resids)
    # Using a standard outlier threshold of 1.5 times the IQR
    outlier_thresh <- 1.5
    # Identifying potential outliers based on the simple threshold*IQR formulat
    tibble(
      'optm' = .y,
      'week' = weeks,
      'residuals' = resids,
      'iqr' = iqr,
      'isextreme' = (resids > outlier_thresh*iqr) | (resids < -outlier_thresh*iqr)
    )
  })
  resid_df %>% 
    mutate(dist = distname, .before=optm)
})

# Plotting residuals distributions and identified outliers
resid_outliers <- drs %>% 
  # filter(dist=='apnormal') %>% 
  ggplot(aes(x = week, y = residuals)) +
  geom_boxplot(colour = 'black', fill=NA, outlier.shape = NA) +
  geom_jitter(aes(colour = isextreme, group=isextreme), alpha=0.8, width=0.2) +
  scale_colour_manual(values=c('darkgrey', 'darkorange')) +
  theme_minimal(11, 'serif') +
  theme(plot.background = element_rect(colour=NA, fill='white')) +
  labs(y = 'Deviance residual', x = 'Week', colour='Outlier') +
  guides(colour='none') +
  facet_wrap(c('dist','optm'), nrow=2, scales='free_y')
resid_outliers

# Saving to disk
ggsave(file.path(hyps_out_dir, 'dev_resid_outliers.png'), 
       plot = resid_outliers, 
       width=7, height=5, dpi=300)

# How many unique models are being examined?
num_models <- drs %>% 
  tidyr::unite(modname, dist, optm) %>% 
  pull(modname) %>% 
  unique() %>% 
  length()

# Ranking identified positive and negative outlier weeks by the number of times
# they are identified by the models. For selecting "significant" hypotheses, we
# only select weeks identified by at least half of the models
resid_candidates <- drs %>% 
  mutate(dist = str_replace_all(dist, 'normal', 'n')) %>% 
  tidyr::unite(mod, dist, optm, sep='-', remove = F) %>% 
  filter(isextreme) %>% 
  group_by(week) %>% 
  summarize('identified' = sum(isextreme),
            'mean_resid' = mean(residuals),
            'mods' = paste0(mod, collapse = ',')) %>% 
  arrange(desc(identified)) %>% 
  # Filtering only weeks identified by at least half of all models
  filter(identified >= num_models/2)

# Separating positive and negative perturbations
pos_weeks <- resid_candidates %>% filter(mean_resid > 0)
neg_weeks <- resid_candidates %>% filter(mean_resid < 0)

# Calling the helper function to clean the identified hypothesis weeks for
# succession, and adding an indicator column to select these weeks in the
# original dataframe
anmdat <- anmdat %>% 
  mutate(pos_pertb_resid = time.step %in% clean_succession_weeks(pos_weeks)$week) %>% 
  mutate(neg_pertb_resid = time.step %in% clean_succession_weeks(neg_weeks)$week)

# ==== Identifying hypotheses from perturbation analysis ====
pos_weeks <- anmdat %>% 
  filter(pos_pertb_cspike | pos_pertb_resid) %>% 
  pull(time.step) %>% 
  unique()
neg_weeks <- anmdat %>% 
  filter(neg_pertb_cspike | neg_pertb_resid) %>% 
  pull(time.step) %>% 
  unique()

# Which ones are identified by both methods - potentially superior candidates
pos_weeks_strong <- anmdat %>% 
  filter(pos_pertb_cspike & pos_pertb_resid) %>% 
  pull(time.step) %>% 
  unique()
neg_weeks_strong <- anmdat %>% 
  filter(neg_pertb_cspike & neg_pertb_resid) %>% 
  pull(time.step) %>% 
  unique()

# ==== Constructing the hypothesis table from candidates ====

# Calling the helper function to construct the hypothesis table from these weeks
pulses <- c(0, 1, 2, -1, -2)
hyp_tab <- make_hyp_table(
  pos_weeks, 
  neg_weeks, 
  pvec=pulses, 
  ts_range = ts_range,
  filterio=F)

# Identifying which ones may be strong hypotheses - those where all candidates
# consist of the strong choices
hyp_tab <- hyp_tab %>% 
  rowwise() %>% 
  mutate(isStrong = all(as.numeric(str_split(in_week, ',')[[1]]) %in% pos_weeks_strong) & 
           all(as.numeric(str_split(out_week, ',')[[1]]) %in% neg_weeks_strong))

# Writing the hypothesis table to disk
xlsx::write.xlsx(hyp_tab, file=hyp_tab_out_path, append=F)
