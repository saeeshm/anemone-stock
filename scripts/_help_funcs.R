# Author: Saeesh Mangwani
# Date: 2023-05-24

# Description: Helper functions for fitting GDM models with CatDyn

# ==== Libraries ====
library(purrr)
library(gtools)

# ==== General-use functions ====

# A symbol function that returns null if an object doesn't exist, otherwise it
# returns that object
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) {
    lhs
  } else {
    rhs
  }
}

# Quick plotting CatDyn model results
qplot_cdyn <- function(model, yrdf, ...){
  args <- list(...)
  plot(
    x=model,
    leg.pos = args$leg.pos %||% "topleft",
    Biom.tstep = args$Biom.tstep %||% 1,
    Biom.xpos = args$Biom.xpos %||% 1.9,
    Biom.ypos = args$Biom.ypos %||% 0.8,
    Cat.tstep = args$Cat.tstep %||% dim(yrdf)[1],
    Cat.xpos = args$Cat.xpos %||% 1.9,
    Cat.ypos = args$Cat.ypos %||% 0.7
  )
}

# ==== Hypothesis selection functions ====

# Helper function that checks, for a given dataframe of candidate weeks
# identified based on model residuals - 1. if any weeks are in succession & 2.
# if any of these succession weeks are identified by the same set of models. For
# these case, it chooses only the week with the largest absolute residual
clean_succession_weeks <- function(weeks){
  weeks <- arrange(weeks, week)
  # Are any weeks in sucesssion?
  inSucc <- (weeks$week) - lag(weeks$week)
  sindex <- which(inSucc == 1)
  # If none, return the weeks as-is
  if(length(sindex) == 0) {
    return(weeks)
    # Otherwise, separating out distinct weeks from sequential weeks so that we
    # can process the sequential ones separately
  }else{
    dist_weeks <- weeks[!(1:nrow(weeks) %in% unique(c(sindex, sindex-1))),]
    seq_weeks <- weeks[(1:nrow(weeks) %in% unique(c(sindex, sindex-1))),]
  }
  # For the sequential weeks, separately managing duplicates per model set (ie
  # if the same set of models identifies a particular mix of weeks, they should
  # be treated as candidate repeats)
  .x = seq_weeks[1:3,]
  sqn_cln <- seq_weeks %>% 
    group_by(mods) %>% 
    group_modify(~{
      # If there is only week for this model, returning as-is
      if(nrow(.x) == 1) return(.x)
      # Adding a column to hold succession status, and setting true for now
      .x <- .x %>% arrange(week)
      .x$inSucc <- c(NA, rep(T, nrow(.x)-1))
      # While there is at least 1 week in succession
      while((nrow(.x) > 1) & any(.x$inSucc)){
        # Getting the first row where the succession flag is true and comparing
        # with the previous row
        currindex <- c(which(.x$inSucc)[1], c(which(.x$inSucc)[1]-1))
        # Selecting only the week with the highest absolute mean residual
        selweek <- .x[currindex, ] %>% 
          filter(abs(mean_resid) == max(abs(mean_resid)))
        # Removing both compared rows from the table and adding back only the
        # selected one
        .x %>% 
          filter(!(1:nrow(.x) %in% currindex)) %>% 
          bind_rows(selweek) %>% 
          # Rechecking for succession
          arrange(week) %>% 
          mutate(inSucc = (week - lag(week) == 1)) %>% 
          assign('.x', ., inherits=T)
      }
      return(.x)
    })
  # Joining the sequential back to the distinct weeks are returning
  out <- bind_rows(dist_weeks, sqn_cln) %>% arrange(week)
  return(out)
}

# Construct a hypothesis table from a given vector of in/out pulses, and a
# vector of pulse scenarios to construct
make_hyp_table <- function(pos_weeks, neg_weeks, pvec, ts_range, edgethresh=1, filterio=T){
  # Removing weeks within a certain number of timesteps near the edge of the
  # timeseries (defaults to 1 timestep)
  start_edge <- head(ts_range, 1) + edgethresh
  end_edge <- tail(ts_range, 1) - edgethresh
  pos_weeks <- pos_weeks[(pos_weeks > start_edge) & (pos_weeks < end_edge)]
  neg_weeks <- neg_weeks[(neg_weeks > start_edge) & (neg_weeks < end_edge)]
  
  # Constructing table of all possible combinations by calling the
  # hyps-from-pulse function
  out_tab <- map_dfr(pvec, ~{hyps_from_pulse(pos_weeks, neg_weeks, .x, filterio)})
  
  # Returning the full hypothesis table
  return(out_tab)
}

# A private helper function for constructing all hypotheses for a single pulse
# scenario
hyps_from_pulse <- function(pos_weeks, neg_weeks, pulse, filterio=T){
  if(abs(pulse) > length(pos_weeks) | abs(pulse) > length(pos_weeks)){
    msg <- paste('Requested pulses exceeds the number of candidate weeks.', 
                 'Select more weeks or limit the number of possible pulses')
    stop(msg)
  }
  if(pulse == 0){
    out_tab <- tibble('type'=0, in_week=NA_character_, out_week=NA_character_)
  }else if(pulse > 0){
    # Getting all combinations of in-weeks for in-pulses
    comb_tab <- combinations(n=length(pos_weeks), r=pulse, v=pos_weeks) %>% 
      as.data.frame()
    # Converting to dataframe, and adding an empty column for out_weeks
    out_tab <- comb_tab %>% 
      tidyr::unite(col='in_week', everything(), sep = ',') %>% 
      mutate(out_week = NA_character_) %>% 
      mutate(type = pulse, .before='in_week')
  }else{
    # Getting all combinations of in-weeks as a vector
    ins <- combinations(n=length(pos_weeks), r=abs(pulse), v=pos_weeks) %>% 
      as.data.frame() %>% 
      rename_with(~str_replace_all(.x, 'V', 'I'))
    # Same for out-weeks
    outs <- combinations(n=length(neg_weeks), r=abs(pulse), v=neg_weeks) %>% 
      as.data.frame() %>% 
      rename_with(~str_replace_all(.x, 'V', 'O'))
    
    # Joining all ins and outs columns row-wise, which gives all combos of ins
    # to outs
    comb_tab <- map_dfr(1:nrow(ins), ~{
      bind_cols(ins[.x, ,drop=F], outs)
    })
 
    # Filtering conditions - only for more than 1-pulse scenarios
    if((abs(pulse) > 1) & filterio){
      for(i in 2:abs(pulse)){
        filtered <- comb_tab %>% 
          # A subsequent in-pulse can't start before the prior out-pulse?
          filter(!(!!sym(paste0('I', i)) < !!sym(paste0('O', i-1))))
        assign('comb_tab', filtered, inherits=TRUE)
      }
    }
    # Converting to an out-table
    out_tab <- comb_tab %>% 
      tidyr::unite(col='in_week', contains('I'), sep = ',') %>% 
      tidyr::unite(col='out_week', contains('O'), sep = ',') %>% 
      mutate(type = pulse, .before='in_week')
  }
  return(out_tab)
}

# ==== Model fitting functions ====

# Creating a timestep range integrating pulse hypotheses, taking a vector of
# total time steps, a vector of in-pulses and a vector of out_pulses
get_ts_range <- function(tsteps, npulses, in_weeks = NA, out_weeks = NA){
  if(npulses == 0){
    out_ts <- range(tsteps)
  }else if (npulses > 0){
    if(abs(npulses) != length(in_weeks)) {
      stop('The number of in-weeks does not equal the requested number of pulses')
    }
    out_ts <- c(head(tsteps, 1), in_weeks, tail(tsteps, 1))
  }else{
    if((abs(npulses) != length(in_weeks)) | length(in_weeks) != length(out_weeks)) {
      stop('The number of items in the week vectors do not match to the requested number of pulses')
    }
    # Getting pulses in an in-out order, from the provided week vectors
    pulse_stops <- unlist(map2(in_weeks, out_weeks, ~{c(.x, .y)}))
    out_ts <- c(head(tsteps, 1), pulse_stops, tail(tsteps, 1))
  }
  return(out_ts)
}

# Pull out the appropriate parameters in the right order, given a full vector of
# initial parameter and a pulse hypothesis
get_relv_params <- function(params, pulse){
  # Separating non-pulse parameters in the right orderstart
  p1 <- params[c('M', 'N0')]
  p2 <- params[c('k', 'alpha', 'beta')]
  # Returning only these if pure depletion
  if(pulse == 0) return(c(p1, p2))
  # Ensuring params are available for the number of pulses requested
  if (is.na(params[paste0('p_in', abs(pulse))])){
    msg <- paste('In-out abundance estimates for', 
                 pulse, 
                 'pulse hypotheses are not available in this parameter set. ')
    stop(msg)
  }
  # Extracting all pulse parameters for the given number of hypotheses
  pulse_vec <- lapply(1:abs(pulse), function(.x){
    params[paste0(c('p_in', 'p_out'), .x)]
  }) %>% unlist()
  # If pulse number is positive, we are only considering in-pulses. Separating
  # these from the full set
  if(pulse > 0) {
    pulse_vec <- pulse_vec[str_detect(names(pulse_vec), 'in')]
  }
  # Returning the pulse parameters in the right order
  return(c(p1, pulse_vec, p2))
}

