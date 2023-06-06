# Author: Saeesh Mangwani
# Date: 2023-05-12

# Description: Building a cleaned daily-timestep dataset for anemone catches and
# effort (used as the input for modelling)

# ==== Libraries ====
library(dplyr)
library(readr)
library(lubridate)

# ==== Paths and global variables ====

# Path to catch and effort data
ceff_path <- 'data/raw/catch_effort_infilled.csv'

# Path to market sale data
sales_path <- 'data/raw/fish_market_data.csv'

# ==== Reading data ====
ceff <- read_csv(ceff_path)
sales <- read_csv(sales_path)

# ==== Calculating and adding fishers per day effort ====

#  How many fishers were active per day, and were any names repeated twice (a
#  data error)
fishers_pday <- ceff %>% 
  group_by(date) %>% 
  group_modify(~{
    fisher_table <- table(.x$fisher_code)
    return(tibble(
      'num_active_fishers' = length(fisher_table),
      'repeat_fishers' = any(fisher_table > 1)
    ))
  }) %>% 
  ungroup()

# Adding number of active fishers per day to the catch/effort
ceff <- ceff %>% 
  left_join(fishers_pday %>% 
              select(-repeat_fishers),
            by = 'date')

# ==== Completing the TS for catch, effort and sale data ====

# Getting a complete time series across the period for which data are available
daily_timestep <- tibble('date' = seq(min(ceff$date), max(ceff$date), by = "day"))

# Filling in the daily time series for ceff and sales data
ceff_complete <- daily_timestep %>% 
  left_join(ceff, by='date', multiple='all') %>% 
  # Replacing all NAs with 0s, to indicate no catch/effort/sales
  mutate(across(c(catch_kg, eff_hours_intp, num_active_fishers), ~{
    ifelse(is.na(.x), 0, .x)
  }))

sales_complete <- daily_timestep %>% 
  left_join(sales, by='date', multiple='all') %>% 
  # Replacing all NAs with 0s, to indicate no catch/effort/sales
  mutate(across(c(kg, price_eur_kg, price_tot_eur), ~{
    ifelse(is.na(.x), 0, .x)
  }))

# ==== Preparing a summarized daily dataset ====

# Catch-effort daily
ceff_daily <- ceff_complete %>% 
  group_by(date) %>% 
  group_modify(~{
    fisher_table <- table(.x$fisher_code)
    return(tibble(
      'catch_kg_tot' = sum(.x$catch_kg),
      'eff_h_tot' = sum(.x$eff_hours_intp),
      'eff_num_fisher_tot' = length(fisher_table)
    ))
  }) %>% 
  ungroup()

ceff_daily <- ceff_complete %>% 
  group_by(date) %>% 
  group_modify(~{
    fisher_table <- table(.x$fisher_code)
    return(tibble(
      'catch_kg_tot' = sum(.x$catch_kg),
      'eff_h_tot' = sum(.x$eff_hours_intp),
      'eff_num_fisher_tot' = length(fisher_table)
    ))
  }) %>% 
  ungroup()

# Summarizing catch and effort per week:
sales_daily <- sales_complete %>% 
  # filter(!isintp) %>% 
  group_by(date) %>% 
  summarize(
    date = min(date),
    sales_kg_tot = sum(kg),
    mean_price_pkg = mean(price_eur_kg),
    revenue_tot_eur = sum(price_tot_eur)
  ) %>% 
  ungroup()

# Summarizing and joining to a single daily dataframe
daily <- left_join(ceff_daily, sales_daily, by='date')

# ==== Preparing a summarized weekly dataset ====

# Summarizing catch and effort per week:
ceff_weekly <- ceff_complete %>% 
  # Separating data by week of the year
  mutate(y = year(date), .after=date) %>% 
  mutate(w = week(date), .after=date) %>% 
  # filter(!isintp) %>% 
  group_by(y, w) %>% 
  group_modify(~{
    fisher_table <- table(.x$fisher_code)
    return(tibble(
      'date' = min(.x$date),
      'catch_kg_tot' = sum(.x$catch_kg),
      'eff_h_tot' = sum(.x$eff_hours_intp),
      'eff_num_fisher_tot' = length(fisher_table)
    ))
  }) %>% 
  ungroup() %>% 
  select(-y, -w)

# Summarizing sales per week
sales_weekly <- sales_complete %>% 
  # Separating data by week of the year
  mutate(y = year(date), .after=date) %>% 
  mutate(w = week(date), .after=date) %>% 
  # filter(!isintp) %>% 
  group_by(y, w) %>% 
  summarize(
    date = min(date),
    sales_kg_tot = sum(kg),
    mean_price_pkg = mean(price_eur_kg),
    revenue_tot_eur = sum(price_tot_eur)
  ) %>% 
  ungroup() %>% 
  select(-y, -w)

# Joining weekly sale and ceff data
weekly <- left_join(ceff_weekly, sales_weekly, by='date')

# ==== Writing complete TS datasets to disk ====
write_csv(daily, 'data/daily_ceff_sale_ts.csv')
write_csv(weekly, 'data/weekly_ceff_sale_ts.csv')







