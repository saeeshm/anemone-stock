# Author: Saeesh Mangwani
# Date: 2023-05-03

# Description: Presentation plots from an exploratory analysis characterizing
# the octopus fishery

# ==== Libraries ====
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(lubridate)
library(patchwork)

# ==== Paths and global variable ====

# Path to output directory for plots
out_dir <- 'output/eda-plots'

# Rapid theme function for quick plots
theme_soosh <- function(base_size = 11, base_family='serif'){
  theme_minimal(base_size, base_family) +
    theme(plot.background = element_rect(fill = 'white', colour=NA))
}

# ==== Reading data ====
ceff <- read_csv('data/catch_effort_infilled.csv') %>% 
  # Column indicating if data was interpolated or not
  mutate(isintp = is.na(eff_hours))
field <- read_csv('data/anemone_field_data.csv')
weight <- read_csv('data/anemone_weight_data.csv')
sales <- read_csv('data/fish_market_data.csv')

# ==== Adding relevant variables and summarizing to a weekly timestep ====

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

# Filling in the daily time series for ceff and sales data
daily <- tibble('date' = seq(min(ceff$date), max(ceff$date), by = "day"))
ceff_daily <- daily %>% 
  left_join(ceff, by='date', multiple='all') %>% 
  # Replacing all NAs with 0s, to indicate no catch/effort/sales
  mutate(across(c(catch_kg, eff_hours_intp, num_active_fishers), ~{
    ifelse(is.na(.x), 0, .x)
  }))

sales_daily <- daily %>% 
  left_join(sales, by='date', multiple='all') %>% 
  # Replacing all NAs with 0s, to indicate no catch/effort/sales
  mutate(across(c(kg, price_eur_kg, price_tot_eur), ~{
    ifelse(is.na(.x), 0, .x)
  }))

# Summarizing catch and effort per week:
ceff_weekly <- ceff_daily %>% 
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
sales_weekly <- sales_daily %>% 
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

# ==== Time series plots ====

# Getting the year-range of the data
yrange <- year(range(weekly$date))
yseq <- yrange[1]:yrange[2]

# Vector with dates of year changes
year_change <- ymd(paste0(yseq[-1], '-01-01'))

# Vector with dates of closed seasons
closed_season <- ymd(c(paste0(yseq[3:length(yseq)], '-04-01'), 
                       paste0(yseq[3:length(yseq)], '-05-31'))) %>% 
  sort()

# Catch, effort and CPUE over time ----------

# Catch
p1 <- weekly %>% 
  ggplot() +
  # geom_vline(xintercept=year_change, colour='darkgrey') +
  geom_line(aes(y = catch_kg_tot, x = date), colour='darkgrey', alpha=0.9) +
  geom_point(aes(y = catch_kg_tot, x = date), colour='black', stroke=0.3,
             alpha=0.2, size=0.3) +
  geom_vline(xintercept = closed_season, colour='red', linetype='dashed') +
  theme_soosh() +
  theme(axis.text.x = element_blank()) +
  labs(
    x = NULL,
    y = 'Catch (kg)'
  )

# Effort - hours
p2 <- weekly %>% 
  ggplot() +
  # geom_vline(xintercept=year_change, colour='darkgrey') +
  geom_line(aes(y = eff_h_tot, x = date), colour='darkgrey', alpha=0.9) +
  geom_point(aes(y = eff_h_tot, x = date), colour='black', stroke=0.3,
             alpha=0.2, size=0.3) +
  geom_vline(xintercept = closed_season, colour='red', linetype='dashed') +
  theme_soosh() +
  labs(
    x = NULL,
    y = 'Effort (hours)'
  )

# CPUE
p3 <- weekly %>% 
  mutate(cpue = catch_kg_tot/eff_h_tot) %>% 
  ggplot() +
  # geom_vline(xintercept=year_change, colour='darkgrey') +
  geom_line(aes(y = cpue, x = date), colour='darkgrey', alpha=0.9) +
  geom_point(aes(y = cpue, x = date), colour='black', stroke=0.3,
             alpha=0.2, size=0.3) +
  geom_vline(xintercept = closed_season, colour='red', linetype='dashed') +
  theme_soosh() +
  labs(
    x = NULL,
    y = 'CPUE'
  )

# Combined plot
mosaic <- p1/p2/p3

# Writing to disk
ggsave(file.path(out_dir, 'cat_eff_cpue_ts.png'), 
       mosaic, 
       width=7, height=5, dpi=300)

# Relationships between each type of effort and catch ----------

# Catch vs hours
p1 <- weekly %>% 
  ggplot(aes(y = catch_kg_tot, x = eff_h_tot)) +
  geom_point(colour='darkgrey', alpha=0.8) +
  geom_smooth(method = 'lm', colour = 'firebrick', alpha=0.6) +
  theme_soosh() +
  ylim(-10, 1250) +
  labs(
    y = 'Catch (kg)',
    x = 'Effort (hours)'
  )

# Catch vs number of fishers
p2 <- weekly %>% 
  ggplot(aes(y = catch_kg_tot, x = eff_num_fisher_tot)) +
  geom_point(colour='darkgrey', alpha=0.8) +
  geom_smooth(method = 'lm', colour = 'firebrick', alpha=0.6) +
  theme_soosh() +
  ylim(-10, 1250) +
  labs(
    y = NULL,
    x = 'Effort (Number of fishers)'
  )

# Saving as mosaic
mosaic <- p1+p2
ggsave(file.path(out_dir, 'cat_eff_smooth.png'), 
       mosaic, 
       width=7, height=5, dpi=300)
