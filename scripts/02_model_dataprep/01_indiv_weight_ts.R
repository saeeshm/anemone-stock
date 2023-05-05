# Author: Saeesh Mangwani
# Date: 2023-05-04

# Description: Building databases for estimating weights over time

# ==== Libraries ====
library(dplyr)
library(readr)
library(ggplot2)
library(lubridate)
library(Runuran)

# ==== Paths and global variable ====

# Rapid theme function for quick plots
theme_soosh <- function(base_size = 11, base_family='serif'){
  theme_minimal(base_size, base_family) +
    theme(plot.background = element_rect(fill = 'white', colour=NULL))
}

# Path to where the output individual weight file will be written
out_path <- 'data/indiv_wt_yd.csv'

# ==== Reading data ====

# Field observations - not individual weights, but aggregated weights and
# numbers over a sampling period
field <- read_csv('data/anemone_field_data.csv')
# Field data observations with individual weights - much more data, but also
# much more restricted in time
weight <- read_csv('data/anemone_weight_data.csv')

# ==== Combining data ====

# For the field data, since it is not individual weights, getting an average
# daily weight from all rows where data are available
field_wts <- field %>% 
  filter(!is.na(weight_g)) %>% 
  mutate(num_indv = ifelse(is.na(num_indv1), num_indv2, num_indv1)) %>% 
  select(date, 'tot_weight' = weight_g, num_indv) %>% 
  mutate(weight_g = tot_weight/num_indv) %>% 
  select(date, weight_g) %>% 
  mutate(data_source = 'field_aggregated')
  
# Joining with the observed weight data to get a weight time series
wtts <- weight %>% 
  select(date, weight_g) %>% 
  mutate(data_source = 'indiv_weighting') %>% 
  bind_rows(field_wts)

# ==== Quick explorations of combined weight data ====

# Plot the total temporal ranges
wtts %>% 
  ggplot(aes(x = date, y = weight_g, colour=data_source)) +
  # A minimum weight threshold (if applicable)?
  geom_hline(aes(yintercept = 10), colour='firebrick', linetype='dashed') +
  geom_point(alpha=0.8) +
  theme_soosh()

# Across the seasonal period
wtts %>% 
  mutate(yd = yday(date)) %>% 
  ggplot(aes(x = yd, y = weight_g, colour=data_source)) +
  geom_point(alpha=0.8) +
  theme_soosh()

# ==== Approach 1 - Using just mean and standard deviation ====

# Filtering for a minimum threshold?
# wtts_cln <- wtts
wtts_cln <- wtts %>% 
  filter(weight_g > 10)

# Getting average and standard deviation across the weight range
mwt <- mean(wtts_cln$weight_g)
sdwt <- sd(wtts_cln$weight_g)

# Setting bound limitations to the weight
lb <- min(wtts_cln$weight_g)
ub <- max(wtts_cln$weight_g)

# Create a vector of re-sampled mean weight for each day of the year and
# introducing random error based on a normal distribution, using this mean and
# standard deviation as the parameters - note that the distribution is truncated
# to the range of the data
set.seed(12)
indiv_wtts <- tibble(
  'yd' = 1:366,
  'weight_g' = urnorm(366, mean=mwt, sd=sdwt, 
                      lb=lb, ub=ub)
)

# Visualizing - no relationship with time or seasonality
indiv_wtts %>% 
  ggplot(aes(x=yd, y=weight_g)) +
  geom_point()

# ==== Writing the weight over time as a CSV ====
write_csv(indiv_wtts, out_path)

