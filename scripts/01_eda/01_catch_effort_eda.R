# Author: Saeesh Mangwani
# Date: 2023-04-28

# Description: Exploratory analysis of fishing catch, effort and weight data

# ==== Libraries ====
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(plotly)

# ==== Paths and global variable ====

# Rapid theme function for quick plots
theme_soosh <- function(base_size = 11, base_family='serif'){
  theme_minimal(base_size, base_family) +
    theme(plot.background = element_rect(fill = 'white', colour=NULL))
}

# ==== Reading data ====
ceff <- read_csv('data/catch_effort_infilled.csv') %>% 
  # Column indicating if data was interpolated or not
  mutate(isintp = is.na(eff_hours))
field <- read_csv('data/anemone_field_data.csv')
weight <- read_csv('data/anemone_weight_data.csv')
sales <- read_csv('data/fish_market_data.csv')

# ==== General explorations ====

# How many observations
nrow(ceff)
# How many missing in each important variable
table(is.na(ceff$catch_kg))
table(is.na(ceff$eff_hours))

# Per year - missing data amounts are largely proportion to available data (i.e
# years with more data tend to have more missing obs)
ceff %>% group_by(y) %>% summarize('num_obs' = n())
ceff %>% group_by(y) %>% summarize('num_missing' = sum(is.na(eff_hours)))

# ==== Analyzing catch and effort data ====

# Total and average effort (number of hours) in each year - averages and
# standard deviations of effort are pretty consistent overall! Though total
# effort varies considerably, potentially linked to just variation in number of
# recorded observations. This can be cross-referenced with the sale data to
# confirm
ceff %>% 
  group_by(y) %>% 
  summarize('tot_hrs' = sum(eff_hours_intp),
            'avg_hrs_pday' = mean(eff_hours_intp),
            'sd_hrs_pday' = sd(eff_hours_intp))

# Hours across the time series
ceff %>% 
  ggplot(aes(x = date, y = eff_hours_intp)) +
  geom_col(show.legend = F) +
  theme_soosh() +
  # facet_wrap('year', scales = 'fixed') +
  labs(x = NULL, y = 'Fishing hours')

# Hours across the seasonality period, split by year
ceff %>% 
  ggplot(aes(x = yd, y = eff_hours_intp)) +
  geom_col(show.legend = F) +
  theme_soosh() +
  facet_wrap('y', scales = 'fixed') +
  labs(x = NULL, y = 'Fishing hours')

# Hours split by council
ceff %>% 
  ggplot(aes(x = yd, y = eff_hours_intp)) +
  geom_col(show.legend = F) +
  theme_soosh() +
  facet_wrap('council', scales = 'fixed') +
  labs(x = NULL, y = 'Fishing hours')

# Hours by fisher
ceff %>% 
  ggplot(aes(x = yd, y = eff_hours_intp)) +
  geom_col(show.legend = F) +
  theme_soosh() +
  facet_wrap('fisher_code', scales = 'fixed') +
  labs(x = NULL, y = 'Fishing hours')

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

# No repeats!
any(fishers_pday$repeat_fishers)

# How many fishers out per day
hist(fishers_pday$num_active_fishers)

# Adding number of active fishers per day to the catch/effort
ceff <- ceff %>% 
  left_join(fishers_pday %>% 
              select(-repeat_fishers),
            by = 'date')

# Distribution of effort as number of fishers across the time series
ceff %>% 
  ggplot(aes(x = date, y = num_active_fishers)) +
  geom_col() +
  theme_soosh() +
  labs(x = NULL, y = 'Number of active fishers')

# Across the seasonality period, split by year
ceff %>% 
  ggplot(aes(x = yd, y = num_active_fishers)) +
  geom_col() +
  facet_wrap('y', scales='fixed') +
  theme_soosh() +
  labs(x = NULL, y = 'Number of active fishers')

# Split by council
ceff %>% 
  ggplot(aes(x = date, y = num_active_fishers, colour = council)) +
  geom_col(show.legend = F) +
  theme_soosh() +
  facet_wrap('council', scales = 'fixed') +
  labs(x = NULL, y = 'Number of active fishers')

# Total distribution by seasonality - clear indication of a closed season
p <- ceff %>% 
  group_by(yd) %>% 
  summarize('total_af' = sum(num_active_fishers, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot(aes(x = yd, y = total_af)) +
  geom_col() +
  theme_soosh() +
  labs(x = NULL, y = 'Number of active fishers')

# Using plotly to precisely identify the date range of a potential closed season
# - days 91 to 151 (Apr 1 to May 31)
p %>% ggplotly()

# Who are the fishers operating in these dates? - many different ones at many
# different times, all in 2017, but never more than one in any given day
ceff %>% 
  filter(yd > 90 & yd < 152) %>% 
  pull(fisher_code) %>% 
  unique()

# Was fishing effort in this period concentrated by region? - nope
ceff %>% 
  filter(yd > 90 & yd < 152) %>% 
  pull(council) %>% 
  unique()

# ==== Daily dynamics of catch, effort and sales ====

# Summarizing catch and effort per day:
ceff_daily <- ceff %>% 
  # filter(!isintp) %>% 
  group_by(date) %>% 
  summarize(
    catch_kg_tot = sum(catch_kg),
    eff_h_tot = sum(eff_hours_intp),
    eff_num_fisher_tot = unique(num_active_fishers)
  )

# Summarizing sales per day
sales_daily <- sales %>% 
  group_by(date) %>% 
  summarize(
    sales_kg_tot = sum(kg),
    mean_price_pkg = mean(price_eur_kg),
    revenue_tot_eur = sum(price_tot_eur)
  )

# Filling in the time series to include missing dates. Setting newly created
# dates to have values of 0 for catch and effort (we assume no trips on these
# days, though it may just be unreported) - first creating an empty dataframe
# containing only all the dates
daily <- tibble('date' = seq(min(ceff$date), max(ceff$date), by = "day"))

# Joining all data to this dataframe
daily <- daily %>% 
  left_join(ceff_daily, by='date') %>% 
  left_join(sales_daily, by='date') %>% 
  # Replacing all NAs with 0s, to indicate no catch/effort/sales
  mutate(across(2:7, ~{
    ifelse(is.na(.x), 0, .x)
  }))

# How many days of no catch/effort
daily %>% filter(catch_kg_tot == 0) %>% nrow()

# How many days of no sales
daily %>% filter(sales_kg_tot == 0) %>% nrow()

# How many days where there were catches but no sales
daily %>% filter((sales_kg_tot == 0) & (catch_kg_tot != 0)) %>% nrow()

# And vice-versa - how many days where there were sales but no catches
daily %>% filter((catch_kg_tot == 0) & (sales_kg_tot != 0)) %>% nrow()

# Days where both were 0 - most 0 days fall in this category
daily %>% filter((catch_kg_tot == 0) & (sales_kg_tot == 0)) %>% nrow()

# Plotting relationship between two types of effort (split by year), pretty
# consistent:
daily %>% 
  mutate(y = year(date)) %>% 
  ggplot(aes(x = eff_h_tot, eff_num_fisher_tot)) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  xlab("Effort (h)") +
  ylab("Catches (kg)") +
  facet_wrap('y') +
  theme_soosh()

# Between catch and effort (again, consistent across years):
daily %>% 
  mutate(y = year(date)) %>% 
  ggplot(aes(x = eff_h_tot, catch_kg_tot)) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  xlab("Effort (h)") +
  ylab("Catches (kg)") +
  facet_wrap('y') +
  theme_soosh()

# Distribution of catch across the time-series
daily %>% 
  ggplot(aes(x = date, y = catch_kg_tot)) +
  geom_col() +
  theme_soosh() +
  labs(x = NULL, y = 'Total catch (kg)')

# Across the seasonality period, split by year
daily %>% 
  mutate(y = year(date)) %>% 
  mutate(yd = yday(date)) %>% 
  ggplot(aes(x = yd, catch_kg_tot)) +
  geom_col() +
  facet_wrap('y') +
  labs(x = NULL, y = 'Total catch (kg)') +
  theme_soosh()

# Catch spike statistic (Roa-Ureta 2015), split by year (calculated at a weekly
# timestep)
daily %>% 
  # filter(!isintp) %>% 
  mutate(y = year(date)) %>% 
  mutate(week = week(date)) %>% 
  group_by(y, week) %>% 
  summarize(
    catch_kg_tot = sum(catch_kg_tot),
    eff_h_tot = sum(eff_h_tot)
  ) %>% 
  mutate(
    catch_spike = 10*(
      (catch_kg_tot/max(.$catch_kg_tot)) - (eff_h_tot/max(.$eff_h_tot))
    )
  ) %>% 
  ggplot(aes(x = week, y = catch_spike)) +
  geom_line() +
  geom_vline(aes(xintercept = 14), colour='red', linetype='dashed') +
  geom_vline(aes(xintercept = 22), colour='red', linetype='dashed') +
  facet_wrap('y') +
  theme_soosh() +
  labs(y = 'Catch Spike', x = NULL)

# CPUE, split by year (calculated at a weekly timestep)
daily %>% 
  # filter(!isintp) %>% 
  mutate(y = year(date)) %>% 
  mutate(week = week(date)) %>%
  group_by(y, week) %>% 
  summarize(
    catch_kg_tot = sum(catch_kg_tot),
    eff_h_tot = sum(eff_h_tot)
  ) %>% 
  mutate(
    cpue = catch_kg_tot/eff_h_tot
  ) %>% 
  ggplot(aes(x = week, y = cpue)) +
  geom_line() +
  geom_vline(aes(xintercept = 14), colour='red', linetype='dashed') +
  geom_vline(aes(xintercept = 22), colour='red', linetype='dashed') +
  facet_wrap('y') +
  theme_soosh() +
  labs(y = 'Catch per unit effort', x = NULL)

# Correlation between catch and sales - not super tight
p <- daily %>% 
  mutate(y = year(date)) %>% 
  ggplot(aes(y = sales_kg_tot, x=catch_kg_tot)) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm") +
  labs(y = "Sales (kg)", x="Catches (kg)") +
  theme_soosh()
p
# Not great by year either, though some years are better than others
p + facet_wrap('y')

# ==== Weight data ====

# Patterns of weight by day in the year - no significant patterns based tentacle, disk
# colour or location
weight %>% 
  mutate(
    yd = lubridate::yday(date)
  ) %>% 
  ggplot(aes(yd, weight_g, colour = location)) +
  geom_point() +
  xlab("") +
  ylab("Indv. weight (g)") +
  theme_bw()


# Any potential relationship between size and price?
weight_price <- weight %>%  
  mutate(yd = lubridate::yday(date)) %>% 
  group_by(yd) %>% 
  summarize(mean_weight_g = mean(weight_g)) %>% 
  inner_join(sales %>% 
              mutate(yd = lubridate::yday(date)) %>% 
              group_by(yd) %>% 
              summarize(mean_price_pkg = mean(price_eur_kg),
                        sd_price_pkg = sd(price_eur_kg)))
weight_price %>% 
  ggplot(aes(x = mean_weight_g, y = mean_price_pkg, colour=sd_price_pkg)) +
  geom_point(alpha = 0.8) +
  scale_colour_distiller(palette='RdYlGn') +
  # stat_smooth(method = "lm") +
  theme_soosh()
