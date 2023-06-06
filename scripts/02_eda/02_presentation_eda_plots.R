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

# Path to weekly summarized dataset
weekly_path <- 'data/weekly_ceff_sale_ts.csv'

# Path to output directory for plots
out_dir <- 'output/eda-plots'

# Rapid theme function for quick plots
theme_soosh <- function(base_size = 11, base_family='serif'){
  theme_minimal(base_size, base_family) +
    theme(plot.background = element_rect(fill = 'white', colour=NA))
}

# ==== Reading data ====
weekly <- read_csv(weekly_path)

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
    y = 'Hours of fishing'
  )

# Effort - number of fishers
p3 <- weekly %>% 
  ggplot() +
  # geom_vline(xintercept=year_change, colour='darkgrey') +
  geom_line(aes(y = eff_num_fisher_tot, x = date), colour='darkgrey', alpha=0.9) +
  geom_point(aes(y = eff_num_fisher_tot, x = date), colour='black', stroke=0.3,
             alpha=0.2, size=0.3) +
  geom_vline(xintercept = closed_season, colour='red', linetype='dashed') +
  theme_soosh() +
  labs(
    x = NULL,
    y = 'Number of fishers'
  )

# CPUE
p4 <- weekly %>% 
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

# Combining plots and exporting
mosaic_hours <- p1/p2/p4
ggsave(file.path(out_dir, 'cat_hrseff_cpue_ts.png'), 
       mosaic_hours, 
       width=7, height=5, dpi=300)

# Combined plot
mosaic_fishers <- p1/p3/p4
ggsave(file.path(out_dir, 'cat_fisheff_cpue_ts.png'), 
       mosaic_fishers, 
       width=7, height=5, dpi=300)

mosaic_all <- p1/p2/p3/p4
ggsave(file.path(out_dir, 'cat_alleff_cpue_ts.png'), 
       mosaic_all, 
       width=7, height=5, dpi=300)

# Relationships between each type of effort and catch ----------

# Catch vs hours
p1 <- weekly %>% 
  ggplot(aes(y = catch_kg_tot, x = eff_h_tot)) +
  geom_point(colour='darkgrey', alpha=0.8) +
  geom_abline(intercept=0, slope=1250/250, linetype='dashed') +
  geom_smooth(method = 'lm', colour = 'firebrick', alpha=0.6) +
  theme_soosh() +
  xlim(0, 250) +
  ylim(-10, 1250) +
  labs(
    y = 'Catch (kg)',
    x = 'Effort (hours)'
  )

# Catch vs number of fishers
p2 <- weekly %>% 
  ggplot(aes(y = catch_kg_tot, x = eff_num_fisher_tot)) +
  geom_point(colour='darkgrey', alpha=0.8) +
  geom_abline(intercept=0, slope=1250/15, linetype='dashed') +
  geom_smooth(method = 'lm', colour = 'firebrick', alpha=0.6) +
  theme_soosh() +
  theme(axis.text.y = element_blank()) +
  xlim(0, 15) +
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
