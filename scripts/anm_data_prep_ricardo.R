# Author: Ricardo González-Gil, adapted by Saeesh Mangwani
# Date: 2023-04-14

# Description: NOTE: This scripts prepares the anemone dataset for exploratory
# analysis + modelling with CatDyn. It is based on Rubén Urea-Rota's script
# AnemonaBD.R

# ==== Libraries ====
library(tidyverse)
library(mgcv)
library(readr)

# ==== Paths and global variables ====

# Path to output plots from this script
pout_dir <- 'output/ricardo-dataprep'

# ==== Loading data ====
load("data/raw/cateff_codes.RData") 
load("data/raw/fish_markets_codes.RData") 
load("data/raw/weight_data.RData") 
load("data/raw/field_data.RData") 

# Writing all raw dataframes to disk
write_csv(cateff_codes, 'data/catch_effort.csv')
write_csv(fish_markets_codes, 'data/fish_market_data.csv')
write_csv(weight_data, 'data/anemone_weight_data.csv')
write_csv(field_data, 'data/anemone_field_data.csv')

# Adding new time related variables
cateff_codes_2 <- cateff_codes %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  relocate(y, m, d, yd, .after = date)

# ==== Exploring data completeness ====

# Once we have all the data and the new time variables, let's inspect whether
# there are several NA effort values that should be estimated. This avoids
# losing information.
cateff_codes_2 %>% 
  filter(is.na(eff_hours))

# How many missing per fisher
cateff_codes_2 %>% 
  group_by(fisher_code) %>% 
  summarize(num_missing = sum(is.na(eff_hours))) %>% 
  arrange(desc(num_missing))

# We can observe that there are many effort data that are NA. To estimate them,
# we can use day of year as predictor variable. In particular, we can fit a
# smooth curve.

# We can be tempted to use catches as predictor of effort too, but then this
# introduces an artificial dependency between effort and catches that will
# derive into an error when fitting the generalized depletion models (GDMs).

# First, let's do a visual exploration of these options.

# Seasonal cycle
seas_eff_hours_anemo <- cateff_codes_2 %>% 
  ggplot(aes(x=yd, y=eff_hours)) +
  geom_point() +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0), name = "Day of year") +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cc")) +
  ylab("Effort (hours)") +
  theme_bw()
seas_eff_hours_anemo

ggsave(file.path(pout_dir, 'seas_eff_hours.png'), 
       plot = seas_eff_hours_anemo, 
       width = 7, height = 5, 
       dpi = 300)


# Any regional effect?
seas_eff_hours_regional <- cateff_codes_2 %>% 
  ggplot(aes(yd, eff_hours, colour=council)) +
  geom_point(alpha=0.6, colour='lightgrey') +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0)) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10, bs = "cc")) +
  facet_wrap(~council) +
  theme_bw() +
  guides(colour = 'none')
ggsave(file.path(pout_dir, 'seas_eff_hours_region.png'), 
       plot = seas_eff_hours_regional,
       width = 7, height = 5, 
       dpi = 300)

# Any effects by fisher?
cateff_codes_2 %>% 
  ggplot(aes(yd, eff_hours, colour=fisher_code)) +
  geom_point(alpha=0.6, colour='lightgrey') +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0)) +
  # stat_smooth(method='loess', formula = y ~ x) +
  stat_smooth(method = "gam", formula = y ~ s(x, bs = "cc")) +
  facet_wrap(~fisher_code) +
  theme_bw() +
  guides(colour = 'none')
ggsave(file.path(pout_dir, 'seas_eff_hours_fisher.png'), 
       width = 7, height = 5, 
       dpi = 300)

# There might be some differences among councils, but I don't think they're
# important. We could formally test this with a selection model process (or
# other stat test), but I'm also afraid that using council as a predictor
# variable could introduce some artificial effect and thus, I'm going to just
# use day of year as predictor, as Rubén did.

# The model
gam_yd <- gam(eff_hours ~ s(yd, bs = "cc", k = 10), 
              knots = list(yd = c(0, 366)), 
              method = "REML", 
              data = cateff_codes_2)

# Checking the model
par(mfrow = c(2, 2))
gam.check(gam_yd)

# For the predictions (average curve):
gam_yd_pred_eff <- data.frame(yd = 1:366) %>% 
  mutate(
    fit_eff = as.numeric(predict(gam_yd, type = "response", 
                                 se.fit = TRUE, 
                                 newdata = .)$fit),
    sd_eff = as.numeric(predict(gam_yd, type = "response", 
                                se.fit = TRUE, 
                                newdata = .)$se.fit)
    # I call it sd (standard deviation) because: "The term "standard error" is
    # used to refer to the standard deviation of various sample statistics, such
    # as the mean or median. For example, the "standard error of the mean"
    # refers to the standard deviation of the distribution of sample means taken
    # from a population. The smaller the standard error, the more representative
    # the sample will be of the overall population." --> see
    # https://www.investopedia.com/terms/s/standard-error.asp#:~:text=The%20standard%20error%20(SE)%20is,known%2C%20or%20accepted%20as%20accurate.
  )

# Adding the average cycle of efforts to original data frame
cateff_codes_3 <- cateff_codes_2 %>% 
  left_join(., gam_yd_pred_eff)

# Create a vector of re-sampled mean efforts using day of year, the GAM model and 
# adding random error to account for using a spline model instead of raw data.
library(Runuran)
set.seed(12) # To always generate the same values

cateff_codes_3$resampled_eff <- NA
pivot <- 2

for(i in 1:length(cateff_codes_3$yd)){
  cateff_codes_3$resampled_eff[i] <- urnorm(
    1, 
    cateff_codes_3$fit_eff[i], 
    cateff_codes_3$sd_eff[i], 
    cateff_codes_3$fit_eff[i] - pivot * cateff_codes_3$sd_eff[i], 
    cateff_codes_3$fit_eff[i] + pivot * cateff_codes_3$sd_eff[i]
  )
}

cateff_codes_3 %>% 
  arrange(date, fisher_code, region) %>% 
  select(-council, -weighting_location, -y, -m, -d) %>%
  View()

# Let's check how the resampled effort looks like:
pred_and_raw_eff <- data.frame(
  date = seq(
    min(cateff_codes_3$date), 
    max(cateff_codes_3$date), 
    by = "day"
  )
) 

pred_and_raw_eff %>% 
  left_join(., cateff_codes_3, multiple = "all") %>% 
  ggplot(aes(date, eff_hours)) +
  geom_point(color = "lightblue") +
  geom_point(aes(y = resampled_eff), color = "red") +
  ylab("Effort (hours)") +
  ggtitle("Predicted + resampled (red lines) and raw (blue dots) effort data") +
  theme_bw()
ggsave(file.path(pout_dir, 'pred_and_raw_eff.png'), 
       width = 7, height = 5, 
       dpi = 300)

# Now, let's replace the effort values in those cases with NA
cateff_codes_replace_final <- cateff_codes_3 %>% 
  mutate(eff_hours_all = ifelse(is.na(eff_hours), 
                                round(resampled_eff, digits = 1), 
                                eff_hours))

# Let's check this:
View(cateff_codes_replace_final %>% filter(is.na(eff_hours)))

# Finally, let's save this data frame
cateff_codes_replace_final %>% 
  rename(eff_hours_intp = eff_hours_all) %>% 
  write_csv('data/catch_effort_infilled.csv')
# save(cateff_codes_replace_final, file = "cateff_codes_replace_final.RData")

# It's important to highlight that we can consider different types of effort,
# the total number of hours per day and the total number of "trips" per day
# (this would be the length of the number of fishers per day).

# Let's estimate the data per day. We could also chose another time step unit
# such as week (I think it's better using isoweek).

# 1) Is any name of a fisher repeated twice in a single date?
fisher_uniq_names_date <- cateff_codes_replace_final %>% 
  group_by(date, fisher_code) %>% 
  summarise(num_uniq_fisher = length(fisher_code)) %>% 
  ungroup()

any(fisher_uniq_names_date$num_uniq_fisher > 1)
# [1] FALSE
# This indicates that there is no repetitions of names in a single date.

# 2) Estimation of variables per day (we could also include other variables).
# It's important to take into account that, for those dates with no catch data,
# we should replace it with 0
cateff_day_1 <- cateff_codes_replace_final %>%
  group_by(date) %>% 
  summarise(
    catch_kg_tot = sum(catch_kg),
    eff_h_tot = sum(eff_hours_all),
    eff_num_fisher_tot = length(fisher_code)
  ) %>% 
  ungroup()

cateff_day <- data.frame(
    date = seq(min(cateff_codes_3$date), max(cateff_codes_3$date), by = "day")
  ) %>% 
  left_join(., cateff_day_1) %>% 
  mutate(
    catch_kg_tot = ifelse(is.na(catch_kg_tot), 0, catch_kg_tot),
    eff_h_tot = ifelse(is.na(eff_h_tot), 0, eff_h_tot),
    eff_num_fisher_tot = ifelse(is.na(eff_num_fisher_tot), 0, eff_num_fisher_tot)
  )

# Let's inspect the relationship between catches and effort per day.
cateff_day %>% 
  ggplot(aes(eff_h_tot, catch_kg_tot)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Effort (h)") +
  ylab("Catches (kg)") +
  theme_bw()

cateff_day %>% 
  ggplot(aes(eff_num_fisher_tot, catch_kg_tot)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Effort (# fihsers)") +
  ylab("Catches (kg)") +
  theme_bw()

# Relationship between 2 types of effort
cateff_day %>% 
  ggplot(aes(eff_h_tot, eff_num_fisher_tot)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Effort (h)") +
  ylab("Effort (# fihsers)") +
  theme_bw()

# Exploration of the relationship between the total captures and total sales per
# day.
fish_market_day_1 <- fish_markets_codes %>%
  group_by(date) %>% 
  summarise(
    sales_tot_kg = sum(kg),
    price_eur_kg = mean(price_eur_kg),
    revenue_tot_eur = sum(price_tot_eur)
  ) %>% 
  ungroup()

fish_market_day <- data.frame(
    date = seq(min(cateff_codes_3$date), max(cateff_codes_3$date), by = "day")
  ) %>% 
  left_join(., fish_market_day_1) %>% 
  mutate(
    sales_tot_kg = ifelse(is.na(sales_tot_kg), 0, sales_tot_kg),
    price_eur_kg = ifelse(is.na(price_eur_kg), 0, price_eur_kg),
    revenue_tot_eur = ifelse(is.na(revenue_tot_eur), 0, revenue_tot_eur)
  )

# Catch and sales data per day together:
cateff_sales_day <- 
  cateff_day %>% 
  left_join(., fish_market_day)

# Visually inspecting both time series.
cateff_sales_day %>% 
  pivot_longer(cols = c(catch_kg_tot, sales_tot_kg), 
               names_to = "var", values_to = "val") %>% 
  ggplot(aes(date, val)) +
  geom_point(aes(color = var)) +
  scale_color_discrete(name = "Variable", labels = c("Catches", "Sales")) +
  scale_x_date(expand = c(0, 0)) +
  xlab("") +
  ylab("Biomass (kg)") +
  theme_bw()

# Changes throught time by month and year
cateff_sales_day %>% 
  pivot_longer(cols = c(catch_kg_tot, sales_tot_kg), names_to = "var", values_to = "val") %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  group_by(m, y, var) %>% 
  summarize(
    val_tot = sum(val)
  ) %>% 
  ggplot(aes(factor(m), val_tot)) +
  geom_bar(aes(fill = var), stat = "identity", color = "black") +
  scale_fill_discrete(name = "Variable", labels = c("Catches", "Sales")) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  xlab("Month") +
  ylab("Biomass (kg)") +
  facet_grid(var ~ y) +
  theme_bw() +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

# Trends in the kg per month:
cateff_sales_day %>% 
  pivot_longer(cols = c(catch_kg_tot, sales_tot_kg), names_to = "var", values_to = "val") %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  group_by(m, y, var) %>% 
  summarize(
    val_tot = sum(val)
  ) %>% 
  ggplot(aes(y, val_tot)) +
  geom_line(aes(group = var)) +
  geom_point(aes(fill = var), stat = "identity", color = "black", shape = 21) +
  scale_fill_discrete(name = "Variable", labels = c("Catches", "Sales")) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  xlab("Month") +
  ylab("Biomass (kg)") +
  facet_wrap(~factor(m)) +
  theme_bw() +
  theme(
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Visually inspecting the relationship.
cateff_sales_day %>% 
  ggplot(aes(catch_kg_tot, sales_tot_kg)) +
  geom_point() +
  stat_smooth(method = "lm") +
  xlab("Catches (kg)") +
  ylab("Fish market sales (kg)") +
  theme_bw()

# Some quick exploration of the individual weights. Data from the field.

# Direct measures of the individual weights
weight_data %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  ggplot(aes(date, weight_g)) +
  geom_point() +
  xlab("") +
  ylab("Indv. weight (g)") +
  theme_bw()

weight_data %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  ggplot(aes(yd, weight_g)) +
  geom_point() +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0)) +
  xlab("Day of year") +
  ylab("Indv. weight (g)") +
  theme_bw()


# Several individuals at the same time per sample
field_data %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  mutate(w_g = weight_g / num_indv2) %>% 
  ggplot(aes(date, w_g)) +
  geom_point() +
  xlab("") +
  ylab("Indv. weight (g)") +
  theme_bw()

field_data %>% 
  mutate(
    y = lubridate::year(date),
    m = lubridate::month(date),
    d = lubridate::day(date),
    yd = lubridate::yday(date)
  ) %>% 
  mutate(w_g = weight_g / num_indv2) %>% 
  ggplot(aes(yd, w_g)) +
  geom_point() +
  scale_x_continuous(limits = c(1, 366), expand = c(0, 0)) +
  xlab("Day of year") +
  ylab("Indv. weight (g)") +
  theme_bw()



