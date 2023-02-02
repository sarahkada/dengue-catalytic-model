# code from doi https://doi.org/10.7554/eLife.45474 continued and adapted by SK, MAJ
rm(list = ls())
gc()
if (!require("pacman")) suppressMessages(install.packages("pacman"))
library(pacman)

# load packages and install if package is absent
##############################
pacman::p_load(
  tidyverse, # includes many packages for tidy data wrangling and presentation
  rstan, 
  coda,
  ggplot2,
  loo
)

options(mc.cores = parallel::detectCores()) # checks number of cores without having later to specify the cores argument
rstan_options(auto_write = TRUE) # extended packages to use stan

source("0-export_code.R")
iter <- 5000 # 50000
yearly_foi_sd <- 0.2 # 20% upper bound for average FOI - decided after comparing 20%, 25% and 30% upper bounds
life_exp <- 77
# average life expectancy in the US # life expectancy average values in the US https://www.cdc.gov/nchs/fastats/life-expectancy.htm
reporting_rate <- 0.1 # average reporting probabilities for all cases (severe + non-severe)
# stan code:
all_stan_codes <- c(
  modelS = "denguefoi_modelS.stan", 
  modelPS = "denguefoi_modelPS.stan", 
  modelP = "denguefoi_modelP.stan")

############################################################
### Catalytic model function and xtract output of interest -- PR data
cases <- read.table("data/cases.txt") # dengue cases
colnames(cases) <- gsub("X", "", colnames(cases)) # remove the X reas.csv adds to column header that have numbers
census <- read.table("data/pop.txt") # population census
colnames(census) <- gsub("X", "", colnames(census)) # remove the X reas.csv adds to column header that have numbers
j_name <- "sim" # jurisdiction name for saving

## informed FOI prior
# calculate median age at secondary infection
cases_p <- cases[1:16, ] %>%
  as_tibble() %>%
  mutate(age_gp = age_gp) %>%
  pivot_longer(!age_gp, names_to = "year", values_to = "cases") %>%
  group_by(age_gp) %>%
  summarize(n = sum(cases))
median_age <- cases_p %>%
  uncount(n) %>%
  pull(age_gp) %>%
  median()

incidence <- calc_incidence(cases_p, census, reporting_rate)
foi_mean <- logitnorm::twCoefLogitnorm(foi_init(life_exp, median_age, incidence), yearly_foi_sd)[1]
foi_sd <- logitnorm::twCoefLogitnorm(foi_init(life_exp, median_age, incidence), yearly_foi_sd)[2]

# reporting parameters
reporting_mean <- -2.2
reporting_sd <- 0.7
#  10% reporting median and 30% upper bound

for (i in 1:3) {

  # run model
  run <- run_catalytic_model(
    stan_code = all_stan_codes[i],
    cases = cases[1:16, ],
    census = census[1:16, ],
    max_age = 80,
    iter = iter,
    foi_mean = foi_mean,
    foi_sd = foi_sd,
    reporting_mean = reporting_mean,
    reporting_sd = reporting_sd
  )

  model_output <- data_extract(run$dat, run$output)
  
  saveRDS(model_output, file = paste0("output/model_", j_name, "_", 
    names(all_stan_codes)[i], ".rds"))
}

# for severe cases, same code but the FOI and reporting probability prior are different
# reporting parameters
reporting_mean <- -4.6
reporting_sd <- 0.8 #  1% reporting median and 5% upper bound
incidence <- calc_incidence(cases_severe, census, reporting_rate = .01) # 1% reporting rate of severe cases
