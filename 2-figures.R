## plots
library(tidyverse)
# some figs using model S as an example
j_name <- "sim"
model <- "modelS"
this_model <- readRDS(file = paste0("output/model_", j_name, "_", model, ".rds"))

# plot cases
epi_data <- this_model$data$secondary_cases %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "cases")
ggplot(epi_data, aes(x = year, y = cases), col = "steelblue", alpha = .6) +
  geom_col(stat = "identity") +
  theme_classic()

# plot cases by age group
epi_data_by_age_group <- this_model$data$secondary_cases %>%
  mutate(age_group = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
    "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79"
  )) %>%
  pivot_longer(cols = -age_group, names_to = "year", values_to = "cases")

ggplot(transform(epi_data_by_age_group,
  age_group = factor(age_group, levels = c(
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
    "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79"
  ))
),
aes(x = age_group, y = cases)) +
  geom_col(fill="steelblue",  alpha = .6) +
  facet_wrap(~year, ncol = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) # 45 degree angle axis text and smaller text

# prep FOI data
max_A <- 80
last_year_of_data <- 2019
prep_foi <- function(data1) {
  # get years
  year <- seq(last_year_of_data, by = -1, length.out = max_A + data1$data$T)
  # date <- colnames(data1$data$secondary_cases)
  date <- paste0(year, "-01-01")

  # lambda
  lambda <- data1$lambda %>%
    mutate(year = rev(date)) %>%
    mutate(location = data1$name)
  #
  return(lambda)
}

# plot FOI
foi <- prep_foi(this_model)
ggplot(foi) +
  geom_line(aes(x = as.Date(year), y = lambda_median)) +
  geom_ribbon(
    aes(x = as.Date(year), ymin = lambda_50_ci_lr, ymax = lambda_50_ci_ur),
    alpha = .6,
    fill = "steelblue"
  ) +
  geom_ribbon(
    aes(x = as.Date(year), ymin = lambda_ci_lr, ymax = lambda_ci_ur),
    alpha = .4,
    fill = "steelblue"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme_classic()

