set.seed(123)

library(tidyverse)

# create age groups
age_groups <- tibble(
  min_age = seq(0, 80, by=5),
  group = paste0(min_age, 'to', min_age+5)) %>%
  mutate(
    idx = 1:nrow(.),
    group = ifelse(min_age == max(min_age), paste0(max(min_age), '+'), group))

# create years
years <- 2010:2019

# create synthetic dengue case data
sim_cases <- function(age_groups, n=100, mu=15, size=10) {
  new_cases <- rnbinom(n, mu = mu, size = size)
  new_cases_int <- findInterval(new_cases, age_groups$min_age, all.inside = T)
  new_cases_df <- tibble(
    idx = as.numeric(names(table(new_cases_int))),
    cases = table(new_cases_int)) 
  new_cases <- left_join(age_groups, new_cases_df) %>%
    mutate(
      cases = ifelse(is.na(cases), 0, cases),
      cases = round(rnorm(n(), cases, 0.1*cases)))
  return(new_cases)
}

case_matrix <- matrix(nrow = nrow(age_groups), ncol=length(years), 
  dimnames=list(age_groups$group, years))
for (i in 1:length(years)) {
  case_matrix[ , i] <- sim_cases(age_groups, n=1000, mu=15, size=5)$cases
}
write.table(case_matrix, "data/cases.txt")


# create synthetic population
sim_pop <- function(age_groups, tot_pop=1e6, meanlog=log(70), sdlog=0.1) {
  age_group_int <- age_groups$min_age[2] - age_groups$min_age[1]
  mutate(age_groups,
    p_surv = (1 - plnorm(min_age + age_group_int/2, meanlog = meanlog, sdlog = sdlog)),
    mu = tot_pop * p_surv/sum(p_surv),
    n = round(rnorm(n(), mu, 0.01*mu))
  )
}

  
pop_matrix <- matrix(nrow = nrow(age_groups), ncol=length(years), 
  dimnames=list(age_groups$group, years))
for (i in 1:length(years)) {
  pop_matrix[ , i] <- sim_pop(age_groups, tot_pop = 1e6, 
    meanlog = log(70), sdlog=0.1)$n
}

write.table(pop_matrix, "data/pop.txt")

