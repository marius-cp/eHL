rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach);library(tidyverse);library(doParallel);library(kableExtra)
source("../nle.R")
source("../eHL.R")
source("../oracle_eHL.R")
MCsim <-
  readRDS("sim_HLsytle_oos.rds")

check <- sapply(MCsim, function(x) inherits(x, 'error'))
out <- MCsim[!check]
MCsim[check]
dat <-do.call("rbind", out)
dat


dt <-
  dat %>%
  #filter(test=="cHL") %>%
  mutate(
    S=as.factor(S),
    value = ifelse(test =="cHL",  (1 - pchisq(value, 10)),value),
    dec = if_else(
      test=="cHL",
      as.numeric(value<=0.05),
      as.numeric(value>=20))
  ) %>%
  group_by(obs,J,test,S,setup) %>%
  summarize(rej = mean(dec)) %>%
  dplyr::filter(
    setup == "quadratic",
    #test == "eHL" && S%in%c("0.2","0.4","0.6","0.8",NA,0,1)
    #!(test=="eHL" && S%in%c("0.3","0.5","0.7")),
    !(test=="oracle_eHL" & S=="0.6"),
    obs<9000,obs>1000,
    J>=.007
  )

sizes <-
  dt %>%
  ungroup() %>%
  filter(J == 0.007) %>%
  select(obs, test, rej, S) %>%
  spread(test, rej)

sizes


tbl <-
  cbind(
    sizes %>%
      select(-"oracle_eHL", -"eHL") %>%
      filter(S%in%c(NA)) %>%
      select(-S) %>%
      mutate(cHL=100*cHL),
    pivot_wider(
      sizes %>%
        select(-"cHL", -"eHL") %>%
        filter(complete.cases(.)),
      names_from= S,
      values_from = oracle_eHL) %>%
      select(-obs)*100,
    pivot_wider(
      sizes %>%
        select(-"cHL", -"oracle_eHL") %>%
        filter(complete.cases(.)),
      names_from= S,
      values_from = eHL
    ) %>%
      select(-obs)*100
  )

kbl(
  data.table::data.table(tbl),
  booktabs = T,
  digits = 2,
  'latex'
) %>%
  kable_styling()%>%
  add_header_above(c(" " = 1, " " = 1, "oracle eHL" = 2, "eHL" = 5))



dat %>%
  mutate(
    S=as.factor(S),
    dec = if_else(
      test=="cHL",
      as.numeric(value<=0.05),
      as.numeric(value>=20))
  ) %>% filter(test=="cHL")
