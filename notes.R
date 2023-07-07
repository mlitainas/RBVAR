rm(list = ls())
devtools::load_all()

usethis::use_r("BVAR_estimation_NIW.R")
load("C:\\Users\\ecs20ml\\OneDrive\\Desktop\\Data_Chapter1\\Dataset.Rdata")
library(expm)
library(tidyverse)
dummy = readxl::read_xlsx("C:\\Users\\ecs20ml\\Desktop\\dummy_proxy.xlsx") %>% mutate(Date = lubridate::ymd(Date))


dataset  = data
Dstart = "1982-01-01"
Dend= "2020-01-01"

df = dataset %>%
  filter(iso2c=="US" & Date < Dend & Date >  Dstart) %>%
  left_join(dummy)

dum =   df %>%
  select(Date) %>%
  mutate(gfc = ifelse( (Date >= "2007-10-01" & Date <= "2009-01-01") ,1,0),
         zlb = ifelse( (Date >= "2010-01-01" ) ,1,0))


#ex =   df %>%  select(lgdp_pc)
en = df %>% select(Date, NT13, lgs_pc,ltax_pc, lgdp_pc, lr,lre)

BVAR_estimation_NIW(data = en, lags = 2, const = T, dummy = NULL, reps = 2000,burn = 1000,tau = 0, lamda = 1000, epsilon = 1/100) %>%
  BVAR_irf_chol(hor = 31)



BVAR_estimation_NIW(data = en, lags = 2, const = T, dummy = dum[,2:3], reps = 2000,burn = 1000,tau = 0, lamda = 1000, epsilon = 1/100) %>%
  BVAR_irf_chol(hor = 31)






