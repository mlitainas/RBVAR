rm(list = ls())
library(tidyverse)
library(ggpmisc)
library(dynlm)
source("G:\\My Drive\\Chapter 3\\GVAR\\call_functions.R")
df = readxl::read_xlsx("C:\\Users\\ecs20ml\\OneDrive\\Desktop\\SSH ProxySVAR\\Dataset.xlsx", sheet = "Dataset")

data = filter(df, Date < "2020-01-01") %>% 
  mutate(ram1 = lag(RAMEY_PCT,1),ram2 = lag(RAMEY_PCT,2), ram3 = lag(RAMEY_PCT,3),ram4 = lag(RAMEY_PCT,4),
         ram5 = lag(RAMEY_PCT,5),ram6 = lag(RAMEY_PCT,6),ram7 = lag(RAMEY_PCT,7),ram8 = lag(RAMEY_PCT,8),
         KILIAN_abs = abs(KILIAN),
         HAMILTON_abs = abs(HAMILTON),
         RAMEY_PCT_abs = abs(RAMEY_PCT),
         RAMEY_abs = abs(RAMEY),
         LEEPER_abs = abs(LEEPER),
         ROMERS_abs = abs(ROMERS),
         MR_abs = abs(MR)
  )


DStart = "1985-01-01"
DEnd = "2016-01-01"
lags = 4

vardata = data %>% 
  filter(Date > DStart & Date  <  DEnd) %>% 
  drop_na(KANZIG) %>% 
  select(Date, s9, UR, LGDP_PC, TBILL, z = KILIAN_abs)


z = vardata %>% select(Date, z ) 
z = as.matrix(z[-c(1:lags),2]) 
z= z %>% abs


mod = vardata %>% 
  select(-z) %>% # print(n = Inf)
  BVAR_estimation_NIW(lags = 4,reps = 2000,burn = 1000, lamda = 1, tau = 4, epsilon = 0.1 )

mod %>% BVAR_irf_chol()

hor = 40
draws = mod$draws
nlags = mod$vardata$number_of_lags
nvar = mod$vardata$number_of_endogenous

HDP = array(0,  list(nvar, 1, hor))
irf = array(0, list(hor, nvar, draws))
pv = rep(0,draws)
lags = 4
for (i in 1:draws) {
  
  res_b =  mod$res[,,i]
  dim(res_b)
  dim(z)
  
  #vc = crossprod(res_b)
  eps_p = res_b[, 1,drop =F]
  eps_q = res_b[, -1]
  
  #regress the reduced form error on the instrument
  fit =  lm(eps_p ~  z +  lag(z,1) +  lag(z,2)  +  lag(z,3) +  lag(z,4) ) 
  u_hat_p = fit %>% fitted() %>% matrix(ncol = 1)
  sss = fit %>% summary()
  pv[i] = sss$fstatistic[1]
  
  eps_q = eps_q[-(1:lags),]
  # this is the correlation between the aproxximated structural shock and the n-1 residuals
  sq_sp = solve(crossprod(u_hat_p)) %*% crossprod(u_hat_p, eps_q)
  s = c(1, sq_sp)

  print(paste("that is the ", i," th draw") )
  
  CM = mod$CM[,,1]
  HDP[, , 1] =  (CM %^% 0)[1:nvar, 1:nvar] %*% s
  HDP[, , 2] =   CM[1:nvar,1:nvar]   %*%  s
  
  CMhix = CM
  
  for (j in 3:hor) {
    CMhix = CMhix %*% CM
    HDP[, , j] = CMhix[1:nvar, 1:nvar]  %*%  s
  }
  
  irf[, , i] = aperm(HDP, perm = c(3, 1, 2))
  
}

plot(pv)
# News
irf_M_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))
irf_U_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))
irf_L_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))
i=1
for (i in 1:nvar) {
  #View(apply(irf[ ,i,], 1, "median"))
  irf_M_n[, i] = as.matrix(apply(irf[, i, 1:draws], 1, "median"))
  names(irf_M_n)[i] = paste(mod$vardata$names_of_endog_variables[i], ",M", sep = "")
  
  irf_U_n[, i] = as.matrix(apply(irf[, i, 1:draws], 1, "quantU"))
  names(irf_U_n)[i] = paste(mod$vardata$names_of_endog_variables[i], ",U", sep = "")
  
  irf_L_n[, i] = as.matrix(apply(irf[, i, 1:draws], 1, "quantL"))
  names(irf_L_n)[i] = paste(mod$vardata$names_of_endog_variables[i], ",L", sep = "")
}

irf_final = as.data.frame(cbind("Horizon" = 1:hor, irf_M_n, irf_L_n, irf_U_n))


irf_final %>% 
  pivot_longer(-Horizon,names_to = "Name", values_to = "Value" ) %>% 
  arrange(Name, Horizon) %>%
  separate(Name,into = c("Variable", "Quant"), sep = ",",remove = F) %>% 
  pivot_wider(id_cols = c("Variable", "Horizon"), names_from = "Quant", values_from = "Value" ) %>% 
  mutate(Horizon = Horizon - 1,
         Variable = factor(Variable, levels = mod$vardata$names_of_endog_variables),
         pID = str_c(Horizon, Variable)) %>% #c("lgr","lge","ramey","lgdp", "eq"))) %>% 
  ggplot()+
  geom_ribbon(aes(ymin = L, ymax = U, x = Horizon), alpha = 0.5)+
  geom_line(aes(x = Horizon, y = M), color = "blue", size = 1)+
  geom_hline(yintercept = 0, color  = "red")+
  facet_wrap(facets = "Variable", ncol = 2, scales = "free_y")+
  theme_minimal()

