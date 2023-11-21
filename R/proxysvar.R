rm(list = ls())
library(tidyverse)
library(dplyr)
library(expm)
source("G:\\My Drive\\Chapter 3\\GVAR\\call_functions.R")

#data = readxl::read_xlsx("C:\\Users\\micha\\OneDrive\\Desktop\\MR2014_data.xlsx", sheet = "Sheet2")
#save(data, file = "D:\\Mike Tree\\RBVAR\\data\\data.rda")
load( file = "D:\\Mike Tree\\RBVAR\\data\\data.rda")

vardata = data %>% 
  select(Date, Tax:Output,m = Tax_m)

lags = 4
m = vardata %>% select(Date, m ) 


mod = vardata %>% filter(Date > "1979-12-01") %>% 
  select(-m) %>% # print(n = Inf)
  BVAR_estimation_NIW(lags = 4,reps = 2000,burn = 1000, lamda = 100000, tau = 0, epsilon = 0.001, trend = T , trend_qua = T)


mod %>% BVAR_irf_chol(shock = 2,normalised = 2)

BVAR_psvar_irf = function(bvar, m = NULL, hor = 20, instrumented = 1 ){
  hor = 20
  if (is.null(m)) {
    stop("Hey, you forgot the intrument 😕")
  }
  
  if (any(is.na(m))) {
    stop("Hey, you instrument contains NAs 😕")
  }
  
  bvar  <- mod 
  date  <- bvar$vardata$TimeID  
  res   <- bvar$res
  nlags <- bvar$vardata$number_of_lags
  nvar  <- mod$vardata$number_of_endogenous
  draws <- mod$draws
  
  mxx = m[m$Date %in% date,   ]
  # Variable, Shock Horizon
  VSH = array(0,  list(nvar, 1, hor))
  
  # Save final IRFs for each posterior 
  irf = array(0, list(hor, nvar, draws))
  
  pv = rep(0,draws)
  
  for (i in 1:draws) {
  i=1
    res_temp =  res[,,i]
    dim(res_temp)
    dim(m)
      
    
    epsilon_p = res_temp[, 1,drop =F]
    epsilon_q = res_temp[, -1]
      
    #regress the reduced form error on the instrument
    fit =  lm(epsilon_p ~  m   - 1 ) 
    u_hat_p = fit %>% fitted() %>% matrix(ncol = 1)
    sss = fit %>% summary()
    pv[i] = sss$fstatistic[1]
      
    epsilon_q = epsilon_q[-(1:nlags),]
    # this is the correlation between the aproxximated structural shock and the n-1 residuals
    sq_sp = solve(crossprod(u_hat_p)) %*% crossprod(u_hat_p, epsilon_q)
    s = c(1, sq_sp)
    
    print(paste("that is the ", i," th draw") )
      
    CM = mod$CM[,,1]
    [, , 1] =  (CM %^% 0)[1:nvar, 1:nvar] %*% s
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

}
