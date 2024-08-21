BVAR_irf_proxy = function(bvar = NULL, m = NULL, hor = 20, instrumented = 1 , pcol = 2){
   bvar = mod
   m = mdf
   hor = 20

   if (any(is.null(bvar))) {
     stop("Provide a bvar object.")
   }

   if ( !is.null(bvar$vardata$TimeID) ) {
     print(paste("Time ID series identified for the BVAR object") )
   } else {
     stop("BVAR_irf_proxy requires a Time ID series during the BVAR estimation. Re-run the model providing a Date series")
   }


  if (is.null(m)) {
    stop("Hey, you forgot the intrument 😕")
  }

  if ( names(m)[1] %in% c("Date", "date", "Time", "time") ) {
     print(paste("Time ID series identified for the instrument") )
  } else {
     stop("BVAR_irf_proxy requires a Time ID series for the instrument.")
  }



  date  <- as.data.frame(bvar$vardata$TimeID)

  date = date %>% separate(col = "Date", into = c("y", "m", "d"),sep = "-")

  date = date %>%
    mutate(Date = paste(y,m,d,sep = "/") ,
           Date = lubridate::ymd(Date)) %>%
    select(Date)

  res   <- bvar$res
  nlags <- bvar$vardata$number_of_lags
  nvar  <- bvar$vardata$number_of_endogenous
  draws <- bvar$draws

  mxx = left_join(date, m )
                # Variable, Shock Horizon
  HDP = array(0,  list(nvar, 1, hor))

  # Save final IRFs for each posterior
  irf = array(0, list(hor, nvar, draws))

  pv = rep(0,draws)

  z = mxx[!is.na(mxx$m),2] %>% as.matrix()

  for (i in 1:draws) {

    res_temp =  res[!is.na(mxx$m),,i]
    # dim(res[,,1])
    # dim(z)
    # dim(bvar$vardata$y_lhs)
    epsilon_p = res_temp[, 1,drop =F]
    epsilon_q = res_temp[, -1]

    #regress the reduced form error on the instrument
    fit =  lm(epsilon_p ~  z   - 1 )
    u_hat_p = fit %>% fitted() %>% matrix(ncol = 1)
    sss = fit %>% summary()
    pv[i] = sss$fstatistic[1]

    #epsilon_q = epsilon_q[-(1:),]
    # this is the correlation between the aproxximated structural shock and the n-1 residuals
    sq_sp = solve(crossprod(u_hat_p)) %*% crossprod(u_hat_p, epsilon_q)
    s = c(1, sq_sp)

    print(paste("Draw:", i, " in ",draws) )

    CM = bvar$CM[,,1]
    HDP[, , 1] =  (CM %^% 0)[1:nvar, 1:nvar] %*% s
    HDP[, , 2] =   CM[1:nvar,1:nvar]   %*%  s

    CMhix = CM

    for (j in 3:hor) {
      CMhix = CMhix %*% CM
      HDP[, , j] = CMhix[1:nvar, 1:nvar]  %*%  s
    }

    irf[, , i] = aperm(HDP, perm = c(3, 1, 2))

  }

  #plot(pv)
  # News
  irf_M_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))
  irf_U_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))
  irf_L_n = as.data.frame(matrix(0, nrow = hor, ncol = nvar))

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


  irf_proxy <- irf_final %>%
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

  return(irf_proxy)

}
