
nvar = bvar$vardata$number_of_variables
p = LAGS
tobs = bvar$vardata$number_of_obs_lags - ident_sample_diff
HDshock_save = array(0,dim = list(tobs, nvar, draws))
HDshock_save %>% dim()

for (i in 1:draws) {
  # shock
  res_temp =  res[-c(1:ident_sample_diff),,i]
  epsilon_p = res_temp[, 1,drop =F]
  epsilon_q = res_temp[, -1]
  
  fit =  lm(epsilon_p ~  z ) 
  u_hat_p = fit %>% fitted() %>% matrix(ncol = 1)
  sss = fit %>% summary()
  pv[i] = sss$fstatistic[1]
  sq_sp = solve(crossprod(u_hat_p)) %*% crossprod(u_hat_p, epsilon_q)
  s = c(1, sq_sp)
  
  varComp = bvar$CM[,,i]
  V  = u_hat_p 
  nvarXeq = nvar * p
  
  # Contribution of the shock
  A0_big = matrix(0, nvarXeq, nvar)
  A0_big[1:nvar, 1] = s
  Icomp = cbind( diag(nvar), matrix(0, nvar, (p-1)*nvar ) )
  HDshock_big = matrix(0, p*nvar, tobs+1 ) 
  HDshock_temp = matrix(0, nvar, tobs+1 ) 
  V_big = matrix(0,nvar, tobs+1 )
  V_big[1,-1] = V
  
  for(j in 2:(tobs+1)) {
    HDshock_big[,j] = A0_big %*% V_big[,j] + varComp %*% HDshock_big[,j-1]
    HDshock_temp[,j] = Icomp %*% HDshock_big[,j]
  }
  
  HDshock_save[,,i] =  HDshock_temp[,-1] %>% t()
}


# initialize 3 data frames that will take the irf percentiles
hd_M = hd_U = hd_L= as.data.frame(matrix(0, nrow = tobs, ncol = nvar))

for (i in 1:nvar){
  hd_M[,i] = as.matrix(apply(HDshock_save[,i,1:draws], 1, "median"))
  names(hd_M)[i] = paste(bvar$vardata$names_of_endogenous[i],",M", sep = "")
  
  hd_U[,i] = as.matrix(apply(HDshock_save[,i,1:draws], 1, "quantU"))
  names(hd_U)[i] = paste(bvar$vardata$names_of_endogenous[i],",U", sep = "")
  
  hd_L[,i] = as.matrix(apply(HDshock_save[,i,1:draws], 1, "quantL"))
  names(hd_L)[i] = paste(bvar$vardata$names_of_endogenous[i],",L", sep = "")
}

hd_final = as.data.frame(cbind("timeHD" = 1:tobs,hd_M, hd_L, hd_U))

hd_date = data.frame("timeHD" = 1:nrow(hd_final), "Date" = as.Date(Date[-c(1:(ident_sample_diff+p)),1]) )
