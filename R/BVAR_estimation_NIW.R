BVAR_estimation_NIW = function(data, lags = 4, reps = 1000, burn = 500, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL, xlags = 0,dummy = dummy, lamda = 1 ,tau = 4  , epsilon = 1) {

  #lags =4; reps = 1000; burn = 500 ;const= T;trend =F; trend_qua = F;  ex = NULL ;  xlags = 4; i=1; ex_lag = NULL;lamda=1;tau=0; epsilon=0.0001
  #data = df[1:100,c("LGSPEN","GDP_OECD_IX")]
  #data =  us %>% filter(Date > "1982-01-01" & Date < "2020-01-01") %>%  select(lgs_pc, NT13, ltax_pc, lgdp_pc, lcon_pc)
  # vardata = df1$dfr[[3]]
  # ex_lag = as.data.frame(vardata) %>% select_if(~ !any(is.na(.))) %>% select(ends_with("_s"), -timeID)
  # data = as.data.frame(vardata) %>% select_if(~ !any(is.na(.))) %>% select(-ends_with("_s"), -timeID)

  vardata = vardata_form(data,lags = lags ,const = const, trend = trend, trend_qua = trend_qua, ex, ex_lag = ex_lag, xlags = xlags, dummy = dummy)

  # Take the arguments from the
  n = vardata$number_of_variables
  ex = vardata$number_of_exogenous
  L = vardata$number_of_lags
  t_lag = vardata$number_of_obs_lags
  draws = reps - burn
  dum = dummy_obs_prior(vardata = vardata, lamdaP = lamda ,tauP = tau , epsilonP = epsilon)

  # conditional mean of the VAR coefficients
  X0 = dum$X_dummy; Y0 = dum$Y_dummy

  mstar= matrix(ols(x = X0,y = Y0))  #ols on the appended data

  ixx= solve(crossprod(X0))  # inv(X0'X0) to be used later in the Gibbs sampling algorithm

  sigma= diag(n); #starting value for sigma


  CM = array(0, list(n * L , n * L, draws))
  Sigma = array(0, list(n, n, draws))
  res = array(0, list(t_lag, n, draws))
  fit = array(0, list(t_lag, n, draws))
  beta_s = array(0, list((n*L + ex),n, draws), dimnames=list( colnames(vardata$x_lhs) , vardata$names_of_endogenous , NULL ))

  for (i in 1:reps){

    if ((i%%20) == 0 ){
      print(paste0("This is my output of iteration ", i, "over",reps,"."))
    }

    vstar = sigma %x% ixx
    beta =  mstar + t(matrix(rnorm(n*(n*L+ex)), nrow = 1) %*% chol(vstar))
    #beta =  t(mvtnorm::rmvnorm(mean = mstar, sigma = vstar ,n = 1))

    # draw covariance
    e = Y0 - X0 %*% matrix(beta,nrow = n*L+ex, ncol = n)

    scale = crossprod(e)

    sigma = MCMCpack::riwish( (t_lag + nrow(Y0)) , scale)  #flag

    if (i > burn) {
      beta_s[, , i-burn ] = matrix(beta,nrow = n*L+ex, ncol = n)
      CM[, , i-burn ]  =   rbind(t(beta_s[(1+ex):nrow(beta_s), , i-burn]), cbind(diag(n*L-n), matrix(0,n*L-n,n)))
      Sigma[, , i-burn] = sigma
      res[, , i-burn ] =  e[1:t_lag,]
      #fit[, , i] = bvar_hat$y_fit
    }


  }

  arguments = list("tau" = tau, "lamda" = lamda, "epsilon" =  epsilon)
  bvar = list(beta_s, Sigma, CM, res, fit, vardata, draws, arguments)
  names(bvar) = c("posterior","Sigma", "CM", "res", "yfit", "vardata", "draws","priors")
  return(bvar)

}
