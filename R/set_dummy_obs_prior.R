dummy_obs_prior <- function(vardata, lamdaP = 1, tauP = 4, epsilonP = 1){
 
  
  # lamdaP = 1  #This controls the tightness of the priors on the first lag 1
  # tauP= 7;  # this controls the tightness of the priors on sum of coefficients 10*lamdaP
  # epsilonP=1;  # this controls tightness of the prior on the constant 1
  # vardata = vardata_form( data = df[,-c(3,7,6)],lags = 3, ex = df[,c(dum) ], ex_lag = df[,c("LGDP")],xlags = 3)
  # 
  X = vardata$x_lhs
  Y = vardata$y_lhs
  L = vardata$number_of_lags
  t_lags = vardata$number_of_obs_lags
  n = vardata$number_of_variables
  ex = vardata$number_of_exogenous
  muP = apply(Y, 2, mean)
  
  s_prior       <-  rep(0,vardata$number_of_variables)
  delta_prior   <-  rep(0,vardata$number_of_variables)
  
  for (i in 1:vardata$number_of_variables) {
    #Run a AR(1) and save the std of the residual for every endogenous variable
    xtemp          <-  dplyr::lag(Y[, i, drop = F])
    xtemp          <-  as.matrix(cbind("cons" = rep(1, nrow(xtemp)), xtemp)[-1,])
    ytemp          <-  as.matrix(Y[-c(1), i, drop = F])
    bols           <-  solve(crossprod(xtemp)) %*% crossprod(xtemp, ytemp)
    restemp        <-  ytemp - xtemp %*% bols
    delta_prior[i] <-  bols[2]
    s_prior[i]     <-  crossprod(restemp)/(nrow(ytemp))
  }
  
  
  deltaP = delta_prior
  sigmaP = s_prior
  
  # Initialise output (necessary for final concatenation to work when tau=0):
  y_dummy1 = x_dummy1 = y_dummy2 = x_dummy2 = NULL
  # Get dummy matrices in equation (5) of Banbura et al. 2007:
  if (lamdaP > 0) {
    
    if (epsilonP >0) { # constant
      
      # Roughly speaking, the first block of dummies imposes prior beliefs on the
      # autoregressive coefficients.
      ar_prior  <- rbind( as.matrix(diag(sigmaP*deltaP)/lamdaP),
                          matrix(0,length(sigmaP)*(L-1), length(sigmaP)) )
      
      # the second block implements the prior for the covariance matrix and the
      cov_prior <- diag(sigmaP)
      
      # third block reflects the uninformative prior for the intercept
      # (is a very small number).
      cons_prior <- matrix(0,ex, length(sigmaP))
      
      
      # the uninformative prior for the intercept (is a very small number)
      
      y_dummy1 <- rbind(cons_prior,
                        ar_prior,
                        cov_prior)
      
      jp <- diag(1:L)
      
      # eq.5 2008(ecb) left matrix
      # autoregressive coefficients.
      
      ar_prior_x  <- cbind( matrix(0,(n*L),ex) ,  jp %x% as.matrix(diag(sigmaP)/lamdaP) )
      #dim(ar_prior_x) # 33 23
      
      # covariance matrix
      cov_prior_x <-  matrix(0, n, ((n*L) + ex))
      #dim(cov_prior_x)
      
      # third block reflects the uninformative prior for the intercept
      cons_prior_x <- cbind(diag(ex)*epsilonP , matrix(0,ex, n*L))
      #dim(cons_prior_x)
      
      x_dummy1 = rbind(cons_prior_x,
                       ar_prior_x,
                       cov_prior_x)
      
      
      
      
    } else {
      
      # Roughly speaking, the first block of dummies imposes prior beliefs on the
      # autoregressive coefficients.
      ar_prior  <- rbind( as.matrix(diag(sigmaP*deltaP)/lamdaP),
                          matrix(0,length(sigmaP)*(L-1), length(sigmaP)) )
      
      # the second block implements the prior for the covariance matrix and the
      cov_prior <- diag(sigmaP)
      
      # third block reflects the uninformative prior for the intercept
      # (is a very small number).
      
      # the uninformative prior for the intercept (is a very small number)
      
      y_dummy1 <- rbind(ar_prior,
                        cov_prior)
      
      jp <- diag(1:L)
      
      # eq.5 2008(ecb) left matrix
      # autoregressive coefficients.
      ar_prior_x  <- rbind(cbind( jp %x% as.matrix(diag(sigmaP)/lamdaP) ),
                           matrix(0, n, ((n*L) + ex)) )
      
      # covariance matrix
      cov_prior_x <- matrix(0, n, ((n*L) ))
      
      # third block reflects the uninformative prior for the intercept
      x_dummy1    <- rbind(ar_prior_x,
                           cov_prior_x)
      
      
    }
    
  }
  
  #dim(x_dummy1);dim(y_dummy1)
  
  # Get additional dummy matrices - see equation (9) of Banbura et al. 2007:
  
  if (tauP > 0) {
    if (epsilonP > 0) {
      
      y_dummy2 <-  diag(deltaP * muP) / tauP
      
      x_dummy2 <- cbind( matrix(0,n,ex),  matrix(1,1,L) %x% y_dummy2 )
      
    } else {
      
      y_dummy2 <-  diag(deltaP * muP) / tauP
      
      x_dummy2 <- cbind( matrix(1,1,L) %x% y_dummy2 )
    }
  }
  
  y_dummy = rbind(y_dummy1, y_dummy2)
  x_dummy = rbind(x_dummy1, x_dummy2)
  
  #dim(x_dummy);dim(y_dummy)
  
  # yd and xd are the dummy data. Append this to actual data
  Y0=rbind(Y,y_dummy)
  X0=rbind(X,x_dummy)
 
  dummy_obs = list("Y_dummy" = Y0, "X_dummy" = X0)
  return(dummy_obs)
  
}