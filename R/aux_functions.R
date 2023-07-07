# a set of auxiliary functions

ols <- function(x,y){
  bols <- solve(crossprod(x)) %*% crossprod(x,y)
  return(bols)
}


quantU = function(df, CI = 16){
  u = quantile(df , (100 - CI) / 100)
  return(u)
}


quantL = function(df, CI = 16){
  l = quantile(df ,  CI / 100)
  return(l)
}

residuals = function(bvar, measure = "median" ){
  
  n = bvar$vardata$number_of_variables
  t = nrow(bvar$vardata$y_lhs)
  draws = bvar$draws
  res = bvar$res
  residuals = as.data.frame(matrix(0,nrow = t, ncol = n ))
  
  for (i in 1:n) {
    residuals[,i] =  apply( as.data.frame(res[,i,1:draws]), 1, median)
    names(residuals)[i] <- bvar$vardata$names_of_endogenous[i]
  }
  
  #residuals = cbind(Date, residuals)
  return(residuals)

}

