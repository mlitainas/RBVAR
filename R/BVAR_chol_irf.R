BVAR_irf_chol = function(bvar, hor = 36, shock = 1, pcol=2, normalised = NULL, sign = 1){
  
  # bvar = bvar1;hor = 36; shock = 1;draws = 1000;i=1;j=1; pcol=2
  CM = bvar$CM
  Sigma = bvar$Sigma
  n = bvar$vardata$number_of_variables
  draws = bvar$draws
  # initialise an array that will take the impulse responses
  #Response, shock,Horizon, Draw
  D = array(0, list(n,n, hor, draws))
  for (i in 1:draws) {
    # Cholesky IRFs
    for (j in 1:(hor)) {
      D[, , j, i] = (CM[, ,i] %^% (j - 1))[1:n, 1:n] %*%  t(chol(Sigma[,,i]))
    }
  }
  
  # permutate to Horiazon, Response, Shock, Draw
  irf = aperm(D , perm = c(3,1,2,4))
  names(irf[,1:n,,]) = bvar$vardata$names_of_endogenous
  
  
  # initialize 3 data frames that will take the irf percentiles
  irf_M = as.data.frame(matrix(0, nrow = hor, ncol = n))
  irf_U = as.data.frame(matrix(0, nrow = hor, ncol = n))
  irf_L = as.data.frame(matrix(0, nrow = hor, ncol = n))
  
  
  for (i in 1:n){
    irf_M[,i] = as.matrix(apply(irf[1:hor,i,shock,1:draws], 1, "median"))
    names(irf_M)[i] = paste(bvar$vardata$names_of_endogenous[i],",M", sep = "")
    
    irf_U[,i] = as.matrix(apply(irf[1:hor,i,shock,1:draws], 1, "quantU"))
    names(irf_U)[i] = paste(bvar$vardata$names_of_endogenous[i],",U", sep = "")
    
    irf_L[,i] = as.matrix(apply(irf[1:hor,i,shock,1:draws], 1, "quantL"))
    names(irf_L)[i] = paste(bvar$vardata$names_of_endogenous[i],",L", sep = "")
  }
  
  irf_final = as.data.frame(cbind("Horizon" = 1:hor,irf_M,irf_L, irf_U))
  irf_final[,-1] = irf_final[,-1]*sign
  
  if (!is.null(normalised)) {
    irf_final[, -1] = irf_final[, -1] / max(irf_M[, normalised])
  }  
  
  cholirf = irf_final %>% 
    pivot_longer(-Horizon,names_to = "Name", values_to = "Value" ) %>% 
    arrange(Name, Horizon) %>%
    separate(Name,into = c("Variable", "Quant"), sep = ",",remove = F) %>% 
    pivot_wider(id_cols = c("Variable", "Horizon"), names_from = "Quant", values_from = "Value" ) %>% 
    mutate(Horizon = Horizon - 1,
           Variable = factor(Variable, levels = bvar$vardata$names_of_endogenous),
           pID = str_c(Horizon, Variable)) %>% #c("lgr","lge","ramey","lgdp", "eq"))) %>% 
    ggplot()+
    geom_ribbon(aes(ymin = L, ymax = U, x = Horizon), alpha = 0.5)+
    geom_line(aes(x = Horizon, y = M), color = "blue", size = 1)+
    geom_hline(yintercept = 0, color  = "red")+
    facet_wrap(facets = "Variable", ncol = pcol, scales = "free_y")+
    theme_minimal()
  
  return(cholirf)
}

