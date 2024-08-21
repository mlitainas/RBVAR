vardata_form = function(data , lags = 1, const = T, trend = F, trend_qua = F,
                        star = NULL, ex= NULL, ex_lag = NULL, dummy = NULL, xlags = 0, feedback_only = F ){
  # Check the consistency of the inputs
  # Test the data input
  # source("test.R")
  # data = health %>% select(y, lemp)
  # star = health %>% select(y_s, lemp_s)
  # lags = 2; reps = 300; burn = 50 ;const= T; trend =T; trend_qua = T; xlags = 2; i=1; lamda= 0.3;tau=0; epsilon=0.0001;dummy = NULL ; ex = NULL; feedback_only = F

  if ( names(data)[1] %in% c("Date", "date", "Time", "time") ) {
    print(paste("Time ID series identified:", names(data)[1], sep = " ") )
    TimeID =  data.frame( "TimeID" = data[,1] )
    data = data[,-1]
  } else {
    print("No Time ID series identified." )
    TimeID = NULL
  }

  if (!is.data.frame(data) & !is.tbl(data)) {
    stop("Please provide the argument 'data' either as a 'data.frame' or a 'tbl' object.")
  }

  if (any(is.na(data))) {
    stop("The argument data contains NAs. Function cannot handle mssing values.")
  }

  if (!is.numeric(as.matrix(data))) {
    stop("The argument data contains non numeric values. Please check you data.")
  }

  # Test the lags input
  if ( !is.numeric(lags)) {
    stop("Object lags must be a positive integer greater or equal to 1")
  }

  if (!is_scalar_atomic(lags) | !lags%%1 == 0 | lags < 1) {
    stop("Object lags must be a positive integer greater or equal to 1")
  }

  # Test the  inclusion of deterministic terms
  if( !is.logical(const) ) {
    stop("argument const takes a logical value. The default value is 'TRUE'")
  }

  if( !is.logical(trend)) {
    stop("argument trend takes a logical value. The default value is 'FALSE'")
  }

  if( !is.logical(trend_qua)) {
    stop("trend_qua takes a logical value. The default value is 'FALSE'")
  }

  # Test the ex input
  if (!is.data.frame(ex) & !is.tbl(ex) & !is.null(ex)) {
    stop("Please provide the argument 'ex' either as a 'data.frame' or a 'tbl' object or a matrix order. The default value is 'NULL'")
  }

  if (!is.null(ex)) {
    if (!is.numeric(as.matrix(ex))) {
      stop("The argument 'ex' does not include numeric values")
    }

    if (any(is.na(ex))) {
      stop("The argument 'ex' contains NAs. Function cannot handle mssing values.")
    }
  }

  # Test the ex_lag input
  #if (!is.data.frame(ex_lag) & !is.tbl(ex_lag) & !is.null(ex_lag)) {
  # stop("Please provide the argument 'ex_lag' either as a 'data.frame' or a 'tbl' object or a matrix order.")
  #}

  # if (!is.null(ex_lag)) {
  #   if (!is.numeric(as.matrix(ex_lag))) {
  #     stop("The argument 'ex_lag' does not include numeric values")
  #   }

  #  if (any(is.na(ex_lag))) {
  #    stop("The argument 'ex_lag' contains NAs. Function cannot handle mssing values.")
  #   }
  #}

  # Test the compatibility of the objects
  if (!is.null(ex)){
    if ( !dim(data)[1] == dim(ex)[1] ){
      stop("Arguments 'data' and 'ex' are not of the same langht")
    }
  }

  #if (!is.null(ex_lag)){
  # if ( !dim(data)[1] == dim(ex_lag)[1] ){
  #  stop("Arument 'data' and 'ex_lag' are not of the same langht")
  #}
  #}

  # Test the xlags input
  if ( !is.numeric(xlags) | !is_scalar_atomic(xlags) | !xlags%%1 == 0 | xlags < 0 ) {
    stop("Object xlags must be a positive integer greater or equal to 1")
  }
  if (xlags > lags) {
    stop("Argument 'xlags' must be lower than argument 'lag'.")
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Start data preparation #
  number_of_deterministic   <- 0 # Sum the number of the deterministic terms
  number_of_exogenous       <- 0 # Sum the number of the exogenous variables
  number_of_exogenous_lags  <- 0
  number_of_star_variables  <- ncol(star)  # number of star variables
  number_of_lagged_star_var <- ncol(star)*xlags  # number of lagged star variables
  number_of_variables       <- ncol(data)  # number of variables
  number_of_observations    <- nrow(data)  # number of observations

  #Names
  names_of_endog_variables  <- names(data)
  names_of_star_variables   <- names(star)



  # separate exogenous from endogenous variables and calculate the lags
  y = data

  # calculate lags for endogenous
  x  <-  dplyr::lag(y, 1)

  names(x) <- paste(names(x), 1, sep = "_")

  if ( lags >1 ) {

    for (i in 2:lags) {
      lag_temp <- dplyr::lag(y, i)
      names(lag_temp) <- paste(names(lag_temp),i, sep = "_")
      x <- cbind(x,  lag_temp)
    }
  }

  # Deterministic terms

  # Initialize the block of the deterministc terms
  deterministic_block <- NULL
  # check for constant
  if ( const == T ) {
    deterministic_block <- cbind(c = rep(1, number_of_observations)) # First column filled with 1s
    number_of_deterministic <- number_of_deterministic + 1
  }

  # check for trend
  if ( trend == T ) {
    deterministic_block <- cbind(deterministic_block, "trend" = 1:nrow(y)) # Second column filled with linear trend
    number_of_deterministic <- number_of_deterministic + 1
  }

  # check for quadratic trend
  if ( trend_qua == T ) {
    deterministic_block     <-  cbind(deterministic_block, "q_trend"  = c(1:nrow(y))^2) # third column filled with quadratic trend
    number_of_deterministic <-  number_of_deterministic + 1
  }

  # Initialize the block of the dummies
  #number_of_dummy           <- NULL
  if( !is.null(dummy) ) {
    deterministic_block     <-  cbind(deterministic_block, dummy) # third column filled with quadratic trend
    number_of_deterministic <-  number_of_deterministic + 1
  }

  star_block <-  NULL
  # Contemporaneous values of star variables
  if( !is.null(star) ) {
    number_of_exogenous = number_of_exogenous + ncol(star)
    temp_star = star
    names(temp_star) = paste0(names_of_star_variables, "_0")

    # Lagged values of star variables
    if( xlags > 0 ) {

      temp_star_lag        <- dplyr::lag(star, 1)
      names(temp_star_lag) <- paste0(names_of_star_variables, "_1")

      if ( xlags > 1 ) {

        for (i in 2:xlags) {

          xlag_temp <- dplyr::lag(star, i)
          names(xlag_temp) <- str_c(names_of_star_variables,i, sep = "_")
          temp_star_lag    <- cbind(temp_star_lag,  xlag_temp)

        }
      }
    }



    # bring contemporaneous and star variables in the same block
    if (feedback_only){
      star_block <- cbind(temp_star_lag)
    } else{
      star_block <- cbind(temp_star, temp_star_lag)
    }
  }

  # bind the exogenous variables with the lagged matrix of endogenous
  ncol_star_block = 0
  if ( !is.null(star_block) ) {
    x = cbind(star_block, x)
    ncol_star_block = ncol(star_block)
  }
  ncol_ex = 0
  if ( !is.null(ex) ) {
    x = cbind(ex, x)
    ncol_ex = ncol(ex)
  }
  ncol_det = 0
  if ( !is.null(deterministic_block) ) {
    x = cbind(deterministic_block, x)
    ncol_det = ncol(deterministic_block)
  }

  # Total number of exogenous variables
  number_of_exogenous = ncol_det + ncol_ex + ncol_star_block

  # I turn the data into matrices to facilitate the step of the calculations
  y_lhs =  as.matrix(y[(lags+1):number_of_observations,]) # LHS
  x_rhs =  as.matrix(x[(lags+1):number_of_observations,]) # RHS

  # Degrees of freedom
  dof  = nrow(data) - ncol(x_rhs)

  # Handle the date series
  if (!is.null(TimeID)) {
    TimeID = TimeID[(lags+1):number_of_observations,1, drop=F]
  }

  # Function will return a list with the VAR data and a number of parameters

  # that will be used in the calculations
  vardata = list(
    "y_lhs"                     = y_lhs,
    "x_lhs"                     = x_rhs,
    "number_of_endogenous"      = number_of_variables,
    "number_of_exogenous"       = number_of_exogenous,
    "number_of_deterministic"   = number_of_deterministic,
    "number_of_observations"    = number_of_observations,
    "number_of_obs_lags"        = nrow(y_lhs),
    "number_of_lags"            = lags,
    "number_of_xlags"           = xlags,
    "dof"                       = dof,
    "TimeID"                    = TimeID,
    # names
    "names_of_endog_variables"  = names_of_endog_variables,
    "names_of_star_variables"   = names_of_star_variables

  )


  return(vardata)
}
