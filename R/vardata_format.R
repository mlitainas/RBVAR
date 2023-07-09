vardata_form = function(data , lags = 1, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL,dummy = NULL, xlags = 0 ){
  # data = df %>% select(NT13, lgs_pc, lgdp_pc);  lags = 2; const = T; trend = F;trend_qua = F; ex = df$TOTAL[,6]; ex_lag = df$TOTAL[,6:7]; xlags = 1
  # Check the consistency of the inputs
  # Test the data input

  if ( names(data)[1] %in% c("Date", "date", "Time", "time") ) {
    print(paste("Time ID series identified:", names(data)[1], sep = " ") )
    TimeID =  data.frame( "TimeID" = data[,1] )
    data = data[,-1]
  } else {
    print("No Time ID series identified." )
    TimeID =  data.frame( "TimeID" = 1:(nrow(data)-lags) )
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
  if (!is.data.frame(ex_lag) & !is.tbl(ex_lag) & !is.null(ex_lag)) {
    stop("Please provide the argument 'ex_lag' either as a 'data.frame' or a 'tbl' object or a matrix order.")
  }

  if (!is.null(ex_lag)) {
    if (!is.numeric(as.matrix(ex_lag))) {
      stop("The argument 'ex_lag' does not include numeric values")
    }

    if (any(is.na(ex_lag))) {
      stop("The argument 'ex_lag' contains NAs. Function cannot handle mssing values.")
    }
  }

  # Test the compatibility of the objects
  if (!is.null(ex)){
    if ( !dim(data)[1] == dim(ex)[1] ){
      stop("Arguments 'data' and 'ex' are not of the same langht")
    }
  }

  if (!is.null(ex_lag)){
    if ( !dim(data)[1] == dim(ex_lag)[1] ){
      stop("Arument 'data' and 'ex_lag' are not of the same langht")
    }
  }

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
  number_of_variables       <- ncol(data)  # number of variables
  number_of_observations    <- nrow(data)  # number of observations
  names_of_exogenous_lagged <- NULL
  temp_ex_lag <- NULL
  # Initialize the block of the deterministc terms
  deterministic_block       <- NULL
  # Initialize the block of the dummies
  number_of_dummy           <- NULL

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

    deterministic_block = cbind(deterministic_block, "q_trend"  = c(1:nrow(y))^2) # third column filled with quadratic trend

    number_of_deterministic <- number_of_deterministic + 1

  }

  if( !is.null(dummy) ) {
    number_of_exogenous = number_of_exogenous + ncol(dummy)
  }


  if( !is.null(ex) ) {
    number_of_exogenous = number_of_exogenous + ncol(ex)
    names(ex) = paste0(names(ex), "_0")
  }

  # calculate lags for the exogenous
  if( !is.null(ex_lag) ) {

    names_of_exogenous_lagged = names(ex_lag)
    temp_ex_lag = dplyr::lag(ex_lag, 1)
    names(temp_ex_lag) = paste0(names_of_exogenous_lagged, "_1")


    if ( xlags > 1 ) {

      for (i in 2:xlags) {

        xlag_temp = dplyr::lag(ex_lag, i)

        names(xlag_temp) = str_c(names(xlag_temp),i, sep = "_")

        temp_ex_lag = cbind(temp_ex_lag,  xlag_temp)

      }


    }

    number_of_exogenous = number_of_exogenous + ncol(temp_ex_lag)
    # bind the exogenous variables with the lagged matrix of endogenous

    if ( !is.null(ex) ) {
      x = cbind(ex,temp_ex_lag,x)
    } else {
      x = cbind(temp_ex_lag,x)
    }
  } else {
    if (!is.null(ex)) {
      x = cbind(ex,x)
    } else {
      x = cbind(x)
    }
  }


  if( !is.null(deterministic_block) ){

    if (!is.null(dummy)) {

      x = cbind(deterministic_block, dummy, x)

    } else {

      x = cbind(deterministic_block, x)

    }
    } else {

      if( !is.null(dummy) ) {

          x = cbind(dummy,x)

        }

  }


  # I turn the data into matrices to facilitate the step of the calculations
  y_lhs =  as.matrix(y[(lags+1):number_of_observations,]) # LHS
  x_rhs =  as.matrix(x[(lags+1):number_of_observations,]) # RHS

  # Total number of exogenous variables
  number_of_exogenous = number_of_deterministic + number_of_exogenous

  #
  names_endogenous  = names(y)

  # Degrees of freedom
  dof  = nrow(data) - ncol(x_rhs)

  # Handle the date series
  TimeID = TimeID[(lags+1):number_of_observations,1, drop=F]

  # Function will return a list with the VAR data and a number of parameters
  # that will be used in the calculations
  vardata = list(
    "y_lhs" = y_lhs,
    "x_lhs" = x_rhs,
    "number_of_variables"       = number_of_variables,
    "number_of_exogenous"       = number_of_exogenous,
    "number_of_deterministic"   = number_of_deterministic,
    "number_of_observations"    = number_of_observations,
    "number_of_obs_lags"        = nrow(y_lhs),
    "number_of_lags"            = lags,
    "number_of_xlags"           = xlags,
    "dof"                       = dof,
    "TimeID"                    = TimeID,
    # names
    "names_of_endogenous"       = names_endogenous,
    "names_of_exogenous_lagged" = names_of_exogenous_lagged

    )


  return(vardata)
}
