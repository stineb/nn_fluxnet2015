niceify <- function( data, data_full, method="exact" ){
  ##------------------------------------------------
  ## Expand data to get all time steps again and fill
  ## missing data with NAs
  ##------------------------------------------------

  ## initialise
  nice <- subset( data_full, select=c( year_dec ) )

  ## get (actual) variables (excluding time variables)
  vars <- names( data )
  vars <- vars[ -which( is.element( vars, names(nice) ) ) ]

  ## expand and initialise 'nice' with NAs
  for (ivar in vars){
    nice[[ ivar ]] <- rep( NA, nrow(data_full) )
  }

  ## fill up data
  for (idx in 1:nrow(data)){

    if (method=="exact"){
      putidx <- which( nice$year_dec==data$year_dec[ idx ] )
    } else if (method=="closest"){
      putidx <- which.min( abs( nice$year_dec-data$year_dec[ idx ]) )      
    }
    for (ivar in vars){
      nice[[ ivar ]][ putidx ] <- data[[ ivar ]][ idx ]
    }

  }

  return( nice )

}
