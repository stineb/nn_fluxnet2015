cleandata_nn <- function( data, varnam ){
  ##------------------------------------------------
  ## Remove cold days and days where GPP is negative
  ##------------------------------------------------
  require( dplyr )

  if (varnam=="gpp_obs"){
    data <- dplyr::filter( data, !is.na(gpp_obs) )
    data <- dplyr::filter( data, gpp_obs > 0.0 )

  } else if (varnam=="et_obs"){
    data <- dplyr::filter( data, !is.na( et_obs ) )    

  } else if (varnam=="wue_obs"){
    data <- dplyr::filter( data, !is.na( wue_obs ) )    

  } else if (varnam=="lue_obs"){
    data <- dplyr::filter( data, !is.na( lue_obs ) )    

  } else if (varnam=="lue_obs_evi"){
    data <- dplyr::filter( data, !is.na( lue_obs_evi ) )    

  } else if (varnam=="lue_obs_fpar"){
    data <- dplyr::filter( data, !is.na( lue_obs_fpar ) )    

  }

  return( data )
}
