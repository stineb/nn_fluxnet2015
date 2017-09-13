##------------------------------------------------
## Get anomaly w.r.t. DOY
##------------------------------------------------
get_anom_bydoy <- function( nice ){

  require( dplyr )

  if (length(unique(nice$year))>3){

    ## aggregate to monthly values and calculate balance
    nice_mean <- nice %>% group_by( doy ) %>% 
                          summarise( soilm_mean_meanbydoy=mean( soilm_mean, na.rm=TRUE ) )

    nice <- nice %>%  left_join( nice_mean, by="doy" ) %>% 
                      mutate( soilm_mean_anom=soilm_mean-soilm_mean_meanbydoy )

  } else {

    nice$soilm_mean_anom <- rep( NA, nrow(nice) )    

  }


  return( nice$soilm_mean_anom )

}
