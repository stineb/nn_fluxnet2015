get_evianomalies_bysite <- function( sitename, quantile, do.plot=TRUE ){

  require(dplyr)
  require(stats)

  source( "get_consecutive.R" )

  # ## debug
  # sitename="FR-Pue"
  # quantile=0.03
  # do.plot=TRUE
  # #################

  df <- read.csv( paste( "sofun/input_fluxnet2015_sofun/evi_modissubsets/evi_modissubset_", sitename, ".csv", sep="" ) )

  df$year_dec <- df$yr + (df$doy-1)/365

  ## get mean seasonal cycle
  df_meanseason <- aggregate( evi ~ doy, data=df, mean )

  ## get EVI anomalies
  df$anom <- rep( NA, nrow(df) )
  for (idx in 1:nrow(df_meanseason)){
    idxs <- which( df$doy==df_meanseason$doy[idx] )
    df$anom[ idxs ] <- df$evi[ idxs ] - df_meanseason$evi[idx]
  }


  ## plot EVI
  years <- unique( df$yr )
  pdf( paste( "fig_fapar_fluxnet/evi_season_", sitename, ".pdf", sep="" ), width=7, heigh=6 )
  par(las=1)
  plot( c(1,365), c(0,1), type="n", xlab="DOY", ylab="EVI" )
  for (iyr in years){
    sub <- dplyr::filter( df, yr==iyr )
    lines( sub$doy, sub$evi, col=rgb(0,0,0,0.3) )
  }
  lines( df_meanseason$doy, df_meanseason$evi, col='red', lwd=2 )
  title( sitename )
  dev.off()

  # ## plot EVI anomalies
  # par(las=1)
  # plot( c(1,365), c(-0.3,0.3), type="n", xlab="DOY", ylab="EVI anomaly" )
  # for (iyr in years){
  #   sub <- dplyr::filter( df, yr==iyr )
  #   lines( sub$doy, sub$anom, col=rgb(0,0,0,0.3) )
  # }

  ## Define whether is drought
  isneg <- ( df$anom < 0 )

  ## replace NA with FALSE
  isneg[ which(is.na(isneg)) ] <- FALSE

  ## get dry periods and attach to list
  instances <- get_consecutive( isneg, anom=df$anom, leng_threshold=2, do_merge=FALSE )

  if (nrow(instances)>0){

    ## define upper 5% quantile and identify extreme instances
    threshold <- quantile( instances$deficit, probs=quantile, na.rm=TRUE )
    # hist( instances$deficit, main=sitename )
    # abline( v=threshold, col='red' )
    extremes <- dplyr::filter( instances, deficit<threshold )

    ## add info (year_dec) to 'extremes' and 'instances' data frames
    instances$year_dec_start <- rep( NA, nrow(instances) )
    instances$year_dec_end   <- rep( NA, nrow(instances) )
    for (idx in 1:nrow(instances)){
      instances$year_dec_start[idx] <- df$year_dec[ instances$idx_start[idx] ]
      instances$year_dec_end[idx]   <- df$year_dec[ (instances$idx_start[idx]+instances$len[idx]-1) ]
    }

    extremes$year_dec_start <- rep( NA, nrow(extremes) )
    extremes$year_dec_end   <- rep( NA, nrow(extremes) )
    for (idx in 1:nrow(extremes)){
      extremes$year_dec_start[idx] <- df$year_dec[ extremes$idx_start[idx] ]
      extremes$year_dec_end[idx]   <- df$year_dec[ (extremes$idx_start[idx]+extremes$len[idx]-1) ]
    }
  
  } else {
    extremes <- instances
  }

  return( list( instances=instances, extremes=extremes ) )

}



