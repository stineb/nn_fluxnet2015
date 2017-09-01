library( dplyr )
library( cgwtools )

source( "remove_outliers.R" )
source( "get_spi_spei.R" )
source( "get_anom_bydoy.R" )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
# do.sites   = do.sites[1:5]
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
overwrite_modis = FALSE
overwrite_mte = FALSE
verbose    = FALSE
##---------------------------------

siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

##------------------------------------------------
## Load overview L2 and initialise additional columns created here
##------------------------------------------------
load( "data/overview_data_fluxnet2015_L2.Rdata" )


## XXX override
avl_data_mte <- FALSE

##------------------------------------------------
## Check availability of MODIS GPP data
##------------------------------------------------
fillist <- list.files( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/" )
if (length(fillist)==0) {
  avl_data_modis <- FALSE
} else {
  avl_data_modis <- TRUE
}

##------------------------------------------------
## Get PRI and CCI for all sites
##------------------------------------------------
filn <- paste( myhome, "data/pri_cci_fluxnet_from_marcos/data_euroflux_modis.Rdata", sep="" )

if ( file.exists(filn) ){

  avl_data_pri <- TRUE
  load( filn )  # loads a dataframe called 'data'
  df_pri <- data  # rename the variable to df_pri
  rm( data )
  df_pri <- df_pri %>%  rename( year=year.modis, moy=month.modis, dom=day.modis, doy_df_pri=julday.modis ) 

  ## Mistake in Marcos' calculation of PRI (before using bands 11 and 3, but it should be calculated using bands 11 and 4)
  ## CCI is calculated as cci = ( BRF_B11 - BRF_B1 ) / ( BRF_B11 + BRF_B1 )
  df_pri <- df_pri %>% mutate( pri = ( BRF_B11 - BRF_B4 ) / ( BRF_B11 + BRF_B4 ) ) # iolanda: 13 instead of 4???

  ## aggregate to daily mean values
  ddf_pri <- df_pri %>% mutate( 
                                year       = as.integer( as.character( year       )), 
                                moy        = as.integer( as.character( moy        )), 
                                dom        = as.integer( as.character( dom        )), 
                                doy_df_pri = as.integer( as.character( doy_df_pri )),  
                                localh     = as.numeric( as.character( localh     ))
                                ) %>%
                        group_by( site, year, moy, dom ) %>%
                        summarise( 
                                    ndvi       = mean( ndvi,    na.rm=TRUE ),
                                    evi_df_pri = mean( evi ,    na.rm=TRUE ),
                                    cci        = mean( cci ,    na.rm=TRUE ),
                                    pri        = mean( pri ,    na.rm=TRUE ),
                                    ndsi       = mean( ndsi,    na.rm=TRUE ),
                                    wateri     = mean( wateri , na.rm=TRUE )
                                  ) %>%
                        mutate( 
                                cci = remove_outliers( cci, coef=3.0 ),
                                pri = remove_outliers( pri, coef=3.0 )
                              )

  ## get list of sites for which we don't have any PRI data
  sitelist_pri <- as.character( unique( ddf_pri$site ) )
  sitelist_my  <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename
  missing_pri  <- sitelist_my[ which( !is.element( sitelist_my, sitelist_pri ) ) ]

} else {

  avl_data_pri <- FALSE

}


##------------------------------------------------
## Initialise aggregated data
##------------------------------------------------
## fvar and soilm data to be complemented with cluster info
nice_agg          <- data.frame()
nice_to_mte_agg   <- data.frame()
nice_to_modis_agg <- data.frame()


## check and override if necessary
if ( nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
  plotlue <- TRUE
  if (nam_target=="lue_obs_evi"){
    fapar_data <- "evi"
  } else if (nam_target=="lue_obs_fpar"){
    fapar_data <- "fpar"
  }
  if (use_fapar){
    print("WARNING: setting use_fapar to FALSE")
    use_fapar <- FALSE
  }
}

## identifier for output files
if (use_fapar){
  if (nam_target=="lue_obs_evi"){
    char_fapar <- "_withEVI"
  } else if (nam_target=="lue_obs_fpar"){
    char_fapar <- "_withFPAR"
  } else {
    print("ERROR: PROVIDE VALID FAPAR DATA!")
  }
} else {
  char_fapar <- ""
}

if (use_weights){
  char_wgt <- "_wgt"
} else {
  char_wgt <- ""
}

print( "Aggregating and complementing data for all sites ..." )

##------------------------------------------------
## Initialise aggregated data
##------------------------------------------------
## fvar and soilm data to be complemented with cluster info
nice_agg          <- data.frame()
nice_to_mte_agg   <- data.frame()
nice_to_modis_agg <- data.frame()
nice_resh         <- data.frame()

## all possible soil moisture datasets
varnams_swc_full <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob", "soilm_obs" )

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 
  if (verbose) print( paste( "opening file", infil ) )

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

    ## remove columns again if they had bee added already before
    if (avl_data_pri){
      if ( is.element( "pri", names( nice ) ) ) nice$pri <- NULL
      if ( is.element( "cci", names( nice ) ) ) nice$cci <- NULL
      if ( is.element( "spri", names( nice ) ) ) nice$spri <- NULL
      if ( is.element( "scci", names( nice ) ) ) nice$scci <- NULL
      if ( is.element( "pri_norm", names( nice ) ) ) nice$pri_norm <- NULL
      if ( is.element( "cci_norm", names( nice ) ) ) nice$cci_norm <- NULL
      if ( is.element( "ndvi", names( nice ) ) ) nice$ndvi <- NULL
      if ( is.element( "ndsi", names( nice ) ) ) nice$ndsi <- NULL
      if ( is.element( "wateri", names( nice ) ) ) nice$wateri <- NULL
      if ( is.element( "evi_df_pri", names( nice ) ) ) nice$evi_df_pri <- NULL
      if ( is.element( "site", names( nice ) ) ) nice$site <- NULL
    }

    ## For SD-Dem, overwrite measured soil moisture - it's probably wrong as it does not fall below 0.25 or so
    if (sitename == "SD-Dem"){
      varnams_swc_mod <- varnams_swc[ !is.element( varnams_swc, varnams_swc_obs ) ]
      nice$soilm_mean <- apply( dplyr::select( nice, one_of(varnams_swc_mod)), 1, FUN=mean, na.rm=TRUE )
      nice$soilm_mean[ is.nan( nice$soilm_mean ) ] <- NA
    }

    nice$bias_pmodel  <-  nice$gpp_pmodel / nice$gpp_obs
    nice$bias_pmodel[ which(is.infinite(nice$bias_pmodel)) ] <- NA

    nice$ratio_obs_mod  <-  nice$gpp_obs / nice$gpp_pmodel
    nice$ratio_obs_mod[ which(is.infinite(nice$ratio_obs_mod)) ] <- NA

    ## Add information of EVI extremes and moving average fvar (preceeding 365 days)
    out_evianomalies <- nn_fluxnet[[ sitename ]]$out_evianomalies; fapar_extremes <- out_evianomalies$extremes
    nice$is_fapar_extreme <- rep( FALSE, nrow(nice) )
    for (iinst in nrow(fapar_extremes)){
      nice$is_fapar_extreme[ fapar_extremes$idx_start[iinst]:(fapar_extremes$idx_start[iinst]+fapar_extremes$len[iinst]-1) ] <- TRUE
    }
    tmp <- approx( nice$year_dec, nice$fvar, xout=nice$year_dec )$y
    nice$fvar_rollmean <- movingAverage( tmp, 365, centered=FALSE )

    ## Add SPI and SPEI data
    load( "./data/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata" )
    mdf_spei <- get_spi_bysite( sitename )
    nice <- nice %>%  dplyr::select( -contains("spi1"), -contains("spi3"), -contains("spei1"), -contains("spei3") ) %>% 
                      left_join( dplyr::select( mdf_spei, year, moy, spi1, spi3, spei1, spei3 ), by=c("year", "moy") )

    ## add soil moisture anomaly
    nice$soilm_mean_anom <- get_anom_bydoy( nice )

    ## add row to aggregated data
    mysitename <- data.frame( mysitename=rep( sitename, nrow(nice) ) )

    if (nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar"){
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act * evi * ppfd, gpp_nn_pot = var_nn_pot * evi * ppfd, gpp_nn_vpd = var_nn_vpd * evi * ppfd )
    } else {
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act, gpp_nn_pot = var_nn_pot, gpp_nn_vpd = var_nn_vpd )
    }

    ## fill with NA if respective soil moisture data was not used
    for (isoilm in varnams_swc_full){
      if ( !is.element( paste( "var_nn_act_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "var_nn_pot_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "var_nn_vpd_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "moist_",      isoilm, sep="" ), names(nice) ) ) nice[[ paste( "moist_",      isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
    }
    
  }

  ##------------------------------------------------
  ## Re-save data after additions to 'nice' dataframe
  ##------------------------------------------------
  # if (verbose) print( paste( "resaving nice into file", infil ) )
  # if (verbose) print( "names:" ); print( names(nice) )
  nn_fluxnet[[ sitename ]]$nice <- nice
  resave( nn_fluxnet, file=infil )

  ##------------------------------------------------
  ## record for aggregated data
  ##------------------------------------------------
  usecols <- c(
                "year_dec",
                "year",
                "doy",
                "gpp_obs", 
                "var_nn_pot", 
                "var_nn_act",
                "var_nn_vpd",
                "gpp_nn_act", 
                "gpp_nn_pot",
                "gpp_pmodel",
                "ppfd", 
                "fvar",
                "fvar_rollmean", 
                "bias_pmodel", 
                "ratio_obs_mod", 
                "soilm_mean",
                "vpd",
                "is_drought_byvar", 
                "is_fapar_extreme",
                "evi", 
                "fpar", 
                "lue_obs_evi", 
                "lue_obs_fpar",
                paste( "var_nn_act_", varnams_swc_full, sep="" ),
                paste( "var_nn_pot_", varnams_swc_full, sep="" ),
                paste( "var_nn_vpd_", varnams_swc_full, sep="" ),
                paste( "moist_", varnams_swc_full, sep="" ),
                "soilm_mean_anom",
                "aet_pmodel",
                "pet_pmodel",
                "spi1",
                "spi3",
                "spei1",
                "spei3",
                )

  # nice_agg <- dplyr::select( nice, one_of( usecols ) ) %>% cbind( mysitename, . ) %>% rbind( . )

  sub <- dplyr::select( nice, one_of( usecols ) )
  nice_agg <- rbind( nice_agg, cbind(  mysitename, sub ) )

  # ##------------------------------------------------
  # ## Reshape dataframe to stack data from different soil moisture datasets along rows
  # ##------------------------------------------------
  # for (isoilm in varnams_swc_full){
  #   if (isoilm=="soilm_obs") nice$soilm_obs <- nice$soilm_obs_mean
  #   if (is.null(nice[[ isoilm ]])){
  #     nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "moist_",      isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ isoilm ]]                                 <- rep( NA, nrow(nice) )
  #   }
  #   addrows <- data.frame( 
  #                           mysitename = sitename,
  #                           year_dec   = nice$year_dec,
  #                           var_nn_act = nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]],
  #                           var_nn_pot = nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]],
  #                           var_nn_vpd = nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]],
  #                           moist      = nice[[ paste( "moist_",      isoilm, sep="" ) ]],
  #                           lue_obs_evi= nice$lue_obs_evi,
  #                           soilm      = nice[[ isoilm ]],
  #                           vpd        = nice$vpd,
  #                           iabs       = nice$evi * nice$ppfd,
  #                           soilm_data = rep( isoilm, nrow(nice) )
  #                           )
  #   nice_resh <- rbind( nice_resh, addrows )
  # }


}

print("... done.")

if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  save( nice_agg, file=paste( "data/nice_agg_", nam_target, char_fapar, ".Rdata", sep="") )

  if (avl_data_mte)   save( nice_to_mte_agg,   file=paste("data/nice_mte_agg_",   nam_target, char_fapar, ".Rdata", sep="") )
  if (avl_data_modis) save( nice_to_modis_agg, file=paste("data/nice_modis_agg_", nam_target, char_fapar, ".Rdata", sep="") )

  if (avl_data_pri) save( missing_pri, file=paste("data/missing_pri_", nam_target, char_fapar, ".Rdata", sep="") )

  save( overview, file="data/overview_data_fluxnet2015_L3.Rdata" )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

