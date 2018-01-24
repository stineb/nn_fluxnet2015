library( dplyr )
library( cgwtools )

source( "remove_outliers.R" )
source( "get_spi_spei.R" )
source( "get_anom_bydoy.R" )
source( "movingAverage.R" )

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
testprofile = TRUE
##---------------------------------

siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

##------------------------------------------------
## Load overview L2 and initialise additional columns created here
##------------------------------------------------
load( "data/overview_data_fluxnet2015_L2.Rdata" )

##------------------------------------------------
## Initialise aggregated data
##------------------------------------------------
## fvar and soilm data to be complemented with cluster info
nice_agg <- data.frame()

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
nice_resh         <- data.frame()

## all possible soil moisture datasets
varnams_swc_full <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob", "soilm_obs" )

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  if (testprofile){
    infil <- paste( "data/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
    if (!file.exists(infil)){
      print(paste( "File", infil, "not found. Are they properly downloaded from bstocker.net and stored in ./data/fvar/ ?"))
    }
  } else {
    infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
  }

  if (verbose) print( paste( "opening file", infil ) )

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
  load( infil ) ## gets list 'nn_fluxnet'
  nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
  varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
  varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

  ## For SD-Dem, overwrite measured soil moisture - it's probably wrong as it does not fall below 0.25 or so
  if (sitename == "SD-Dem"){
    varnams_swc_mod <- varnams_swc[ !is.element( varnams_swc, varnams_swc_obs ) ]
    nice$soilm_mean <- apply( dplyr::select( nice, one_of(varnams_swc_mod)), 1, FUN=mean, na.rm=TRUE )
    nice$soilm_mean[ is.nan( nice$soilm_mean ) ] <- NA
  }

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
  mdf_spei <- get_spi_bysite( sitename, dplyr::filter( siteinfo, mysitename==sitename)$lat )
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
    
  ##------------------------------------------------
  ## Re-save data after additions to 'nice' dataframe
  ##------------------------------------------------
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
                "spei3"
                )

  # nice_agg <- dplyr::select( nice, one_of( usecols ) ) %>% cbind( mysitename, . ) %>% rbind( . )

  sub <- dplyr::select( nice, one_of( usecols ) )
  nice_agg <- rbind( nice_agg, cbind(  mysitename, sub ) )

  ##------------------------------------------------
  ## Reshape dataframe to stack data from different soil moisture datasets along rows
  ##------------------------------------------------
  for (isoilm in varnams_swc_full){
    if (isoilm=="soilm_obs") nice$soilm_obs <- nice$soilm_obs_mean
    if (is.null(nice[[ isoilm ]])){
      nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      nice[[ paste( "moist_",      isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      nice[[ isoilm ]]                                 <- rep( NA, nrow(nice) )
    }
    addrows <- data.frame( 
                            mysitename = sitename,
                            year_dec   = nice$year_dec,
                            var_nn_act = nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]],
                            var_nn_pot = nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]],
                            var_nn_vpd = nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]],
                            moist      = nice[[ paste( "moist_",      isoilm, sep="" ) ]],
                            lue_obs_evi= nice$lue_obs_evi,
                            soilm      = nice[[ isoilm ]],
                            vpd        = nice$vpd,
                            iabs       = nice$evi * nice$ppfd,
                            soilm_data = rep( isoilm, nrow(nice) )
                            )
    nice_resh <- rbind( nice_resh, addrows )
  }

}

print("... done.")

if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  save( nice_agg, file=paste( "data/nice_agg_", nam_target, char_fapar, ".Rdata", sep="") )
  save( nice_resh, file="data/nice_resh_lue_obs_evi.Rdata" )
  save( overview, file="data/overview_data_fluxnet2015_L3.Rdata" )               

} else {

  print("WARNING: NO SAVING AT THE END!")

}

