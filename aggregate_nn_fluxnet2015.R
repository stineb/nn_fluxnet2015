library( dplyr )
library( cgwtools )

source( "remove_outliers.R" )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
do.sites   = "FR-Pue"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
nrep       = 5
dotrain    = FALSE
overwrite_modis = FALSE
overwrite_mte = FALSE
##---------------------------------

siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

##------------------------------------------------
## Load overview L2 and initialise additional columns created here
##------------------------------------------------
load( "data/overview_data_fluxnet2015_L2.Rdata" )

##------------------------------------------------
## Get MTE-GPP for all sites
##------------------------------------------------
filn <- paste( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_8Days_4Beni.csv", sep="" )
if ( file.exists( filn ) ){
  mte_dl <- read.csv( paste( myhome, "/data/gpp_mte_rf_fluxnet_tramontana/GPP_Daily_4Beni.csv", sep = "" ), as.is=TRUE )
  mte_8d <- read.csv( filn, as.is=TRUE )
  avl_data_mte <- TRUE
} else {
  avl_data_mte <- FALSE
}

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

## XXX override
avl_data_modis <- FALSE

##------------------------------------------------
## Get PRI and CCI for all sites
##------------------------------------------------
filn <- paste( myhome, "data/pri_cci_fluxnet_from_marcos/data_euroflux_modis.Rdata" )

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

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 

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

    ## add row to aggregated data
    mysitename <- data.frame( mysitename=rep( sitename, nrow(nice) ) )

    if (nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar"){
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act * evi * ppfd, gpp_nn_pot = var_nn_pot * evi * ppfd )
    } else {
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act, gpp_nn_pot = var_nn_pot )
    }


  ##------------------------------------------------
  ## Add PRI and CCI from Marcos' dataset
  ##------------------------------------------------    
  if (avl_data_pri){

    if ( !is.element( sitename, missing_pri ) ) {
  
      ## add PRI and CCI data to 'nice'
      nice <- nice %>%  left_join( dplyr::filter( ddf_pri, site==sitename ), by=c( "year", "moy", "dom" ) ) %>%
                        mutate( 
                                # cci_norm = remove_outliers( cci / ( evi * ppfd ), coef=3.0 ),
                                # pri_norm = remove_outliers( pri / ( evi * ppfd ), coef=3.0 )
                                cci_norm = remove_outliers( cci / ( fpar ), coef=3.0 ),
                                pri_norm = remove_outliers( pri / ( fpar ), coef=3.0 )
                                )
      avl_pri <- TRUE
  
    } else {

      nice <- nice %>% mutate( ndvi   = rep( NA, nrow(.) ), 
                               evi_df_pri = rep( NA, nrow(.) ), 
                               cci    = rep( NA, nrow(.) ), 
                               pri    = rep( NA, nrow(.) ), 
                               scci   = rep( NA, nrow(.) ), 
                               spri   = rep( NA, nrow(.) ), 
                               cci_norm = rep( NA, nrow(.) ), 
                               pri_norm = rep( NA, nrow(.) ), 
                               ndsi   = rep( NA, nrow(.) ), 
                               wateri = rep( NA, nrow(.) )
                              )
      avl_pri <- FALSE

    }

    ##------------------------------------------------
    ## Add scaled PRI and CCI to within 0 and 1 by range
    ##------------------------------------------------
    range_cci <- range( nice$cci, na.rm=TRUE )
    range_pri <- range( nice$pri, na.rm=TRUE )

    nice <- nice %>% mutate(
                              scci = as.vector( scale( cci, center=range_cci[1], scale=( range_cci[2] - range_cci[1] ) ) ),
                              spri = as.vector( scale( pri, center=range_pri[1], scale=( range_pri[2] - range_pri[1] ) ) )
                            )

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
                "gpp_nn_act", 
                "gpp_nn_pot",
                "gpp_pmodel",
                "ppfd", 
                "fvar", 
                "bias_pmodel", 
                "ratio_obs_mod", 
                "soilm_mean", 
                "is_drought_byvar", 
                "evi", 
                "fpar", 
                "lue_obs_evi", 
                "lue_obs_fpar"
                )

  if (avl_data_pri) usecols <- c( usecols, 
                                  c(
                                    "ndvi", 
                                    "evi_df_pri", 
                                    "cci", 
                                    "pri", 
                                    "scci", 
                                    "spri", 
                                    "ndsi", 
                                    "wateri", 
                                    "pri_norm", 
                                    "cci_norm" 
                                    ) 
                                  )

  nice_agg <- rbind( 
                    nice_agg, 
                    cbind( 
                          mysitename, 
                          dplyr::select( nice, 
                                         one_of( usecols )
                                        )
                          ) 
                    )

  ##------------------------------------------------
  ## Get MTE-GPP for this site
  ##------------------------------------------------
  if (avl_data_mte){

    if (is.element( sitename, mte_8d$Site.code)){

      filn <- paste( "data/mte_", sitename, ".Rdata", sep="" )

      if ( file.exists(filn) && !overwrite_mte ){

        ## load 'nice_to_mte' from file
        load( filn )

      } else {

        ## prepare dataframe 'nice_to_mte'
        mte  <- dplyr::filter( mte_8d, Site.code==sitename )
        dmte <- dplyr::filter( mte_dl, Site.code==sitename )
        missing_mte <- FALSE

        mte$year_dec_start <- mte$StartYear + ( mte$StartDoY - 1 ) / 365 
        mte$year_dec_end   <- mte$EndYear   + ( mte$EndDoY   - 1 ) / 365 
        mte$year_dec       <- mte$StartYear + ( (mte$StartDoY + mte$EndDoY)/2 - 1 ) / 365 ## year_dec is the start of the period. approx( ..., method="linear") holds value [i] constant from year_dec[1] to year_dec[i+1]

        dmte$year_dec <- dmte$StartYear + ( dmte$StartDoY - 1 ) / 365 

        mte  <- dplyr::rename( mte,  mysitename=Site.code, gpp_mte=MTE, gpp_mte_m=MTE_M, gpp_mte_viterbo=MTE_Viterbo, gpp_rf=RF, gpp_rf_fromdaily=RF_from_daily )
        dmte <- dplyr::rename( dmte, mysitename=Site.code, doy=StartDoY, year=StartYear, gpp_rf=RFdaily )
        # mte <- dplyr::rename( mte, mysitename=Site.code, gpp_mte=MTE )

        ## replace with NA
        for (ivar in names(mte)){
          mte[[ ivar ]][ which( mte[[ ivar ]]==-9999 ) ] <- NA
        }
        for (ivar in names(dmte)){
          dmte[[ ivar ]][ which( dmte[[ ivar ]]==-9999 ) ] <- NA
        }


        ## Make 'nice' dataframe conform with 'mte'
        nice_to_mte <- c()
        for (idx in 1:nrow(mte)){
          sub <- dplyr::filter( nice, year_dec>=mte$year_dec_start[idx] & year_dec<=mte$year_dec_end[idx] )
          year_dec_save <- sub$year_dec[1]
          addline <- unlist( unname( apply( sub, 2, FUN=mean, na.rm=TRUE )))
          # addline[1] <- year_dec_save
          nice_to_mte <- rbind( nice_to_mte, addline )
        }

        nice_to_mte <- as.data.frame( nice_to_mte )
        colnames( nice_to_mte ) <- names( sub )
        rownames( nice_to_mte ) <- NULL

        mycolnames <- c( "gpp_mte", "gpp_mte_m", "gpp_mte_viterbo", "gpp_rf", "gpp_rf_fromdaily" )

        for (ivar in mycolnames){
          nice_to_mte[[ ivar ]] <- mte[[ ivar ]]
        }

        nice_to_mte$is_drought_byvar <- with( nice_to_mte, ifelse( is_drought_byvar<0.5, FALSE, TRUE ) )

        ## get bias measure: ( mod / obs )
        ## MTE
        nice_to_mte$bias_mte <- nice_to_mte$gpp_mte / nice_to_mte$gpp_obs
        nice_to_mte$bias_mte[ which( is.infinite( nice_to_mte$bias_mte ) ) ] <- NA

        nice_to_mte$ratio_obs_mod_mte  <-  nice_to_mte$gpp_obs / nice_to_mte$gpp_mte
        nice_to_mte$ratio_obs_mod_mte[ which(is.infinite(nice_to_mte$ratio_obs_mod_mte)) ] <- NA        

        ## RF
        nice_to_mte$bias_rf <- nice_to_mte$gpp_rf / nice_to_mte$gpp_obs
        nice_to_mte$bias_rf[ which( is.infinite( nice_to_mte$bias_rf ) ) ] <- NA

        nice_to_mte$ratio_obs_mod_rf  <-  nice_to_mte$gpp_obs / nice_to_mte$gpp_rf
        nice_to_mte$ratio_obs_mod_rf[ which(is.infinite(nice_to_mte$ratio_obs_mod_rf)) ] <- NA        


        ## add dmte to 'nice' dataframe
        nice$gpp_rf_daily <- rep( NA, nrow(nice) )
        for (idx in 1:nrow(nice)){
          idx_use <- which.min( abs( nice$year_dec[idx] - dmte$year_dec ) )
          nice$gpp_rf_daily[idx] <- dmte$gpp_rf[ idx_use ]

          ## add tolerance
          if ( (nice$year_dec[idx] - dmte$year_dec[idx_use]) > 1.5/365 ) nice$gpp_rf_daily[idx] <- NA
        }

        ## get bias measure: ( mod / obs )
        nice$bias_rf <- nice$gpp_rf_daily / nice$gpp_obs
        nice$bias_rf[ which( is.infinite( nice$bias_rf ) ) ] <- NA

        nice$ratio_obs_mod_dmte  <-  nice$gpp_obs / nice$gpp_rf_daily
        nice$ratio_obs_mod_dmte[ which(is.infinite(nice$ratio_obs_mod_dmte)) ] <- NA          

        ## save to file
        save( nice_to_mte, file=filn )

      }

      ## add row to aggregated data
      mysitename <- data.frame( mysitename=rep( sitename, nrow(nice_to_mte) ) )
      nice_to_mte_agg <- rbind( nice_to_mte_agg, cbind( mysitename, dplyr::select( nice_to_mte, bias_mte, is_drought_byvar, gpp_mte, gpp_obs ) ) )

      idxs_drought <- which( nice_to_mte$is_drought_byvar )
      overview[ which(overview$mysitename==sitename), ]$perr_mte <- 100 * sum( nice_to_mte$gpp_mte[idxs_drought] - nice_to_mte$gpp_obs[idxs_drought], na.rm=TRUE ) / sum( nice_to_mte$gpp_obs[idxs_drought], na.rm=TRUE )

    } else {

      missing_mte <- TRUE

    }

  }

  ##------------------------------------------------
  ## Get MODIS-GPP for this site
  ##------------------------------------------------
  if (avl_data_modis){

    filn <- paste( "data/modis_", sitename, ".Rdata", sep="" )

    if ( file.exists(filn) && !overwrite_modis ){

      ## load 'nice_to_modis' from file
      load( filn )
      avl_modisgpp <- TRUE

    } else {

      ## prepare 'nice_to_modis'
      modis <- try( read.csv( paste( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/", sitename, "/gpp_8d_modissubset_", sitename, ".csv", sep="" ), as.is=TRUE ))
      if (class(modis)!="try-error"){
        avl_modisgpp <- TRUE
        modis <- dplyr::rename( modis, gpp_modis=data )
        modis$gpp_modis <- modis$gpp_modis / 8

        ## Make 'nice' dataframe conform with 'modis'. Take mean of all variables across days for which year_dec is within four +/- 4 days of the modis date
        nice_to_modis <- c()
        for (idx in 1:nrow(modis)){
          sub <- dplyr::filter( nice, year_dec>=(modis$year_dec[idx] - 4/365 ) & year_dec<=(modis$year_dec[idx] + 4/365) )
          addline <- unlist( unname( apply( sub, 2, FUN=mean, na.rm=TRUE ) ) )
          nice_to_modis <- rbind( nice_to_modis, addline )
        }

        nice_to_modis <- as.data.frame(nice_to_modis)
        colnames( nice_to_modis ) <- names( sub )
        rownames( nice_to_modis ) <- NULL

        nice_to_modis$year_dec  <- modis$year_dec
        nice_to_modis$gpp_modis <- modis$gpp_modis

        nice_to_modis$is_drought_byvar <- with( nice_to_modis, ifelse( is_drought_byvar<0.5, FALSE, TRUE ) )

        # plot(  nice_to_modis$year_dec, nice_to_modis$gpp_obs, type='l' )
        # lines( nice_to_modis$year_dec, nice_to_modis$gpp_modis, col='red' )

        # plot( nice_to_modis$year_dec, (nice_to_modis$gpp_modis - nice_to_modis$gpp_obs) / nice_to_modis$gpp_obs, type='l' )

        ## get bias measure: log( mod / obs )
        nice_to_modis$bias_modis <- nice_to_modis$gpp_modis / nice_to_modis$gpp_obs
        nice_to_modis$bias_modis[ which(is.infinite(nice_to_modis$bias_modis)) ] <- NA

        nice_to_modis$ratio_obs_mod_modis  <-  nice_to_modis$gpp_obs / nice_to_modis$gpp_modis
        nice_to_modis$ratio_obs_mod_modis[ which(is.infinite(nice_to_modis$ratio_obs_mod_modis)) ] <- NA

        ## save to file
        save( nice_to_modis, file=filn )

      } else {
        avl_modisgpp <- FALSE
      }

    }

    ## add row to aggregated data
    if (avl_modisgpp){
      idxs_drought <- which( nice_to_modis$is_drought_byvar )
      overview[ which(overview$mysitename==sitename), ]$perr_modis <- 100 * sum( nice_to_modis$gpp_modis[idxs_drought] - nice_to_modis$gpp_obs[idxs_drought], na.rm=TRUE ) / sum( nice_to_modis$gpp_obs[idxs_drought], na.rm=TRUE )

      mysitename <- data.frame( mysitename=rep( sitename, nrow(nice_to_modis) ) )
      nice_to_modis_agg <- rbind( nice_to_modis_agg, cbind( mysitename, dplyr::select( nice_to_modis, bias_modis, is_drought_byvar, gpp_modis, gpp_obs ) ) )
    }

  }

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

