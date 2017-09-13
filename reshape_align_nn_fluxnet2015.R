reshape_align_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", bysm=FALSE, use_fapar=FALSE, use_weights=FALSE, overwrite=TRUE, verbose=FALSE ){

  # ## debug-------------------
  # sitename = "FR-Pue"
  # nam_target="lue_obs_evi"
  # bysm=FALSE
  # use_fapar=FALSE
  # use_weights=FALSE
  # overwrite=TRUE
  # verbose=TRUE
  # #--------------------------

  require( dplyr )
  require( tidyr )
  require( cgwtools )

  source( "get_consecutive.R" ) 

  ## check and override if necessary
  if ( nam_target=="lue_obs" || nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
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

  if (bysm){
    char_bysm <- "_bysm"
  } else {
    char_bysm <- ""
  }

  before <- 30
  after  <- 100

  ## Bins for different variables
  fvarbins  <- seq( from=-20, to=40, by=20 )
  faparbins <- seq( from=-20, to=40, by=20 )
  iwuebins  <- seq( from=-30, to=60, by=30 )

  bincentres_fvar  <- fvarbins[1:(length(fvarbins)-1)]   + (fvarbins[2]-fvarbins[1])/2
  bincentres_fapar <- faparbins[1:(length(faparbins)-1)] + (faparbins[2]-faparbins[1])/2
  bincentres_iwue  <- iwuebins[1:(length(iwuebins)-1)]   + (iwuebins[2]-iwuebins[1])/2

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
  if (verbose) print("loading nn_fVAR file ...")
  dir <- "./data/fvar/"
  infil <- paste( dir, "nn_fluxnet2015_", sitename, "_", nam_target, char_fapar, ".Rdata", sep="" ) 
  if (verbose) print( paste( "reading file", infil ) )
  load( infil ) ## gets list 'nn_fluxnet'
  df <- as.data.frame( nn_fluxnet[[ sitename ]]$nice ) %>% dplyr::select( year_dec, gpp_obs, var_nn_pot, var_nn_act, temp, ppfd, fvar, soilm_mean, vpd, evi, fpar, wue_obs, is_drought_byvar, gpp_pmodel, gpp_obs_gfd, iwue )

  droughts <- nn_fluxnet[[ sitename ]]$droughts        

  names_alg <- c( names(df), "dday")

  if (nrow(droughts)>1){
    ##--------------------------------------------------------
    ## re-arrange data, aligning by beginning of drought
    ##--------------------------------------------------------
    filn <- paste( "data/aligned_", sitename, char_bysm, ".Rdata", sep="" )
    if (!file.exists(filn)||overwrite){
      if (verbose) print("aligning df ...")
      # after  <- min( max(droughts$len), 200 )
      # before <- floor( after/2.0 )
      data_alg_dry  <- array( NA, dim=c( before+after+1, ncol(df)+1, nrow(droughts) ) )

      for ( iinst in 1:nrow(droughts) ){
        data_alg_dry[,(ncol(df)+1),iinst] <- (-before:after)+1 ## calling this 'dday' = drought day
        for (idx in -before:min( (droughts$len[iinst]-1), after) ){
          if ( (droughts$idx_start[iinst]+idx)>0 ){
            for (icol in 1:ncol(df)){
              data_alg_dry[ idx+before+1, icol, iinst ] <- df[ droughts$idx_start[iinst]+idx, icol ]
            }
          }
        }
        ## remove data after drought onset that is no longer classified as drought
        dropidxs <- which( data_alg_dry[,which(names_alg=="is_drought_byvar"),iinst]==1 & data_alg_dry[,which(names_alg=="dday"),iinst]<0 )
        data_alg_dry[ dropidxs,,iinst ] <- NA
      }
      save( data_alg_dry, names_alg, fvarbins, faparbins, iwuebins, before, after, bincentres_fvar, bincentres_fapar, bincentres_iwue, file=filn )
    
    } else {
      
      load( filn )
    
    }
    if (verbose) print("done")

    ##--------------------------------------------------------
    ## Bin aligned data and expand from 3D array to dataframe
    ##--------------------------------------------------------
    ## expand 'data_alg_dry' to get a data frame that now has 'dday' in it
    df_dday <- data.frame()
    for (iinst in seq(dim(data_alg_dry)[3])){
      add <- as.data.frame( data_alg_dry[,,iinst] )
      colnames(add) <- names_alg
      add$mysitename <- rep( sitename, dim( data_alg_dry[,,iinst] )[1])
      add$inst <- rep( iinst, dim( data_alg_dry[,,iinst] )[1] )
      df_dday <- rbind( df_dday, add )
    }
    # df_dday <- df_dday[ !is.na(df_dday$year_dec), ]

    ## add bin information based on dday to expanded df
    df_dday <- df_dday %>% mutate( infvarbin  = cut( as.numeric(dday), breaks = fvarbins ) )
    df_dday <- df_dday %>% mutate( infaparbin = cut( as.numeric(dday), breaks = faparbins ) )
    df_dday <- df_dday %>% mutate( iniwuebin  = cut( as.numeric(dday), breaks = iwuebins ) )

    ## add row: normalised VPD
    tmp <- df_dday %>% group_by( infvarbin ) %>% 
                       summarise( vpd  = median( vpd , na.rm=TRUE ) ) %>%
                       complete( infvarbin, fill = list( vpd  = NA ) ) %>% 
                       dplyr::select( vpd )
    tmp <- unlist( tmp )[1:(length(fvarbins)-1)]
    df_dday$dvpd = df_dday$vpd / tmp[1]

    ## add row: normalised soil moisture
    tmp <- df_dday %>% group_by( infvarbin ) %>% 
                       summarise( soilm_mean  = median( soilm_mean , na.rm=TRUE ) ) %>%
                       complete( infvarbin, fill = list( soilm_mean  = NA ) ) %>% 
                       dplyr::select( soilm_mean )
    tmp <- unlist( tmp )[1:(length(fvarbins)-1)]
    df_dday$soilm_norm = df_dday$soilm_mean / tmp[1]

    ## add row: normalised fAPAR (EVI)
    tmp <- df_dday %>% group_by( infaparbin ) %>% 
                       summarise( evi  = median( evi , na.rm=TRUE ) ) %>%
                       complete( infaparbin, fill = list( evi  = NA ) ) %>% 
                       dplyr::select( evi )
    tmp <- unlist( tmp )[1:(length(faparbins)-1)]
    df_dday$evi_norm = df_dday$evi / tmp[1]


    ## aggregate by 'dday'
    df_dday_aggbydday <- df_dday %>%  group_by( dday ) %>% 
                                      summarise(
                                                ## soil moisture
                                                soilm_med=median( soilm_mean, na.rm=TRUE ), soilm_upp=quantile( soilm_mean, 0.75, na.rm=TRUE ), soilm_low=quantile( soilm_mean, 0.25, na.rm=TRUE ),
                                                soilm_norm_med=median( soilm_norm, na.rm=TRUE ), soilm_norm_upp=quantile( soilm_norm, 0.75, na.rm=TRUE ), soilm_norm_low=quantile( soilm_norm, 0.25, na.rm=TRUE ),

                                                ## VPD
                                                vpd_med=median( vpd, na.rm=TRUE ), vpd_upp=quantile( vpd, 0.75, na.rm=TRUE ), vpd_low=quantile( vpd, 0.25, na.rm=TRUE ),

                                                ## relative VPD change
                                                dvpd_med=median( dvpd, na.rm=TRUE ), dvpd_upp=quantile( dvpd, 0.75, na.rm=TRUE ), dvpd_low=quantile( dvpd, 0.25, na.rm=TRUE ),

                                                ## fLUE
                                                fvar_med=median( fvar, na.rm=TRUE ), fvar_upp=quantile( fvar, 0.75, na.rm=TRUE ), fvar_low=quantile( fvar, 0.25, na.rm=TRUE ),

                                                ## EVI and FPAR
                                                evi_med=median( evi, na.rm=TRUE ), evi_upp=quantile( evi, 0.75, na.rm=TRUE ), evi_low=quantile( evi, 0.25, na.rm=TRUE ),
                                                fpar_med=median( fpar, na.rm=TRUE ), fpar_upp=quantile( fpar, 0.75, na.rm=TRUE ), fpar_low=quantile( fpar, 0.25, na.rm=TRUE )

                                                ) %>%
                                      mutate( mysitename=sitename )

    ## drop rows dday=NA
    df_dday_aggbydday <- df_dday_aggbydday[ which( !is.na(df_dday_aggbydday$dday)), ]
    save( df_dday_aggbydday, file=paste( "data/df_dday_aggbydday_", sitename, char_bysm, ".Rdata", sep="" ) )
    
    ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
    resave( df_dday, file=paste( "data/aligned_", sitename, char_bysm, ".Rdata", sep="" ) )

  } else {

    df_dday <- NA
    df_dday_aggbydday <- NA

  }

  out <- list( df_dday=df_dday, df_dday_aggbydday=df_dday_aggbydday )
  return( out )

}
