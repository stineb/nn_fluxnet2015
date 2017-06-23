reshape_align_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", use_fapar=FALSE, use_weights=FALSE, overwrite=TRUE, verbose=FALSE ){

  # ## debug-------------------
  # sitename = "IT-Cp2"
  # nam_target="lue_obs_evi"
  # use_fapar=FALSE
  # use_weights=FALSE
  # overwrite=TRUE
  # verbose=FALSE
  # #--------------------------

  require( dplyr )

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

	before <- 30
	after  <- 100

	load( paste( "data/missing_pri_", nam_target, char_fapar, ".Rdata", sep="") )

	## Bins for different variables
	fvarbins  <- seq( from=-20, to=40, by=20 )
	faparbins <- seq( from=-20, to=40, by=20 )
	iwuebins  <- seq( from=-30, to=60, by=30 )

	bincentres_fvar  <- fvarbins[1:(length(fvarbins)-1)]   + (fvarbins[2]-fvarbins[1])/2
	bincentres_fapar <- faparbins[1:(length(faparbins)-1)] + (faparbins[2]-faparbins[1])/2
	bincentres_iwue  <- iwuebins[1:(length(iwuebins)-1)]   + (iwuebins[2]-iwuebins[1])/2

  if ( !is.element( sitename, missing_pri ) ) { avl_pri <- TRUE } else { avl_pri <- FALSE }

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
  if (verbose) print("loading nn_fVAR file ...")
  infil <- paste( "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_fapar, ".Rdata", sep="" ) 
  load( infil ) ## gets list 'nn_fluxnet'
  df <- as.data.frame( nn_fluxnet[[ sitename ]]$nice ) %>% dplyr::select( year_dec, gpp_obs, var_nn_pot, var_nn_act, ppfd, fvar, soilm_mean, evi, fpar, wue_obs, is_drought_byvar, gpp_pmodel, gpp_obs_gfd, iwue, pri, cci, spri, scci )

  droughts <- nn_fluxnet[[ sitename ]]$droughts        

  df$bias_pmodel  <-  df$gpp_pmodel / df$gpp_obs
  df$bias_pmodel[ which(is.infinite(df$bias_pmodel)) ] <- NA

  df$ratio_obs_mod  <-  df$gpp_obs / df$gpp_pmodel
  df$ratio_obs_mod[ which(is.infinite(df$ratio_obs_mod)) ] <- NA

  names_alg <- c( names(df), "dday")

  if (nrow(droughts)>1){
    ##--------------------------------------------------------
    ## re-arrange data, aligning by beginning of drought
    ##--------------------------------------------------------
      filn <- paste( "data/aligned_", sitename, ".Rdata", sep="" )
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

      ## normalise PRI and CCI to mean by fvar-sized-bin number 2 (zerobin)
      if ( is.element( sitename, missing_pri ) ) { avl_pri <- FALSE } else { avl_pri <- TRUE }

      if ( avl_pri ){

        df_dday_meanbybin <- df_dday %>% group_by( infvarbin ) %>% summarise( cci  = mean( cci, na.rm=TRUE ), 
                                                                              pri  = mean( pri, na.rm=TRUE ), 
                                                                              scci = mean( scci, na.rm=TRUE ), 
                                                                              spri = mean( spri, na.rm=TRUE ) 
                                                                            )

        pri0 <- dplyr::filter( df_dday_meanbybin, infvarbin=='(-20,0]')$pri
        cci0 <- dplyr::filter( df_dday_meanbybin, infvarbin=='(-20,0]')$cci

        spri0 <- dplyr::filter( df_dday_meanbybin, infvarbin=='(-20,0]')$spri
        scci0 <- dplyr::filter( df_dday_meanbybin, infvarbin=='(-20,0]')$scci

        df_dday <- df_dday %>% mutate( 
                                      dpri  = pri  / pri0,
                                      dcci  = cci  / cci0,
                                      dspri = spri / spri0,
                                      dscci = scci / scci0
                                      )

      } else {

        df_dday$dpri  <- rep( NA, nrow(df_dday) )
        df_dday$dcci  <- rep( NA, nrow(df_dday) )
        df_dday$spri  <- rep( NA, nrow(df_dday) )
        df_dday$scci  <- rep( NA, nrow(df_dday) )
        df_dday$dspri <- rep( NA, nrow(df_dday) )
        df_dday$dscci <- rep( NA, nrow(df_dday) )

      }

      ## aggregate by 'dday'
      df_dday_aggbydday <- df_dday %>%  group_by( dday ) %>% 
                                        summarise(
                                                  ## soil moisture
                                                  soilm_med=median( soilm_mean, na.rm=TRUE ), soilm_upp=quantile( soilm_mean, 0.75, na.rm=TRUE ), soilm_low=quantile( soilm_mean, 0.25, na.rm=TRUE ),

                                                  ## fLUE
                                                  fvar_med=median( fvar, na.rm=TRUE ), fvar_upp=quantile( fvar, 0.75, na.rm=TRUE ), fvar_low=quantile( fvar, 0.25, na.rm=TRUE ),

                                                  ## EVI and FPAR
                                                  evi_med=median( evi, na.rm=TRUE ), evi_upp=quantile( evi, 0.75, na.rm=TRUE ), evi_low=quantile( evi, 0.25, na.rm=TRUE ),
                                                  fpar_med=median( fpar, na.rm=TRUE ), fpar_upp=quantile( fpar, 0.75, na.rm=TRUE ), fpar_low=quantile( fpar, 0.25, na.rm=TRUE ),

                                                  ## P-model bias
                                                  bias_pmodel_med=median( bias_pmodel, na.rm=TRUE ), bias_pmodel_upp=quantile( bias_pmodel, 0.75, na.rm=TRUE ), bias_pmodel_low=quantile( bias_pmodel, 0.25, na.rm=TRUE ),

                                                  ## PRI and CCI
                                                  dcci_med=median( dcci, na.rm=TRUE ),   dcci_upp=quantile( dcci, 0.75, na.rm=TRUE ),   dcci_low=quantile( dcci, 0.25, na.rm=TRUE ),
                                                  dpri_med=median( dpri, na.rm=TRUE ),   dpri_upp=quantile( dpri, 0.75, na.rm=TRUE ),   dpri_low=quantile( dpri, 0.25, na.rm=TRUE ),
                                                  dscci_med=median( dscci, na.rm=TRUE ), dscci_upp=quantile( dscci, 0.75, na.rm=TRUE ), dscci_low=quantile( dscci, 0.25, na.rm=TRUE ),
                                                  dspri_med=median( dspri, na.rm=TRUE ), dspri_upp=quantile( dspri, 0.75, na.rm=TRUE ), dspri_low=quantile( dspri, 0.25, na.rm=TRUE ),
                                                  scci_med=median( scci, na.rm=TRUE ),   scci_upp=quantile( scci, 0.75, na.rm=TRUE ),   scci_low=quantile( scci, 0.25, na.rm=TRUE ),
                                                  spri_med=median( spri, na.rm=TRUE ),   spri_upp=quantile( spri, 0.75, na.rm=TRUE ),   spri_low=quantile( spri, 0.25, na.rm=TRUE )

                                                  ) %>%
                                        mutate( mysitename=sitename )

      ## drop rows dday=NA
      df_dday_aggbydday <- df_dday_aggbydday[ which( !is.na(df_dday_aggbydday$dday)), ]
      save( df_dday_aggbydday, file=paste( "data/df_dday_aggbydday_", sitename, ".Rdata", sep="" ) )
      

      ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
      resave( df_dday, file=paste( "data/aligned_", sitename, ".Rdata", sep="" ) )


    ##--------------------------------------------------------
    ## re-arrange MODIS dataframe
    ##--------------------------------------------------------
      filn <- paste( "data/aligned_", sitename, ".Rdata", sep="" )
      if (!file.exists(filn)||overwrite){
        if (verbose) print("aligning MODIS ...")
        filn <- paste( "data/modis_", sitename, ".Rdata", sep="" )
        error <- try( load( filn ) ) # loads 'nice_to_modis', file prepared in 'plot_nn_fVAR_fluxnet2015.R'
        if (class(error)!="try-error"){
          
          if ( !is.element( "ratio_obs_mod_modis", names(nice_to_modis) ) ) {
            nice_to_modis$ratio_obs_mod_modis  <-  nice_to_modis$gpp_obs / nice_to_modis$gpp_modis
            nice_to_modis$ratio_obs_mod_modis[ which(is.infinite(nice_to_modis$ratio_obs_mod_modis)) ] <- NA
          }
          
          nice_to_modis <- nice_to_modis %>% dplyr::select( year_dec, gpp_obs, fvar, soilm_mean, is_drought_byvar, gpp_modis, bias_modis, ratio_obs_mod_modis )
          names_alg_modis <- c( names(nice_to_modis), "dday")

          droughts_modis <- get_consecutive( nice_to_modis$is_drought_byvar, leng_threshold=2, do_merge=FALSE )

          before_modis <- floor( max(droughts_modis$len) / 2 )
          after_modis  <- max(droughts_modis$len)
          data_alg_dry_modis  <- array( NA, dim=c( before_modis+after_modis+1, ncol(nice_to_modis)+1, nrow(droughts_modis) ) )
          
          for ( iinst in 1:nrow(droughts_modis) ){
            data_alg_dry_modis[,(ncol(nice_to_modis)+1),iinst] <- ((-before_modis:after_modis)+1)*8 ## calling this 'dday' = drought day
            for (idx in -before_modis:after_modis){
              if ( (droughts_modis$idx_start[iinst]+idx)>0 ){
                for (icol in 1:ncol(nice_to_modis)){
                  data_alg_dry_modis[ idx+before_modis+1, icol, iinst ] <- nice_to_modis[ droughts_modis$idx_start[iinst]+idx, icol ]
                }
              }
            }
            ## remove data after drought onset that is no longer classified as drought
            dropidxs <- which( data_alg_dry_modis[,which(names_alg_modis=="is_drought_byvar"),iinst]==1 & data_alg_dry_modis[,which(names_alg_modis=="dday"),iinst]<0 )
            data_alg_dry_modis[ dropidxs,,iinst ] <- NA
          }

          ##--------------------------------------------------------
          ## Bin aligned data and expand from 3D array to dataframe
          ##--------------------------------------------------------
          ## expand 'data_alg_dry_modis' to get a data frame that now has 'dday' in it
          df_dday_modis <- data.frame()
          for (iinst in seq(dim(data_alg_dry_modis)[3])){
            add <- as.data.frame( data_alg_dry_modis[,,iinst] )
            colnames(add) <- names_alg_modis
            add$mysitename <- rep( sitename, dim( data_alg_dry_modis[,,iinst] )[1])
            add$inst <- rep( iinst, dim( data_alg_dry_modis[,,iinst] )[1] )
            df_dday_modis <- rbind( df_dday_modis, add )
          }
          df_dday_modis <- df_dday_modis[ !is.na(df_dday_modis$year_dec), ]

          ## aggregate by 'dday'
          df_dday_aggbydday_modis <- df_dday_modis %>%  group_by( dday ) %>% 
            summarise(
                      ## MODIS bias
                      bias_modis_med=median( bias_modis, na.rm=TRUE ), bias_modis_upp=quantile( bias_modis, 0.75, na.rm=TRUE ), bias_modis_low=quantile( bias_modis, 0.25, na.rm=TRUE )
                      ) %>%
            mutate( mysitename=sitename )

          ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
          save( data_alg_dry_modis, df_dday_modis, df_dday_aggbydday_modis, names_alg_modis, before_modis, after_modis, file=paste( "data/aligned_modis_", sitename, ".Rdata", sep="" ) )

        } else {

          df_dday_modis <- NA

        }

      } else {

        load(filn)

      }


    ##--------------------------------------------------------
    ## re-arrange MTE dataframe
    ##--------------------------------------------------------
      filn <- paste( "data/aligned_", sitename, ".Rdata", sep="" )
      if (!file.exists(filn)||overwrite){
        if (verbose) print("aligning MTE ...")
        filn <- paste( "data/mte_", sitename, ".Rdata", sep="" )
        avl_mte <- FALSE

        if (file.exists(filn)){

          avl_mte <- TRUE
          load( filn )

          nice_to_mte <- nice_to_mte %>% dplyr::select( year_dec, gpp_obs, fvar, soilm_mean, is_drought_byvar, gpp_mte, bias_mte, ratio_obs_mod_mte )
          names_alg_mte <- c( names(nice_to_mte), "dday")

          droughts_mte <- get_consecutive( nice_to_mte$is_drought_byvar, leng_threshold=2, do_merge=FALSE )

          if (nrow(droughts_mte)>0){
            before_mte <- floor( max(droughts_mte$len) / 2 )
            after_mte  <- max(droughts_mte$len)
            data_alg_dry_mte  <- array( NA, dim=c( before_mte+after_mte+1, ncol(nice_to_mte)+1, nrow(droughts_mte) ) )
            
            for ( iinst in 1:nrow(droughts_mte) ){
              data_alg_dry_mte[,(ncol(nice_to_mte)+1),iinst] <- ((-before_mte:after_mte)+1)*8 ## calling this 'dday' = drought day
              for (idx in -before_mte:after_mte){
                if ( (droughts_mte$idx_start[iinst]+idx)>0 ){
                  for (icol in 1:ncol(nice_to_mte)){
                    data_alg_dry_mte[ idx+before_mte+1, icol, iinst ] <- nice_to_mte[ droughts_mte$idx_start[iinst]+idx, icol ]
                  }
                }
              }
              ## remove data after drought onset that is no longer classified as drought
              dropidxs <- which( data_alg_dry_mte[,which(names_alg_mte=="is_drought_byvar"),iinst]==1 & data_alg_dry_mte[,which(names_alg_mte=="dday"),iinst]<0 )
              data_alg_dry_mte[ dropidxs,,iinst ] <- NA
            }        
          } else {
            avl_mte <- FALSE
          }

          ##--------------------------------------------------------
          ## Bin aligned data and expand from 3D array to dataframe
          ##--------------------------------------------------------
          ## expand 'data_alg_dry_mte' to get a data frame that now has 'dday' in it
          df_dday_mte <- data.frame()
          for (iinst in seq(dim(data_alg_dry_mte)[3])){
            add <- as.data.frame( data_alg_dry_mte[,,iinst] )
            colnames(add) <- names_alg_mte
            add$mysitename <- rep( sitename, dim( data_alg_dry_mte[,,iinst] )[1])
            add$inst <- rep( iinst, dim( data_alg_dry_mte[,,iinst] )[1] )
            df_dday_mte <- rbind( df_dday_mte, add )
          }
          df_dday_mte <- df_dday_mte[ !is.na(df_dday_mte$year_dec), ]

          ## aggregate by 'dday'
          df_dday_aggbydday_mte <- df_dday_mte %>%  group_by( dday ) %>% 
            summarise(
                      ## mte bias
                      bias_mte_med=median( bias_mte, na.rm=TRUE ), bias_mte_upp=quantile( bias_mte, 0.75, na.rm=TRUE ), bias_mte_low=quantile( bias_mte, 0.25, na.rm=TRUE )
                      ) %>%
            mutate( mysitename=sitename )

          ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
          save( data_alg_dry_mte, df_dday_mte, df_dday_aggbydday_mte, names_alg_mte, before_mte, after_mte, file=paste( "data/aligned_mte_", sitename, ".Rdata", sep="" ) )

        } else {

        	df_dday_mte <- NA

        }


      } else {

        load( filn )

      }

  } else {

    df_dday <- NA
    df_dday_aggbydday <- NA
    df_dday_modis <- NA
    df_dday_mte <- NA

  }

  out <- list( df_dday=df_dday, df_dday_aggbydday=df_dday_aggbydday, df_dday_modis=df_dday_modis, df_dday_mte=df_dday_mte )
  return( out )

}
