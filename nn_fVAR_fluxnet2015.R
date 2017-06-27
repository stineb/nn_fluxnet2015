nn_fVAR_fluxnet <- function( sitename, nam_target="lue_obs_evi", use_weights=ifelse( nam_target=="lue_obs_evi", TRUE, FALSE ), use_fapar=FALSE, nrep=5, dotrain=FALSE, package="nnet", makepdf=TRUE, testprofile=FALSE, do_evianom=FALSE ){

  # ## XXX debug----------------
  # sitename   = "FR-Pue"
  # nam_target = "lue_obs_evi"
  # use_weights= FALSE    
  # use_fapar  = FALSE
  # package    = "nnet"
  # nrep       = 5
  # dotrain    = FALSE
  # makepdf    = FALSE
  # testprofile= FALSE
  # ##--------------------------

  require( dplyr )
  require( abind )
  
  source( paste( workingdir, "/analyse_modobs.R", sep="" ) )
  source( paste( workingdir, "/add_alpha.R", sep="" ) )
  source( paste( workingdir, "/get_consecutive.R", sep="" ) )
  source( paste( workingdir, "/predict_nn.R", sep="" ) )
  source( paste( workingdir, "/remove_outliers.R", sep="" ) )
  source( paste( workingdir, "/niceify.R", sep="" ) )
  source( paste( workingdir, "/get_evianomalies_fluxnet.R", sep="" ) )
  source( paste( workingdir, "/cutna_headtail.R", sep="" ) )
  source( paste( workingdir, "/prune_droughts.R", sep="" ) )
  source( paste( workingdir, "/get_iwue.R", sep="" ) )
  source( paste( workingdir, "/cleandata_nn.R", sep="" ) )

  ## IMPORTANT: USE SOILMOISTURE FROM S13 FOR NN-TRAINING
  load( paste( workingdir, "/data/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata", sep="" ) ) # "new data" with s13

  nn_fluxnet <- list()

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

  if (use_weights){
    char_wgt <- "_wgt"
  } else {
    char_wgt <- ""
  }

  ## Get NN-soilm-profile data
  if (testprofile){
    filn <- paste( workingdir, "/data/profile/profile_", nam_target, char_wgt, char_fapar, "_nn_", sitename, ".Rdata", sep="" )
  } else {
    filn <- paste( myhome, "/data/nn_fluxnet/profile/profile_", nam_target, char_wgt, char_fapar, "_nn_", sitename, ".Rdata", sep="" )
  }

  deleteagain <- FALSE
  print( paste( "loading profile file:", filn ) )
  if (!file.exists(filn)){
    print( "downloading profile file ..." )
    system( paste("rsync -avz bstocker@login.cx1.hpc.ic.ac.uk:/home/bstocker/data/nn_fluxnet/profile/profile_", nam_target, char_wgt, char_fapar, "_nn_", sitename, ".Rdata"," /alphadata01/bstocker/data/nn_fluxnet/profile/", sep="" ) )
    deleteagain <- TRUE
    print("done.")
  }

  out <- try( load( filn ) )
  if ( class(out)!="try-error" ){

    ##------------------------------------------------
    ## Extract data
    ##------------------------------------------------
    if ( !is.na( profile_nn[[ sitename ]][1] ) ){
      best_soilm     <- profile_nn[[ sitename ]]$best_smdata[1]
      best_soilm_trh <- profile_nn[[ sitename ]][[ best_soilm ]]$nnet$best_soilm_trh[1]
      idxs_good      <- profile_nn[[ sitename ]][[ best_soilm ]]$nnet[[ paste("smtrh_", best_soilm_trh, sep="" ) ]]$idxs_good
      var_nn_bad     <- apply( profile_nn[[ sitename ]][[ best_soilm ]]$nnet[[ paste("smtrh_", best_soilm_trh, sep="" ) ]]$var_nn_bad, 1, FUN=mean )
      var_nn_all     <- apply( profile_nn[[ sitename ]][[ best_soilm ]]$nnet$var_nn_all, 1, FUN=mean )
    }

    ##------------------------------------------------
    ## Get neural network for this site, each target individually and looping over repetitions to get statistics
    ##------------------------------------------------
    print( "===============================================")
    print( paste( "train NN at", nam_target, "for", sitename, "..." ) )

    ##------------------------------------------------
    ## Get model output from simulation without temperature or soil moisture limitation (s10)
    ##------------------------------------------------
    data <- fluxnet[[ sitename ]]$ddf$s13

    if (!is.null(data)){

      data <- data %>% dplyr::select( year_dec, year, doy, moy, dom, soilm_splash150=wcont, gpp_pmodel=gpp, aet_pmodel=aet, pet_pmodel=pet )

      ##------------------------------------------------
      ## Get alternative soil moisture data
      ##------------------------------------------------
      data$soilm_splash220 <- fluxnet[[ sitename ]]$ddf$s11$wcont
      data$soilm_swbm      <- fluxnet[[ sitename ]]$ddf$s12$wcont
      data$soilm_etobs     <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et
      data$soilm_etobs_ob  <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et_orthbucket

      varnams_swc <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob" )

      ##------------------------------------------------
      ## Get observational soil moisture data (variable availability!)
      ##------------------------------------------------
      relevant <- names(fluxnet[[ sitename ]]$ddf$swc_obs)[(is.element( substr(names(fluxnet[[ sitename ]]$ddf$swc_obs), start=1, stop=3), "SWC" ))]
      varnams_swc_obs <- c()
      if (length(relevant)>0){
        for (iobs in seq(length(relevant))){
          varnam <- paste( "soilm_obs_", iobs, sep="" )
          data[[ varnam ]] <- fluxnet[[ sitename ]]$ddf$swc_obs[[ relevant[iobs] ]]
          varnams_swc_obs <- c( varnams_swc_obs, varnam )
        }
        varnams_swc <- c( varnams_swc, "soilm_obs" )
      }

      ##------------------------------------------------
      ## Get observational data and add to 'data'
      ##------------------------------------------------
      obs <- dplyr::select( fluxnet[[ sitename ]]$ddf$obs, year, moy, dom, gpp_obs2015_GPP_NT_VUT_REF, gpp_obs2015_GPP_NT_VUT_REF_gfd, le_f_mds )

      data <- data  %>% left_join( obs, by=c( "year", "moy", "dom" ) ) %>%
                        rename( gpp_obs=gpp_obs2015_GPP_NT_VUT_REF, gpp_obs_gfd=gpp_obs2015_GPP_NT_VUT_REF_gfd, et_obs=le_f_mds ) %>%
                        mutate( wue_obs=gpp_obs/(et_obs*1e-6) )

      ##------------------------------------------------
      ## Get input data
      ##------------------------------------------------
      inp <- fluxnet[[ sitename ]]$ddf$inp %>% dplyr::select( year, moy, dom, temp, ppfd, vpd, prec, evi, fpar )

      data <- data %>% left_join( inp, by=c( "year", "moy", "dom" ) )

      ##------------------------------------------------
      ## Get IWUE* and remove outliers 
      ##------------------------------------------------
      data <- data %>% mutate( iwue = remove_outliers( get_iwue( et_obs, gpp_obs, vpd, prec ), coef=3.0 ) )

      ##------------------------------------------------
      ## Remove outliers in WUE
      ##------------------------------------------------
      data <- data %>% mutate( wue_obs = remove_outliers( wue_obs, coef=3.0 ) )

      ##------------------------------------------------
      ## get LUE and remove outliers
      ##------------------------------------------------
      data <- data %>% mutate( lue_obs_evi  = remove_outliers( gpp_obs / ( ppfd * evi  ), coef=3.0 ) )
      data <- data %>% mutate( lue_obs_fpar = remove_outliers( gpp_obs / ( ppfd * fpar ), coef=3.0 ) )      

      ##------------------------------------------------
      ## normalise soil moisture
      ##------------------------------------------------
      data$soilm_splash150 <- data$soilm_splash150 / max( data$soilm_splash150, na.rm=TRUE )
      data$soilm_splash220 <- data$soilm_splash220 / max( data$soilm_splash220, na.rm=TRUE )
      data$soilm_swbm      <- data$soilm_swbm      / max( data$soilm_swbm     , na.rm=TRUE )
      data$soilm_etobs     <- data$soilm_etobs     / max( data$soilm_etobs    , na.rm=TRUE )
      data$soilm_etobs_ob  <- data$soilm_etobs_ob  / max( data$soilm_etobs_ob , na.rm=TRUE )
      for (ivar in varnams_swc_obs){
        data[[ ivar ]]  <- data[[ ivar ]]  / max( data[[ ivar ]] , na.rm=TRUE )
      }

      data_full_save <- data

      # ##------------------------------------------------
      # ## Use best-performing soil moisture data for NN-fGPP algorithm
      # ##------------------------------------------------
      # isoilm_data <- profile_nn[[ sitename ]]$best_smdata[1]
      # # isoilm_data <- "soilm_splash220"

      ##------------------------------------------------
      ## Use all soil moisture data to get the full uncertainty range
      ##------------------------------------------------
      var_nn_pot_by_smdata <- c()
      var_nn_act_by_smdata <- c()
      fvar_by_smdata       <- c()

      do_smdata <- profile_nn[[ sitename ]]$varnams_swc
      # do_smdata <- c("soilm_splash220", "soilm_swbm")

      # # xxx try:
      # do_smdata <- do_smdata[-which(do_smdata=="soilm_obs")]

      lab_smdata <- do_smdata

      for ( isoilm_data in do_smdata ){
        
        print(paste("Using soil moisture data source:", isoilm_data))
        
        ##------------------------------------------------
        ## Use only days where temperature is above 5 degrees
        ##------------------------------------------------
        df_nona <- data
        # df_nona <- dplyr::filter( df_nona, fapar > quantile( data[[ fapar_data ]], probs=0.25 ) )
        df_nona <- dplyr::filter( df_nona, temp > 5.0 )

        ##------------------------------------------------
        ## do additional data cleaning, removing NA in target variable, necessary for NN training
        ##------------------------------------------------
        df_nona <- cleandata_nn( df_nona, nam_target ) # doesn't produce big additional gap

        if ( isoilm_data=="soilm_obs" ){
      
          ## use obs soilm data only if of sufficient length
          lengths <- apply( subset( df_nona, select=varnams_swc_obs ), 2, function(x) sum(!is.na(x)) )

          ## drop layer swc obs data if length of data is less than 75% of legth of maximum
          idx <- 0
          drop_idx <- c()
          for (ivar in varnams_swc_obs){
            idx <- idx + 1
            if (lengths[ivar]<0.75*max(lengths)){
              df_nona[[ ivar ]] <- NULL
              drop_idx <- c(drop_idx, idx)
            }
          }
          if ( length(drop_idx)>0 ) { varnams_swc_obs <- varnams_swc_obs[-drop_idx] }
          
          ## remove NAs in observed soil moisture data
          for (ivar in varnams_swc_obs){
            df_nona <- df_nona[ which(!is.na(df_nona[[ivar]])), ]  ## XXX THIS CAUSES THE PROBLEM OF REMOVING THE ENTIRE DATA WHEN SOILM_OBS DATA IS MISSING XXX
          }

        }

        ## get weights for NN training: absorbed light
        if (use_weights) {
          weights <- df_nona$ppfd * df_nona$evi
        } else {
          weights <- rep( 1.0, nrow(df_nona) )
        }

        ##------------------------------------------------
        ## Use best-performing soil moisture threshold for NN-fGPP algorithm
        ##------------------------------------------------
        ## get soilm_thrsh_avl where best fit is achieved for good-days model (lowest variance of fVAR during good days)
        df_stat_good <- profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$df_stat
        diff_good_bad <- profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$diff_good_bad
        if (is.null(diff_good_bad)){
          isoilm_trh <- NA
        } else {
          n_best <- 5
          list_best_soilm_trh <- diff_good_bad$soilm_threshold[ order(-diff_good_bad$diff)][1:min(n_best,nrow(diff_good_bad))] ## N best thresholds based on magnitude of difference
          sub <- df_stat_good[ is.element( df_stat_good$soilm_thrsh_avl, list_best_soilm_trh ), ]
          best_soilm_trh <- sub$soilm_thrsh_avl[ order(sub$var)]
          isoilm_trh <- best_soilm_trh[1]
        }

        print( paste( "Using soil moisture threshold:", isoilm_trh ) )

        if ( dotrain ) {

          print("Doing NN training again ...")

          ##------------------------------------------------
          ## Get best-performing number of hidden layers for this site and soil moisture data source
          ##------------------------------------------------
          hidden_good <-  profile_nn[[ sitename ]][[ isoilm_data ]]$nnet[[ paste( "smtrh_", as.character( isoilm_trh ), sep="" ) ]]$hidden_best_good
          hidden_all  <-  profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$hidden_best_all


          ##------------------------------------------------
          ## Determine "good days", i.e. where soil moisture is abover threshold
          ##------------------------------------------------
          if ( isoilm_data=="soilm_obs" ){
            ## for observational data, do subset w.r.t. soil layer with highest value
            tmp <- subset( df_nona, select=varnams_swc_obs )
            tmp[ is.na(tmp) ] <- 0.0
            vec <- apply( tmp, 1, FUN = max, na.rm=TRUE )
          } else {
            vec <- df_nona[[ isoilm_data ]]
          }

          idxs_good <- which( vec > isoilm_trh )

          if ( length(idxs_good)>30 && (nrow(df_nona) - length(idxs_good))>30 ){

            ## Initialise
            nn_act_vals <- array( NA, dim=c(nrow(df_nona),nrep) )
            nn_pot_vals <- array( NA, dim=c(nrow(df_nona),nrep) )
            nn_fXX_vals <- array( NA, dim=c(nrow(df_nona),nrep) )

            for (irep in 1:nrep){
              print(paste("NN repetition", irep))
              ##------------------------------------------------
              ## Good days GPP (LUE)
              ##------------------------------------------------
              ## Get good-days-NN based on subset of data where soil moisture is above threshold
              # print("training good-days model ...")

              predictors <- c( "ppfd", "temp", "vpd" )
              if ( use_fapar ){ predictors <- c( predictors, fapar_data ) }

              out_nn_good <- predict_nn( 
                data       = df_nona[ idxs_good, ],
                weights    = weights[ idxs_good ],
                predictors = predictors,
                nam_target = nam_target,
                do_predict = TRUE,
                package    = package,
                lifesign   = "full",
                seed       = irep,
                hidden     = hidden_good
                )

              ## Evaluate predictions of good days model
              stats_nn_good <- analyse_modobs( 
                                              apply( var_nn_good, 1, FUN=mean ) * weights[ idxs_good ], 
                                              df_nona$gpp_obs[idxs_good], 
                                              plot.title=paste( "Good, trh=", isoilm_trh, "data=", isoilm_data ), 
                                              do.plot=FALSE
                                              )

              ##------------------------------------------------
              ## Bad days GPP (LUE), predicted using the model trained on good days data
              ##------------------------------------------------
              # print("predicting bad days with good-days model ...")
              out_nn_pot <- predict_nn( 
                data       = df_nona, 
                weights    = weights,
                predictors = predictors, 
                nam_target = nam_target, 
                nn         = out_nn_good$nn, 
                do_predict = TRUE, 
                package    = package
                )

              ## Evaluate predictions of good days model
              stats_nn_bad <- analyse_modobs( 
                                              (apply( var_nn_bad, 1,  FUN=mean ) * weights)[-idxs_good], 
                                              df_nona$gpp_obs[-idxs_good], 
                                              plot.title=paste( "Bad, trh=", isoilm_trh, "data=", isoilm_data ), 
                                              do.plot=FALSE
                                              )

              print( paste( "R2 of NN-good with soil moisture threshold ", isoilm_trh, "is", stats_nn_good$rsq ) )


              ##------------------------------------------------
              ## Get full GPP-model, including soil moisture as a predictor ("act")
              ##------------------------------------------------
              # print("training full GPP model ...")
              if ( isoilm_data=="soilm_obs" ){
                ## include the full list of soil moisture variables
                soilm_predictors <- varnams_swc_obs
              } else {
                soilm_predictors <- isoilm_data
              }

              predictors <- c( "ppfd", "temp", "vpd", soilm_predictors )
              if ( use_fapar ){ predictors <- c( predictors, fapar_data ) }

              out_nn_act <- predict_nn( 
                data       = df_nona, 
                weights    = weights,
                predictors = predictors, 
                nam_target = nam_target, 
                do_predict = TRUE, 
                package    = package,
                seed       = irep,
                hidden     = hidden_all
                )

              ## get statistics of mod vs. obs of all-days full model
              stats_nn_all <- analyse_modobs( 
                                              apply( var_nn_all, 1, FUN=mean ) * weights, 
                                              df_nona$gpp_obs, 
                                              plot.title=paste( "All, data=", isoilm_data ), 
                                              do.plot=FALSE
                                              )

              ## keep only values
              nn_pot_vals[,irep] <- as.vector( out_nn_pot$vals )
              nn_act_vals[,irep] <- as.vector( out_nn_act$vals )
              nn_fXX_vals[,irep] <- as.vector( out_nn_act$vals ) / as.vector( out_nn_pot$vals )

            }

          } else {

            print(paste("too few data points with soil moisture threshold", isoilm_trh))

            # print(paste("too few data points with soil moisture threshold", isoilm))
            nn_pot_vals[,] <- NA
            nn_act_vals[,] <- NA
            nn_fXX_vals[,] <- NA

          }

          ##------------------------------------------------
          ## Aggregate over repetitions
          ##------------------------------------------------
          if (nrep>1){

            df_nona$var_nn_pot    <- apply( nn_pot_vals[,], 1, FUN = mean )
            df_nona$var_nn_act    <- apply( nn_act_vals[,], 1, FUN = mean )

            df_nona$var_nn_pot_sd <- apply( nn_pot_vals[,], 1, FUN = sd )
            df_nona$var_nn_act_sd <- apply( nn_act_vals[,], 1, FUN = sd )

            df_nona$fvar          <- apply( nn_fXX_vals[,], 1, FUN = mean )
            df_nona$fvar_sd       <- apply( nn_fXX_vals[,], 1, FUN = sd )
         
          } else {

            df_nona$var_nn_pot    <- nn_pot_vals[,1]
            df_nona$var_nn_act    <- nn_act_vals[,1]

            df_nona$var_nn_pot_sd <- nn_pot_vals[,1]
            df_nona$var_nn_act_sd <- nn_act_vals[,1]

            df_nona$fvar          <- nn_fXX_vals[,1]
            df_nona$fvar_sd       <- nn_fXX_vals[,1]
          }

        } else if (!is.na(isoilm_trh)) {

          print("reading from profile data...")
          
          ## take mean over N repetitions
          var_nn_pot <- apply( profile_nn[[ sitename ]][[ isoilm_data ]]$nnet[[ paste( "smtrh_", as.character( isoilm_trh ), sep="" ) ]]$var_nn_bad, 1, FUN=mean )
          var_nn_act <- apply( profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$var_nn_all, 1, FUN=mean )
          fvar       <- apply( profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$var_nn_all / profile_nn[[ sitename ]][[ isoilm_data ]]$nnet[[ paste( "smtrh_", as.character( isoilm_trh ), sep="" ) ]]$var_nn_bad, 1, FUN=mean )
          idxs_good  <- profile_nn[[ sitename ]][[ isoilm_data ]]$nnet[[ paste( "smtrh_", as.character( isoilm_trh ), sep="" ) ]]$idxs_good

          tmp <- data.frame( 
            var_nn_pot = var_nn_pot,
            var_nn_act = var_nn_act,
            fvar       = fvar,
            year_dec   = profile_nn[[ sitename ]][[ isoilm_data ]]$nnet$year_dec
          )

          if (nam_target=="gpp_obs"){
            tmp$gpp_obs <- profile_nn[[ sitename ]][[isoilm_data]]$nnet$gpp_obs
          } else if (nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar"){
            tmp$lue_obs <- profile_nn[[ sitename ]][[isoilm_data]]$nnet$lue_obs
          }

          ## expand to all days, filling with NA
          tmp <- niceify( tmp, df_nona )
          df_nona$var_nn_pot    <- tmp$var_nn_pot
          df_nona$var_nn_act    <- tmp$var_nn_act
          df_nona$fvar          <- tmp$fvar

          print("... done.")

        }

        if (!is.na(isoilm_trh)) {
          ##------------------------------------------------
          ## Remove outliers
          ##------------------------------------------------
          df_nona$fvar <- remove_outliers_fXX( df_nona$fvar, coef=3.0 )    

          ##------------------------------------------------
          ## Niceify and select only relevant
          ##------------------------------------------------
          tmp <- dplyr::select( niceify( df_nona, data ), year_dec, var_nn_pot, var_nn_act, fvar )

          ##------------------------------------------------
          ## record for later aggregation across soilm. datasets
          ##------------------------------------------------
          var_nn_pot_by_smdata <- cbind( var_nn_pot_by_smdata, tmp$var_nn_pot )
          var_nn_act_by_smdata <- cbind( var_nn_act_by_smdata, tmp$var_nn_act )
          fvar_by_smdata       <- cbind( fvar_by_smdata,       tmp$fvar       )

        } else {

          lab_smdata <- lab_smdata[ -which(lab_smdata==isoilm_data) ]

        }

      } ## END OF DO_SMDATA LOOP

      ##------------------------------------------------
      ## Create nice dataframe that has same dimensions as original data
      ##------------------------------------------------
      nice <- data

      ##------------------------------------------------
      ## Remove NAs at head and tail based on 'nam_target'
      ##------------------------------------------------
      drop_idx <- cutna_headtail( nice[[ nam_target ]] )
      if (length(drop_idx)>0){
        var_nn_pot_by_smdata   <- var_nn_pot_by_smdata[ -drop_idx, ] 
        var_nn_act_by_smdata   <- var_nn_act_by_smdata[ -drop_idx, ] 
        fvar_by_smdata         <- fvar_by_smdata[ -drop_idx, ]    
        nice                   <- nice[ -drop_idx, ]
      }

      ##------------------------------------------------
      ## Aggregate over nn_pot, nn_act, and fVAR derived using different soil moisture datasets
      ##------------------------------------------------
      print("aggregate across soil moisture datasets ...")

      if (length(lab_smdata)>1){

        colnames(var_nn_pot_by_smdata) <- lab_smdata
        colnames(var_nn_act_by_smdata) <- lab_smdata
        colnames(fvar_by_smdata)       <- lab_smdata

        nice$var_nn_pot <- apply( var_nn_pot_by_smdata, 1, FUN=mean, na.rm=TRUE ); nice$var_nn_pot[ is.nan(nice$var_nn_pot) ] <- NA
        nice$var_nn_act <- apply( var_nn_act_by_smdata, 1, FUN=mean, na.rm=TRUE ); nice$var_nn_act[ is.nan(nice$var_nn_act) ] <- NA
        nice$fvar       <- apply( fvar_by_smdata, 1, FUN=mean, na.rm=TRUE )      ; nice$fvar[ is.nan(nice$fvar) ]             <- NA

        ## take that of best-performing soilm. dataset as standard
        if (is.element( "soilm_obs", lab_smdata )){
          nice$fvar_obs <- fvar_by_smdata[ , "soilm_obs" ]  
        } else {
          nice$fvar_obs <- rep(NA, nrow(nice))
        }

        var_nn_pot_median <- apply( var_nn_pot_by_smdata, 1, FUN=median, na.rm=TRUE )
        var_nn_act_median <- apply( var_nn_act_by_smdata, 1, FUN=median, na.rm=TRUE )
        fvar_med          <- apply( fvar_by_smdata, 1, FUN=median, na.rm=TRUE )

        var_nn_pot_max <- apply( var_nn_pot_by_smdata, 1, FUN=max, na.rm=TRUE ); var_nn_pot_max[ is.infinite(var_nn_pot_max) ] <- NA
        var_nn_act_max <- apply( var_nn_act_by_smdata, 1, FUN=max, na.rm=TRUE ); var_nn_act_max[ is.infinite(var_nn_act_max) ] <- NA
        fvar_max       <- apply( fvar_by_smdata, 1, FUN=max, na.rm=TRUE );       fvar_max[ is.infinite(fvar_max) ] <- NA

        var_nn_pot_min <- apply( var_nn_pot_by_smdata, 1, FUN=min, na.rm=TRUE ); var_nn_pot_min[ is.infinite(var_nn_pot_min) ] <- NA
        var_nn_act_min <- apply( var_nn_act_by_smdata, 1, FUN=min, na.rm=TRUE ); var_nn_act_min[ is.infinite(var_nn_act_min) ] <- NA
        fvar_min       <- apply( fvar_by_smdata, 1, FUN=min, na.rm=TRUE );       fvar_min[ is.infinite(fvar_min) ] <- NA

      } else {

        nice$var_nn_pot <- var_nn_pot_by_smdata
        nice$var_nn_act <- var_nn_act_by_smdata
        nice$fvar       <- fvar_by_smdata

        ## take that of best-performing soilm. dataset as standard
        if (is.element( "soilm_obs", lab_smdata )){
          nice$fvar_obs <- fvar_by_smdata[ , "soilm_obs" ]  
        } else {
          nice$fvar_obs <- rep(NA, nrow(nice))
        }

        var_nn_pot_median <- var_nn_pot_by_smdata
        var_nn_act_median <- var_nn_act_by_smdata
        fvar_med          <- fvar_by_smdata

        var_nn_pot_max <- var_nn_pot_by_smdata
        var_nn_act_max <- var_nn_act_by_smdata
        fvar_max       <- fvar_by_smdata

        var_nn_pot_min <- var_nn_pot_by_smdata
        var_nn_act_min <- var_nn_act_by_smdata
        fvar_min       <- fvar_by_smdata

      }


      print("done ...")


      ##------------------------------------------------
      ## Get cutoff and plot histogram on non-filled fGPP
      ##------------------------------------------------
      ## define drought by finding the lower 5% quantile of the hypothetical symetrical distribution
      print("get cutoff ...")
      positive <- nice$fvar[ which(nice$fvar>=1.0) ]
      even <- c( positive, 1.0-(positive-1.0) )
      cutoff <- quantile( even, 0.05, na.rm=TRUE )

      ## consider only instances where fvar falls below 0.97 as drought
      cutoff <- min( 0.97, cutoff )

      print("done ...")

      if (length(positive)>0){

        if (makepdf) pdf( paste( workingdir, "/fig_nn_fluxnet2015/hist/hist_fgpp_", package, "_", sitename, "_", as.character(isoilm_trh), ".pdf", sep="" ), width=8, height=6 )
          par( xaxs="i", las=1 )
          hist( even,      
                breaks=c(min(even),seq(0, 1.5, 0.01),max(even)),           
                xlim=c(0,1.5), 
                col=add_alpha("blue", 0.4), 
                freq=FALSE, 
                xlab="fLUE", 
                main=paste(sitename) 
                )
          hist( nice$fvar, 
                breaks=c(min(nice$fvar, na.rm=TRUE),seq(0, 1.5, 0.01),max(nice$fvar, na.rm=TRUE)), 
                xlim=c(0,1.5), 
                col=add_alpha("red", 0.4),  
                freq=FALSE, 
                add=TRUE 
                )
          abline( v=cutoff, col="red" )
        if (makepdf) dev.off()


        ##------------------------------------------------
        ## Fill and smooth fGPP (fLUE) in case gaps are smaller than 'fill_threshold' days
        ##------------------------------------------------
        print("fill and smooth ...")
        fill_threshold <- 20

        ## identify and store data gaps bigger than 'fill_threshold'
        na_instances <- get_consecutive( is.na(nice$fvar), leng_threshold=fill_threshold, do_merge=FALSE )

        ## "fill" data: interpolate to missing values
        nice$fvar_filled <- approx( nice$year_dec, nice$fvar, xout=nice$year_dec )$y

        ## take running median to smooth data
        idxs <- which( !is.na(nice$fvar_filled) )
        tmp <- data.frame( year_dec=nice$year_dec[idxs], fvar_smooth=runmed( nice$fvar_filled[idxs], 5 ) )
        nice$fvar_smooth <- niceify( tmp, nice )$fvar_smooth

        ## re-create gaps defined by 'na_instances'
        if (nrow(na_instances)>0){
          for (idx in 1:nrow(na_instances)){
            nice$fvar_filled[ na_instances$idx_start[idx]:(na_instances$idx_start[idx]+na_instances$len[idx]-1) ] <- NA
            nice$fvar_smooth[ na_instances$idx_start[idx]:(na_instances$idx_start[idx]+na_instances$len[idx]-1) ] <- NA
          }
        }

        ## "fill" NAs and create 'minmax' data frame for minimum and maximum fvar by soil moisture
        fvar_min_filled <- approx( nice$year_dec, fvar_min, xout=nice$year_dec )$y
        fvar_max_filled <- approx( nice$year_dec, fvar_max, xout=nice$year_dec )$y
        fvar_med_filled <- approx( nice$year_dec, fvar_med, xout=nice$year_dec )$y
        minmax <- data.frame( 
          year_dec        = nice$year_dec  [ which(!is.na(fvar_min_filled ) ) ],
          fvar_min_filled = fvar_min_filled[ which(!is.na(fvar_min_filled ) ) ],
          fvar_max_filled = fvar_max_filled[ which(!is.na(fvar_max_filled ) ) ],
          fvar_med_filled = fvar_med_filled[ which(!is.na(fvar_med_filled ) ) ]
          )

        print("done ...")

        ##------------------------------------------------
        ## Get observational soil moisture data (availability different for each site)
        ##------------------------------------------------
        ## all
        relevant <- names( nice )[ is.element( substr(names(nice), start=1, stop=6), "soilm_" ) ]
        varnams_swc <- c()
        if (length(relevant)>0){
          for (iobs in seq(length(relevant))){
            varnam <- relevant[iobs]
            varnams_swc <- c( varnams_swc, varnam )
          }
        }

        ## only obs
        relevant <- names( nice )[ is.element( substr(names(nice), start=1, stop=10), "soilm_obs_" ) ]
        varnams_swc_obs <- c()
        if (length(relevant)>0){
          for (iobs in seq(length(relevant))){
            varnam <- relevant[iobs]
            varnams_swc_obs <- c( varnams_swc_obs, varnam )
          }
        }

        ## Average over soil moisture datasets; all and obs-only
        nice$soilm_mean <- apply( dplyr::select( nice, one_of(varnams_swc)), 1, FUN=mean, na.rm=TRUE )
        nice$soilm_mean[ is.nan( nice$soilm_mean ) ] <- NA
        if (length(varnams_swc_obs)>0){
          nice$soilm_obs_mean <- apply( dplyr::select( nice, one_of(varnams_swc_obs)), 1, FUN=mean, na.rm=TRUE )
          nice$soilm_obs_mean[ is.nan( nice$soilm_obs_mean ) ] <- NA
        }


        ##------------------------------------------------
        ## Get fLUE droughts
        ##------------------------------------------------
        print("get drought events ...")
        nice$fvar_smooth_filled <- approx( nice$year_dec, nice$fvar_smooth, xout=nice$year_dec )$y
        nice$is_drought_byvar <- ( nice$fvar_smooth_filled < cutoff )
        droughts <- get_consecutive( 
                                    nice$is_drought_byvar, 
                                    leng_threshold = 10, 
                                    do_merge       = FALSE
                                    )

        ##------------------------------------------------
        ## Prune identified drought events and update is_drought_byvar
        ##------------------------------------------------
        soilm     <- unlist( ifelse( is.element("soilm_obs_mean",names(nice)), dplyr::select( nice, soilm_obs_mean ), dplyr::select( nice, soilm_mean ) ) )
        soilm_mod <- unlist( dplyr::select( nice, soilm_mean ) )
        soilm[ is.nan( soilm ) ] <- NA
        soilm_mod[ is.nan( soilm_mod ) ] <- NA
        out <- prune_droughts(  
                              droughts,   
                              nice$is_drought_byvar,
                              nice$fvar_smooth,
                              nice$fvar_smooth_filled,                                  
                              mod_pot   = nice$var_nn_pot, 
                              mod_act   = nice$var_nn_act, 
                              obs       = nice[[ nam_target ]],
                              soilm     = soilm,
                              soilm_mod = soilm_mod,
                              apar      = nice$ppfd * nice[[ fapar_data ]]
                              )
        droughts <- out$instances
        nice$is_drought_byvar_recalc <- out$is_drought_byvar

        ##------------------------------------------------
        ## get EVI anomalies
        ##------------------------------------------------
        if (do_evianom) out_evianomalies <- get_evianomalies_bysite( sitename, quantile=0.03, do.plot=FALSE ) 

        ##------------------------------------------------
        ## Determine wether algorithm failed
        ## Classify as failed depending on the following criteria:
        ## * during droughts NNgood is lower than NNall
        ## * R2 of NNall vs. observed
        ## * R2 of NNgood vs. observed during good days
        ## * RMSE of NNgood vs. NNall during good days
        ##------------------------------------------------
        # failed <- FALSE

        droughtdays <- nice %>% dplyr::filter(  is_drought_byvar )
        nondrgtdays <- nice %>% dplyr::filter( !is_drought_byvar )

        # ## * during droughts NNgood is lower than NNall
        # if ( nrow(droughtdays) > 0 ){
        #   if ( mean( droughtdays$var_nn_pot, na.rm=TRUE ) < mean( droughtdays$var_nn_act, na.rm=TRUE ) ) failed <- TRUE 
        # }

        ## * R2 of NNall
        mod <- nice$var_nn_act * nice$ppfd * nice$evi
        obs <- nice[[ nam_target ]] * nice$ppfd * nice$evi
        out_NNall <- analyse_modobs( 
                                    mod, 
                                    obs, 
                                    do.plot=TRUE, 
                                    plot.title=paste( "NN all", sitename ), 
                                    nrcol=1, 
                                    plot.fil=paste( workingdir, "/fig_nn_fluxnet2015/modobs_alldays/modobs_alldays", sitename, ".pdf", sep="") 
                                    )
        # if ( out_NNall$prmse > 60 || out_NNall$rsq < 0.3 ) failed <- TRUE

        ## * R2 of NNgood vs. observed during good days
        mod <- nondrgtdays$var_nn_pot * nondrgtdays$ppfd * nondrgtdays$evi
        obs <- nondrgtdays[[ nam_target ]] * nondrgtdays$ppfd * nondrgtdays$evi
        if (makepdf){
          plot.fil <- paste( workingdir, "/fig_nn_fluxnet2015/modobs_gooddays/modobs_gooddays", sitename, ".pdf", sep="") 
        } else {
          plot.fil <- NA
        }
        out_NNgood <- analyse_modobs( 
                                      mod, 
                                      obs, 
                                      do.plot=TRUE, 
                                      plot.title=paste( "LUE of NN_pot, moist days", sitename ), 
                                      nrcol=1, 
                                      plot.fil=plot.fil
                                      )
        # if ( out_NNgood$prmse > 60 || out_NNgood$rsq < 0.3 ) failed <- TRUE

        ## * R2 of NNall during bad days has no big bias
        if ( nrow(droughtdays) > 0 && nrow(droughtdays)>0.05*nrow(nondrgtdays) ){
          mod <- droughtdays$var_nn_act * droughtdays$ppfd * droughtdays$evi
          obs <- droughtdays[[ nam_target ]] * droughtdays$ppfd * droughtdays$evi
          if (makepdf){
            plot.fil <- paste( workingdir, "/fig_nn_fluxnet2015/modobs_alldays/modobs_alldays_bad", sitename, ".pdf", sep="") 
          } else {
            plot.fil <- NA
          }
          out_NNallbad <- analyse_modobs( 
                                          mod, 
                                          obs, 
                                          do.plot=TRUE, 
                                          plot.title=paste( "LUE of NN_act, all days", sitename ), 
                                          nrcol=1, 
                                          plot.fil=plot.fil
                                          )
          # if ( out_NNallbad$prmse > 70 || out_NNallbad$rsq < 0.3 ) failed <- TRUE
        } else {
          out_NNallbad <- NA
        }

        ## * RMSE of NNgood vs. NNall during good days
        if (makepdf){
          plot.fil <- paste( workingdir, "/fig_nn_fluxnet2015/pot_vs_act_NN/pot_vs_act_NN", sitename, ".pdf", sep="") 
        } else {
          plot.fil <- NA
        }
        out_NNgoodall <- analyse_modobs( 
                                        nondrgtdays$var_nn_act, 
                                        nondrgtdays$var_nn_pot, 
                                        do.plot=TRUE, 
                                        plot.title=paste( "LUE of NN_pot vs. NN_act, moist days", sitename ), 
                                        nrcol=1, 
                                        plot.fil=plot.fil
                                        )
        # if ( out_NNgoodall$prmse > 15 ) failed <- TRUE

        stats <- list( out_NNall=out_NNall, out_NNgood=out_NNgood, out_NNgoodall=out_NNgoodall, out_NNallbad=out_NNallbad )

        ##------------------------------------------------
        ## Save list of variables
        ##------------------------------------------------     
        nn_fluxnet[[ sitename ]]$nice             <- nice
        nn_fluxnet[[ sitename ]]$minmax           <- minmax
        nn_fluxnet[[ sitename ]]$droughts         <- droughts
        nn_fluxnet[[ sitename ]]$cutoff           <- cutoff
        nn_fluxnet[[ sitename ]]$varnams_swc      <- varnams_swc
        nn_fluxnet[[ sitename ]]$varnams_swc_obs  <- varnams_swc_obs
        nn_fluxnet[[ sitename ]]$stats            <- stats
        if (do_evianom) nn_fluxnet[[ sitename ]]$out_evianomalies <- NA

      } else {

        print("Drought identification failed.")
        nn_fluxnet[[ sitename ]]$nice             <- NA
        nn_fluxnet[[ sitename ]]$minmax           <- NA
        nn_fluxnet[[ sitename ]]$droughts         <- NA
        nn_fluxnet[[ sitename ]]$cutoff           <- NA
        nn_fluxnet[[ sitename ]]$varnams_swc      <- NA
        nn_fluxnet[[ sitename ]]$varnams_swc_obs  <- NA
        nn_fluxnet[[ sitename ]]$stats            <- NA
        if (do_evianom) nn_fluxnet[[ sitename ]]$out_evianomalies <- out_evianomalies
        
      }

      if (testprofile){
        outfil <- paste( workingdir, "/data/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
      } else {
        outfil <- paste( myhome, "/data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
      }

      print( paste( "writing file with fLUE and all other data into:", outfil ) )
      save( nn_fluxnet, file=outfil )

    }

  } else {

    print( paste( "No profile data available for site", sitename ) )

  }

}



