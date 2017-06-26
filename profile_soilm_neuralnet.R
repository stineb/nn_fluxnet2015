profile_soilm_neuralnet <- function( sitename, nam_target="gpp_obs", use_weights=ifelse( nam_target=="lue_obs_evi", TRUE, FALSE ), use_fapar=FALSE, varnams_swc=NA, soilm_threshold=NA, packages="nnet", overwrite_profile=TRUE, nrep=1, testprofile=FALSE, makepdf=TRUE ){

  # ## XXX debug----------------
  # sitename          = "FR-Pue"
  # nam_target        = "lue_obs_evi"
  # soilm_threshold   = 0.5
  # varnams_swc       = c("soilm_swbm")
  # packages          = "nnet"
  # nrep              = 2
  # overwrite_profile = TRUE  
  # use_fapar         = FALSE
  # makepdf           = FALSE
  # ##--------------------------

  ## default is sampling full range
  if (is.na(soilm_threshold)){
    soilm_threshold <- seq( 0.1, 0.60, 0.05 )
  }

  require( dplyr )
  require( abind )

  # source( "repeat_neuralnet.R" )
  source( paste( workingdir, "/analyse_modobs.R", sep="" ) )
  source( paste( workingdir, "/add_alpha.R", sep="" ) )
  source( paste( workingdir, "/predict_nn.R", sep="" ) )
  source( paste( workingdir, "/clean_fluxnet.R", sep="" ) )
  source( paste( workingdir, "/remove_outliers.R", sep="" ) )
  source( paste( workingdir, "/cleandata_nn.R", sep="" ) )

  load( paste( workingdir, "/data/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata", sep="" ) )

  packages <- "nnet"

  ## check and override if necessary
  if ( nam_target=="lue_obs" || nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
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

  if (use_weights){
    char_wgt <- "_wgt"
  } else {
    char_wgt <- ""
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

  if (testprofile){
    outdir <- "./data/profile/"
  } else {
    outdir <- paste( roothome, "/data/nn_fluxnet/profile/", sep="" )
  }

  if (!file.exists(outdir)) system( "mkdir -p ", outdir )

  outfilnam       <- paste( outdir, "profile_", nam_target, char_wgt, char_fapar, "_nn_", sitename, ".Rdata", sep="" )
  outfilnam_light <- paste( outdir, "profile_light_", nam_target, char_wgt, char_fapar, "_nn_", sitename, ".Rdata", sep="" )

  if ( !file.exists(outfilnam) || overwrite_profile ){

    if (is.na(varnams_swc)){
      varnams_swc <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob" )
    } else {
      print("USING THE FOLLOWING SOIL MOISTURE DATA:")
      print(varnams_swc)
    }

    ## Get neural network for this site, each target individually and looping over repetitions to get statistics
    print("=============================================")
    print( paste( "get profile for", sitename, "..." ) )

    profile_nn <- list()
    profile_nn_light <- list()
    rsq_all_by_smdata <- data.frame()
    rmse_good_by_smdata <- data.frame()

    ## Get model output from simulation without temperature or soil moisture limitation (s13)
    data <- fluxnet[[ sitename ]]$ddf$s13 ## wcont is SPLASH WITH 150 mm
    data <- subset( data, select=c(year_dec, wcont, gpp, aet, pet) )
    data <- data %>% dplyr::rename( soilm_splash150=wcont, gpp_pmodel=gpp, aet_pmodel=aet, pet_pmodel=pet )
    # data <- rename( data, c("wcont"="soilm", "gpp"="gpp_pmodel", "aet"="aet_pmodel", "pet"="pet_pmodel" ) )

    ## Get alternative soil moisture data
    data$soilm_splash220 <- fluxnet[[ sitename ]]$ddf$s11$wcont
    data$soilm_swbm      <- fluxnet[[ sitename ]]$ddf$s12$wcont
    data$soilm_etobs     <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et
    data$soilm_etobs_ob  <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et_orthbucket

    ## Get observational soil moisture data (variable availability!)
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

    ## Get observational data
    # obs <- subset( fluxnet[[ sitename ]]$ddf$obs, select=c( year_dec, gpp_obs2015_GPP_NT_VUT_REF, le_f_mds ) )
    obs <- dplyr::select( fluxnet[[ sitename ]]$ddf$obs, year_dec, gpp_obs2015_GPP_NT_VUT_REF, le_f_mds )

    ## Get input data
    inp <- fluxnet[[ sitename ]]$ddf$inp %>% dplyr::select( temp, ppfd, match( fapar_data, names(.) ), vpd, prec )

    if (dim(data)[1]==dim(obs)[1] && dim(data)[1]==dim(inp)[1] && data$year_dec[1]==obs$year_dec[1] && data$year_dec[dim(data)[1]]==obs$year_dec[dim(obs)[1]]){

      ## attach to data 
      data$gpp_obs <- obs$gpp_obs2015_GPP_NT_VUT_REF
      data$et_obs  <- obs$le_f_mds
      data$wue_obs <- data$gpp_obs / (data$et_obs*1e-6)
      data <- cbind( data, inp )

    }

    ##------------------------------------------------
    ## Remove outliers in WUE
    ##------------------------------------------------
    data$wue_obs <- remove_outliers( data$wue_obs, coef=3.0 )    

    ##------------------------------------------------
    ## get LUE and remove outliers
    ##------------------------------------------------
    if (!is.null(data$evi)){
      data <- data %>% mutate( lue_obs_evi  = remove_outliers( gpp_obs / ( ppfd * evi  ), coef=3.0 ) )
    }
    if (!is.null(data$fpar)){
      data <- data %>% mutate( lue_obs_fpar = remove_outliers( gpp_obs / ( ppfd * fpar ), coef=3.0 ) )      
    }

    ##------------------------------------------------
    ## normalise soil moisture
    ##------------------------------------------------
    data <- data %>% mutate( soilm_splash150 = soilm_splash150 / max( soilm_splash150, na.rm=TRUE ) )
    data <- data %>% mutate( soilm_splash220 = soilm_splash220 / max( soilm_splash220, na.rm=TRUE ) )
    data <- data %>% mutate( soilm_swbm      = soilm_swbm      / max( soilm_swbm     , na.rm=TRUE ) )
    data <- data %>% mutate( soilm_etobs     = soilm_etobs     / max( soilm_etobs    , na.rm=TRUE ) )
    data <- data %>% mutate( soilm_etobs_ob  = soilm_etobs_ob  / max( soilm_etobs_ob , na.rm=TRUE ) )
    for (ivar in varnams_swc_obs){
      data[[ ivar ]]  <- data[[ ivar ]]  / max( data[[ ivar ]] , na.rm=TRUE )
    }

    for ( isoilm_data in varnams_swc ){
      print("---------------------------------------------")
      print( paste( "soilm data source:", isoilm_data, "..." ) )

      ##------------------------------------------------
      ## Use only days where EVI is above 25% quantile of all days values and temperature is above 5 degrees
      ##------------------------------------------------
      df_nona <- data

      # df_nona <- dplyr::filter( df_nona, fapar > quantile( df_nona[[ fapar_data ]], probs=0.25 ) )
      df_nona <- dplyr::filter( df_nona, temp > 5.0 )

      ##------------------------------------------------
      ## do additional data cleaning, removing NAs, necessary for NN training
      ##------------------------------------------------
      df_nona <- cleandata_nn( df_nona, nam_target )

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
          df_nona <- df_nona[ which(!is.na(df_nona[[ivar]])), ]
        }

      }

      ## get weights for NN training: absorbed light
      if (use_weights) {
        weights <- df_nona$ppfd * df_nona$evi
      } else {
        weights <- rep( 1.0, nrow(df_nona) )
      }

      ##------------------------------------------------
      ## Get NN-predictions and its performance statistics for each soil moisture threshold
      ##------------------------------------------------
      for (isoilm_trh in soilm_threshold){

        ## subset of supposedly non-soil moisture limited days
        if ( isoilm_data=="soilm_obs" ){
          ## for observational data, do subset w.r.t. soil layer with highest value
          tmp <- subset( df_nona, select=varnams_swc_obs )
          tmp[ is.na(tmp) ] <- 0.0
          vec <- apply( tmp, 1, FUN = max, na.rm=TRUE )
        } else {
          vec <- df_nona[[ isoilm_data ]]
        }
        idxs_good <- which( vec > isoilm_trh )

        for (ipackage in packages){

          if ( length(idxs_good)>30 && (nrow(df_nona) - length(idxs_good))>30 ){
            ##------------------------------------------------
            ## Good days GPP (LUE)
            ##------------------------------------------------
            ## Get good-days-NN based on subset of data where soil moisture is above threshold
            # print("training good-days model ...")

            predictors <- c( "ppfd", "temp", "vpd" )
            if ( use_fapar ){ predictors <- c( predictors, fapar_data ) }

            ## Initialise
            var_nn_good <- array( NA, dim=c(length(idxs_good),nrep) )
            var_nn_bad  <- array( NA, dim=c(nrow(df_nona),nrep) )
            hidden_best <- NULL

            for (irep in 1:nrep){
              print(paste("NN repetition", irep))

              out_nn_good <- predict_nn( 
                data       = df_nona[ idxs_good, ],
                weights    = weights[ idxs_good ],
                predictors = predictors,
                nam_target = nam_target,
                do_predict = TRUE,
                package    = ipackage,
                lifesign   = "full",
                seed       = irep,
                hidden     = hidden_best
                )

              hidden_best <- out_nn_good$hidden_best

              ## keep only values
              var_nn_good[,irep] <- as.vector( out_nn_good$vals )
              
              ##------------------------------------------------
              ## Bad days GPP (LUE), predicted using the model trained on good days data
              ##------------------------------------------------
              # print("predicting bad days with good-days model ...")
              out_nn_bad <- predict_nn( 
                data       = df_nona, 
                weights    = weights,
                predictors = predictors, 
                nam_target = nam_target, 
                nn         = out_nn_good$nn, 
                do_predict = TRUE, 
                package    = ipackage
                )

              ## keep only values
              var_nn_bad[,irep] <- as.vector( out_nn_bad$vals )

            }     

            ## Evaluate predictions of good days model
            stats_nn_good <- analyse_modobs( 
                                            apply( var_nn_good, 1, FUN=mean ) * weights[ idxs_good ], 
                                            df_nona[[ nam_target ]][idxs_good], 
                                            plot.title=paste( "Good, trh=", isoilm_trh, "data=", isoilm_data ), 
                                            do.plot=FALSE
                                            )
  
            ## record number of hidden layers of best-performing NN (for quick NN training after this is known)
            profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$hidden_best_good <- hidden_best
            profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$hidden_best_good <- hidden_best

            ## Evaluate predictions of good days model during bad days
            stats_nn_bad <- analyse_modobs( 
                                            (apply( var_nn_bad, 1,  FUN=mean ) * weights)[-idxs_good], 
                                            df_nona[[ nam_target ]][-idxs_good], 
                                            plot.title=paste( "Bad, trh=", isoilm_trh, "data=", isoilm_data ), 
                                            do.plot=FALSE
                                            )

            print( paste( "R2 of NN-good with soil moisture threshold ", isoilm_trh, "is", stats_nn_good$rsq ) )


          } else {

            print(paste("too few data points with soil moisture threshold", isoilm_trh))

            idxs_good <- c()
            var_nn_good   <- NULL
            var_nn_bad    <- NULL
            stats_nn_bad  <- NULL
            stats_nn_good <- NULL

          }

          ## attach to nn_fuxnet list
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$idxs_good <- idxs_good

          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad  <- var_nn_bad
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_good <- var_nn_good

          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$stats_nn_bad  <- stats_nn_bad
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$stats_nn_good <- stats_nn_good

        }
      }

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

      ## Initialise
      var_nn_all <- array( NA, dim=c(nrow(df_nona),nrep) )
      hidden_best <- NULL

      for (irep in 1:nrep){
        print(paste("NN repetition", irep))

        out_nn_all <- predict_nn( 
          data       = df_nona, 
          weights    = weights,
          predictors = predictors, 
          nam_target = nam_target, 
          do_predict = TRUE, 
          package    = ipackage,
          seed       = irep,
          hidden     = hidden_best
          )

        hidden_best <- out_nn_good$hidden_best

        ## keep only values
        var_nn_all[,irep] <- as.vector( out_nn_all$vals )

      }      

      ## get statistics of mod vs. obs of all-days full model
      stats_nn_all <- analyse_modobs( 
                                      apply( var_nn_all, 1, FUN=mean ) * weights, 
                                      df_nona[[ nam_target ]], 
                                      plot.title=paste( "All, data=", isoilm_data ), 
                                      do.plot=FALSE
                                      )

      ##------------------------------------------------
      ## attach to nn_fuxnet list
      ##------------------------------------------------
      profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]]$var_nn_all      <- var_nn_all
      profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]]$stats_nn_all    <- stats_nn_all
      profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]] <- df_nona[[ nam_target ]]
      profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]]$year_dec        <- df_nona$year_dec
      profile_nn      [[ sitename ]][[ isoilm_data ]][[ ipackage ]]$hidden_best_all <- hidden_best      
      profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$hidden_best_all <- hidden_best      

      print( paste( "R2 of NN-all is", stats_nn_all$rsq ) )

    } 

    for ( isoilm_data in varnams_swc ){

      for (ipackage in packages){

        ##------------------------------------------------
        ## Reorganise data and  evaluate sensitivity to soil moisture threshold
        ##------------------------------------------------
        ## re-organise data into new data frames:
        ## df_good contains all predicted data of good days based on NN trained on "good days"

        soilm_thrsh_avl <- as.numeric(substr( ls(profile_nn[[sitename]][[isoilm_data]][[ipackage]]), start=7, stop=10 )[ is.element( substr( ls(profile_nn[[sitename]][[isoilm_data]][[ipackage]]), start=1, stop=6 ), "smtrh_" ) ])
        
        ##------------------------------------------------
        ## Do profiling only if it worked for more than zero soil moisture thresholds
        ##------------------------------------------------
        if (length(soilm_thrsh_avl)>0){

          ## df_good contains all predicted data of good days based on NN trained on "good days"
          df_good <- data.frame()
          for (isoilm_trh in soilm_thrsh_avl){
            if (!is.null(profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_good) &&
                !is.null(  profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad)){

              var_nn_good <- apply( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_good, 1, FUN=mean )
              idxs_good   <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$idxs_good

              tmp <- data.frame( 
                var_obs         = profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ idxs_good ],
                var_nn          = var_nn_good,
                ratio           = var_nn_good 
                                / profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ idxs_good ],  # ratio is mod / obs
                soilm_threshold = rep( isoilm_trh, length( idxs_good ) )
               )
              df_good <- rbind( df_good, tmp ) 

            }
          }

          df_good$goodbad <- rep("good", nrow(df_good))

          # ## get available (=tested) soil moisture levels
          # soilm_thrsh_avl <- unique( df_good$soilm_threshold )

          ## df_bad contains all predicted data of bad days based on NN trained on "good days"
          df_bad <- data.frame()
          for (isoilm_trh in soilm_threshold){
            if (!is.null(profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad)){

              idxs_good  <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$idxs_good
              var_nn_bad <- apply( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad, 1, FUN=mean )

              tmp <- data.frame( 
                var_obs         =profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ -idxs_good ],
                var_nn          =var_nn_bad[ -idxs_good ],
                ratio           =var_nn_bad[ -idxs_good ] 
                               / profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ -idxs_good ],  # ratio is mod / obs
                soilm_threshold =rep( isoilm_trh, length( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ -idxs_good ] ) )
               )
              df_bad <- rbind( df_bad, tmp ) 
            }
          }
          df_bad$goodbad <- rep("bad", nrow(df_bad))

          ##------------------------------------------------
          ## Determine N best soil moisture cutoffs w.r.t. best split based on difference in median of good-model vs. all-model predicted values of bad days
          ## criterium: median( NN-good(bad) / obs ) - median( NN-good(good) / obs ) = max!
          ##------------------------------------------------
          n_best <- 5
          diff_good_bad <- data.frame( soilm_threshold=c(), diff=c() )
          for (isoilm_trh in soilm_thrsh_avl){
            
            ## add row to statistics data frame
            addrow <-  data.frame( 
                                  soilm_threshold=isoilm_trh, 
                                  diff=median( dplyr::filter( df_bad, soilm_threshold==isoilm_trh )$ratio, na.rm=TRUE ) - median( dplyr::filter( df_good, soilm_threshold==isoilm_trh )$ratio, na.rm=TRUE ) 
                                  ) #, p_val=tmp_ks$p.value

            diff_good_bad <- rbind( diff_good_bad, addrow )

          }  
            
          list_best_soilm_trh <- diff_good_bad$soilm_threshold[ order(-diff_good_bad$diff)][1:min(n_best,nrow(diff_good_bad))]
          
          print("list_best_soilm_trh:")
          print( diff_good_bad[ order(-diff_good_bad$diff), ] )

          ##------------------------------------------------
          ## Determine best split based on variability of fVAR during good days
          ##------------------------------------------------
          ## re-organise statistics into new data frame
          df_stat_good <- data.frame()
          for (isoilm_trh in soilm_thrsh_avl){
            if (!is.null(profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$stats_nn_good)){
              idxs_good <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$idxs_good
              var_nn_bad <- apply( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad, 1, FUN=mean )
              fvar <- apply( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$var_nn_all / profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", isoilm_trh, sep="" ) ]]$var_nn_bad, 1, FUN=mean )
              mod <- var_nn_bad[ idxs_good ]
              obs <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]][ idxs_good ]
              linmod <- lm( mod ~ obs )
              tmp <- data.frame( 
                soilm_thrsh_avl = isoilm_trh,
                rmse            = rmse( obs, mod ),
                rsq             = summary( linmod )$r.squared,
                var             = var( fvar[idxs_good] )
               )
              df_stat_good <- rbind( df_stat_good, tmp ) 
            }
          }

          ## get soilm_thrsh_avl where best fit is achieved for good-days model (lowest variance of fVAR during good days)
          sub <- df_stat_good[ is.element( df_stat_good$soilm_thrsh_avl, list_best_soilm_trh ), ]
          best_soilm_trh <- sub$soilm_thrsh_avl[ order(sub$var)]

          ##------------------------------------------------
          ## Get RMSE-good for this dataset at best soil moisture threshold
          ##------------------------------------------------      
          stats_nn_good <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste( "smtrh_", best_soilm_trh[1], sep="" ) ]]$stats_nn_good
          stats_nn_all  <- profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$stats_nn_all

          rsq_all_by_smdata   <- rbind( rsq_all_by_smdata,   data.frame( smdata=isoilm_data, rsq_all=stats_nn_all$rsq ) )
          rmse_good_by_smdata <- rbind( rmse_good_by_smdata, data.frame( smdata=isoilm_data, rmse_good=stats_nn_good$rmse ) )

          ##------------------------------------------------
          ## save best soil moisture threshold (for this soil moisture dataset) and available soil moisture data sources
          ##------------------------------------------------
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$best_soilm_trh  <- best_soilm_trh
          profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$best_soilm_trh  <- best_soilm_trh
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$df_stat <- df_stat_good
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$diff_good_bad <- diff_good_bad

          ##------------------------------------------------
          ## Plot profile
          ##------------------------------------------------
          # library(sm)
          # library(vioplot)
          # vioplot( dplyr::filter( df_good, soilm_threshold=="0.01" )$ratio, col=c("springgreen"), ylim=c(0,3), las=1, xlab="soil moisture threshold", ylab="ratio of NN-modelled vs. observed GPP during good days", yaxs="i" )

          # maxy <- max( quantile( df_bad$ratio, probs=0.90 ), quantile( df_good$ratio, probs=0.90 ) )
          # maxy <- min( 10, max( max( df_bad$ratio ), max( df_good$ratio ) ) )
          maxy <- min( 10, max( quantile( df_bad$ratio, probs=0.90 ), quantile( df_good$ratio, probs=0.90 ) ) )
          if (makepdf) pdf( paste( workingdir, "/fig_nn_fluxnet2015/ratio_vs_threshold/ratio_vs_threshold_", nam_target, char_fapar, "_", isoilm_data, "_", sitename, ".pdf", sep="" ), width=8, height=6 )      
            bp1 <- boxplot( 
              ratio ~ soilm_threshold, 
              data =df_good, 
              col  =c("springgreen"), ylim=c(0,maxy), 
              # xlim=c(min(soilm_threshold), max(soilm_threshold)+0.025), 
              las=1, 
              xlab ="soil moisture threshold", 
              ylab ="ratio of NN-modelled vs. observed GPP", 
              yaxs ="i",
              at   =soilm_thrsh_avl-0.01,
              xlim =range(soilm_thrsh_avl)+c(-0.02,0.02),
              boxwex = 0.01,
              outline=FALSE
              )
            bp2 <- boxplot( 
              ratio ~ soilm_threshold, 
              data  =df_bad, 
              col   =c("grey70"), ylim=c(0,maxy), las=1, 
              xlab  ="", 
              ylab  ="", 
              yaxs  ="i",
              at    =soilm_thrsh_avl+0.01,
              # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
              boxwex = 0.01,
              add   = TRUE,
              axes  = FALSE,
              outline=FALSE
              )
            abline( h=1.0, col='red' ) 
            rect( soilm_thrsh_avl-0.02, rep(0, length(soilm_thrsh_avl)), soilm_thrsh_avl+0.02, rep(maxy,length(soilm_thrsh_avl)), 
              border=NA, col=c(rgb(0,0,0,0.2),rgb(0,0,0,0.0)) 
              )
            imax <- which( soilm_thrsh_avl==best_soilm_trh[1] )
            if (!is.na(best_soilm_trh[1])){
              rect( soilm_thrsh_avl[imax]-0.02, 0.01, soilm_thrsh_avl[imax]+0.02, maxy-0.01, 
                border='red', col=NA 
                )
            }
            title( sitename )
          if (makepdf) dev.off()


          ##------------------------------------------------
          ## Plot R-squared versus soil moisture threshold
          ##------------------------------------------------
          if (!is.na(best_soilm_trh[1])){
            if (makepdf) pdf( paste( workingdir, "/fig_nn_fluxnet2015/plot_rsq_vs_smtrh/plot_rsq_vs_smtrh_", nam_target, char_fapar, "_", isoilm_data, "_", sitename, ".pdf", sep="" ) )
              par( las=1, yaxs="i", xaxs="i" )
              plot( df_stat_good$soilm_thrsh_avl, df_stat_good$rsq, pch=16, xlab="soil moisture threshold (rel. units)", ylab=expression(paste("R"^2)), xlim=c(0,1), ylim=c(0,1), col='red' )
              abline( v=df_stat_good$soilm_thrsh_avl, col="grey70")
              abline( v=best_soilm_trh[1], col='green' )
            if (makepdf) dev.off()
          }

          ##------------------------------------------------
          ## Analyse modelled vs. observed on good days
          ##------------------------------------------------
          if (makepdf){ 
            plot.fil <- NA 
          } else {
            plot.fil <- paste( workingdir, "/fig_nn_fluxnet2015/modobs_gooddays/profile/modobs_gooddays_", nam_target, char_fapar, "_", isoilm_data, "_", sitename, ".pdf", sep="" )
          }
          if (!is.na(best_soilm_trh[1])){
            out <- analyse_modobs( 
                dplyr::filter( df_good, soilm_threshold==best_soilm_trh[1] )$var_nn, dplyr::filter( df_good, soilm_threshold==best_soilm_trh[1] )$var_obs,
                plot.xlab  = expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), 
                plot.ylab  = expression(paste("modelled GPP (gC m"^{-2}, " d"^{-1}, ")")), 
                plot.title = paste( sitename, ", good NN,", nam_target ), plot.col=rgb(0,0,0,0.3),
                plot.fil   = plot.fil
              ) 
          }

          ##------------------------------------------------
          ## Analyse modelled vs. observed on all days (full model)
          ##------------------------------------------------
          if (makepdf){ 
            plot.fil <- NA 
          } else {
            plot.dir <- paste( workingdir, "/fig_nn_fluxnet2015/modobs_alldays/profile/", sep="" )
            if (!file.exists(plot.dir)) system( paste( "mkdir -p", plot.dir ) )
            plot.fil <- paste( plot.dir, "modobs_alldays_", nam_target, char_fapar, "_", isoilm_data, "_", sitename, ".pdf", sep="" )
          }
          out <- analyse_modobs( 
              apply( profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$var_nn_all, 1, FUN=mean ), 
              profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ nam_target ]],
              plot.xlab  = expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), 
              plot.ylab  = expression(paste("modelled GPP (gC m"^{-2}, " d"^{-1}, ")")), 
              plot.title = paste( sitename, ", full NN,", nam_target ), plot.col=rgb(0,0,0,0.3),
              plot.fil   = plot.fil, 
              do.plot    = FALSE
            ) 


        } else {

          ##------------------------------------------------
          ## save best soil moisture threshold (for this soil moisture dataset) and available soil moisture data sources
          ##------------------------------------------------
          profile_nn[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$best_soilm_trh        <- NA
          profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$best_soilm_trh  <- NA

        }


      }

    }

    profile_nn[[ sitename ]]$varnams_swc           <- varnams_swc
    profile_nn[[ sitename ]]$varnams_swc_obs       <- varnams_swc_obs
    profile_nn_light[[ sitename ]]$varnams_swc     <- varnams_swc
    profile_nn_light[[ sitename ]]$varnams_swc_obs <- varnams_swc_obs

    ##------------------------------------------------
    ## Determine best soil moisture data w.r.t. RMSE of good days (given best soil moisture threshold)
    ##------------------------------------------------
    profile_nn      [[ sitename ]]$best_smdata <- as.character( rmse_good_by_smdata$smdata[ order( rmse_good_by_smdata$rmse_good ) ] )
    profile_nn_light[[ sitename ]]$best_smdata <- as.character( rmse_good_by_smdata$smdata[ order( rmse_good_by_smdata$rmse_good ) ] )

    ##------------------------------------------------
    ## Save per site
    ##------------------------------------------------
    print( paste( "saving profile file:", outfilnam ) )
    save( profile_nn, file=outfilnam )
    save( profile_nn_light, file=outfilnam_light )

  } else {

    print("=============================================")
    print( paste( "profile already available for site", sitename, "..." ) )

  }

}


