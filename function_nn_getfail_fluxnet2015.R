nn_getfail_fluxnet <- function( sitename, recalc=TRUE, nam_target="lue_obs_evi", use_weights=ifelse( nam_target=="lue_obs_evi", TRUE, FALSE ), use_fapar=FALSE, testprofile=FALSE, makepdf=TRUE, verbose=FALSE ){

  # ## XXX debug----------------
  # sitename   = "FR-Pue"
  # nam_target = "lue_obs_evi"
  # use_weights= FALSE    
  # use_fapar  = FALSE
  # recalc     = TRUE
  # testprofile= TRUE
  # ##--------------------------

  require(dplyr)

  source( paste( myhome, "/analyse_modobs.R", sep="" ) )
  source( paste( myhome, "/plot_panel_nn.R", sep="" ) )

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

  filn  <- paste( "nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
  if (testprofile){
    dir <- "./data/"
  } else {
    dir <- paste( myhome, "data/nn_fluxnet/fvar/", sep="" )
  }

  avl <- TRUE
  failed_before <- FALSE
  path <- paste( dir, filn, sep="" )
  if (verbose) print( paste( "reading file", path ) )
  if (file.exists( path ) ){
    load( path )
  } else {
    avl <- FALSE
  }

  if (avl){

    ##------------------------------------------------
    ## "detatch"
    ##------------------------------------------------     
    nice <- nn_fluxnet[[ sitename ]]$nice             
    minmax <- nn_fluxnet[[ sitename ]]$minmax           
    droughts <- nn_fluxnet[[ sitename ]]$droughts         
    out_evianomalies <- nn_fluxnet[[ sitename ]]$out_evianomalies 
    cutoff <- nn_fluxnet[[ sitename ]]$cutoff           
    varnams_swc <- nn_fluxnet[[ sitename ]]$varnams_swc      
    varnams_swc_obs <- nn_fluxnet[[ sitename ]]$varnams_swc_obs  

    if (recalc){
      ##------------------------------------------------
      ## Determine wether algorithm failed
      ## Classify as failed depending on the following criteria:
      ## * during droughts NNgood is lower than NNall
      ## * R2 of NNall vs. observed
      ## * R2 of NNgood vs. observed during good days
      ## * RMSE of NNgood vs. NNall during good days
      ##------------------------------------------------
      droughtdays <- nice %>% dplyr::filter(  is_drought_byvar_recalc )
      nondrgtdays <- nice %>% dplyr::filter( !is_drought_byvar_recalc )

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
                                  plot.fil=paste( "fig_nn_fluxnet2015/modobs_alldays/modobs_alldays", sitename, ".pdf", sep="") 
                                  )
      # if ( out_NNall$prmse > 60 || out_NNall$rsq < 0.3 ) failed <- TRUE

      ## * R2 of NNgood vs. observed during good days
      mod <- nondrgtdays$var_nn_pot * nondrgtdays$ppfd * nondrgtdays$evi
      obs <- nondrgtdays[[ nam_target ]] * nondrgtdays$ppfd * nondrgtdays$evi
      out_NNgood <- analyse_modobs( 
                                    mod, 
                                    obs, 
                                    do.plot=TRUE, 
                                    plot.title=paste( "NN good, good days", sitename ), 
                                    nrcol=1, 
                                    plot.fil=paste( "fig_nn_fluxnet2015/modobs_gooddays/modobs_gooddays", sitename, ".pdf", sep="") 
                                    )
      # if ( out_NNgood$prmse > 60 || out_NNgood$rsq < 0.3 ) failed <- TRUE

      ## * R2 of NNall during bad days has no big bias
      if ( nrow(droughtdays) > 0 && nrow(droughtdays)>0.02*nrow(nondrgtdays) ){
        mod <- droughtdays$var_nn_act * droughtdays$ppfd * droughtdays$evi
        obs <- droughtdays[[ nam_target ]] * droughtdays$ppfd * droughtdays$evi
        out_NNallbad <- analyse_modobs( 
                                        mod, 
                                        obs, 
                                        do.plot=TRUE, 
                                        plot.title=paste( "NN all", sitename ), 
                                        nrcol=1, 
                                        plot.fil=paste( "fig_nn_fluxnet2015/modobs_alldays/modobs_alldays_bad", sitename, ".pdf", sep="") 
                                        )
        # if ( out_NNallbad$prmse > 70 || out_NNallbad$rsq < 0.3 ) failed <- TRUE
      } else {
        out_NNallbad <- list( prmse=NA, rsq=NA, ptoohi=NA, pbias=NA )
      }

      ## * RMSE of NNgood vs. NNall during good days
      out_NNNN <- analyse_modobs( 
                                  nondrgtdays$var_nn_act, 
                                  nondrgtdays$var_nn_pot, 
                                  do.plot=TRUE, 
                                  plot.title=paste( "NN pot vs NN act, good days", sitename ), 
                                  nrcol=1, 
                                  plot.fil=paste( "fig_nn_fluxnet2015/pot_vs_act_NN/pot_vs_act_NN", sitename, ".pdf", sep="") 
                                  )
      # if ( out_NNNN$prmse > 15 ) failed <- TRUE

      stats <- list( out_NNall=out_NNall, out_NNgood=out_NNgood, out_NNNN=out_NNNN, out_NNallbad=out_NNallbad )

    } else {

      stats <- nn_fluxnet[[ sitename ]]$stats            

    }


    ##------------------------------------------------
    ## Determine wether algorithm failed
    ## Classify as failed depending on the following criteria:
    ## * during droughts NNgood is lower than NNall
    ## * R2 of NNall vs. observed
    ## * R2 of NNgood vs. observed during good days
    ## * RMSE of NNgood vs. NNall during good days
    ##------------------------------------------------
    failed <- FALSE

    droughtdays <- nice %>% dplyr::filter(  is_drought_byvar_recalc )
    nondrgtdays <- nice %>% dplyr::filter( !is_drought_byvar_recalc )

    ## * during droughts NNgood is lower than NNall
    if ( nrow(droughtdays) > 0 ){
      if ( mean( droughtdays$var_nn_pot, na.rm=TRUE ) < mean( droughtdays$var_nn_act, na.rm=TRUE ) ) failed <- TRUE 
    }

    ## * R2 of NNall
    if ( stats$out_NNall$prmse > 60 || stats$out_NNall$rsq < 0.3 ) failed <- TRUE

    ## * R2 of NNgood vs. observed during good days
    if ( stats$out_NNgood$prmse > 60 || stats$out_NNgood$rsq < 0.3 ) failed <- TRUE

    ## * R2 of NNall during bad days has no big bias
    if ( nrow(droughtdays)>0.02*nrow(nondrgtdays) ){
      # if ( stats$out_NNallbad$ptoohi > 0.7 || stats$out_NNallbad$rsq < 0.2 ) failed <- TRUE
      if ( stats$out_NNallbad$rsq < 0.17 || stats$out_NNallbad$ptoohi > 0.88 || stats$out_NNallbad$prmse > 76 ) failed <- TRUE
    }

    ## * RMSE of NNgood vs. NNall during good days
    if ( stats$out_NNNN$rsq > 50  ) failed <- TRUE

    # ## delete panel file
    # dir <- paste( "/fig_nn_fluxnet2015/panel_potentialgpp/", sep="")
    # filn <- paste( "panel_crude_potentialgpp_", sitename, "_", nam_target, char_wgt, char_fapar, ".pdf", sep="" )
    # path <- paste( dir, filn, sep="" )
    # if (file.exists( path )) {
    #   system( paste( "rm", path ) )
    # } else {
    #   path <- paste( dir, "failed/", filn, sep="" )
    #   system( paste( "rm", path ) )
    # }

    # ## move NN fvar data file if it's re-classified now
    # dir <- paste( myhome, "data/nn_fluxnet/fvar/", sep="")
    # filn <- paste( "nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" )
    # path <- paste( dir, filn, sep="" )
    # if ( file.exists( path ) && failed ) {
    #   system( paste( "mv", path, paste( dir, "failed/" ) ) )
    # } else if ( !file.exists( path ) && !failed ){
    #   system( paste( "mv", paste( dir, "failed/", filn, sep="" ), dir ) )
    # }
    

    ##------------------------------------------------
    ## Return code for success
    ##------------------------------------------------
    if (failed){
      successcode <- 3
    } else if ( sum(droughts$len) < 0.02 * nrow(nice) ){
      successcode <- 2
    } else {
      successcode <- 1
    }

    out <- data.frame(
                      mysitename=sitename, successcode=successcode,

                      NNallbad_prmse=stats$out_NNallbad$prmse,
                      NNallbad_rsq=stats$out_NNallbad$rsq,
                      NNallbad_ptoohi=stats$out_NNallbad$ptoohi,
                      NNallbad_pbias=stats$out_NNallbad$pbias,

                      NNgood_prmse=stats$out_NNgood$prmse,
                      NNgood_rsq=stats$out_NNgood$rsq,
                      NNgood_ptoohi=stats$out_NNgood$ptoohi,
                      NNgood_pbias=stats$out_NNgood$pbias,

                      NNall_prmse=stats$out_NNall$prmse,
                      NNall_rsq=stats$out_NNall$rsq,
                      NNall_ptoohi=stats$out_NNall$ptoohi,
                      NNall_pbias=stats$out_NNall$pbias,

                      NNNN_prmse=stats$out_NNNN$prmse,
                      NNNN_rsq=stats$out_NNNN$rsq,
                      NNNN_ptoohi=stats$out_NNNN$ptoohi,
                      NNNN_pbias=stats$out_NNNN$pbias

                      )

    ## create new panel file and put to failed folder if necessary
    dir <- paste( "/fig_nn_fluxnet2015/panel_potentialgpp/", paste( "level_", as.character(successcode), sep="" ), sep="" )
    if (!file.exists( dir ) ) system( paste( "mkdir", dir ) )
    panelfiln <- paste( dir, "/panel_crude_potentialgpp_", sitename, "_", nam_target, char_wgt, char_fapar, ".pdf", sep="")
    # system( paste( "rm ", "/fig_nn_fluxnet2015/panel_potentialgpp/level_*/*", sep="" ) )

    ##------------------------------------------------
    ## plot time series with identified droughts
    ##------------------------------------------------
    plot_panel_nn( 
                  sitename, nice, out_nn_all$nn,
                  minmax         = minmax,
                  dolek          = FALSE, 
                  filename       = panelfiln,
                  droughts       = droughts,
                  cutoff         = cutoff,
                  fapar_extremes = out_evianomalies$extremes,
                  nam_target     = nam_target, 
                  makepdf        = makepdf
                  )


  } else {
    
    out <- data.frame(
                      mysitename=sitename, 
                      successcode=0,

                      NNallbad_prmse=NA,
                      NNallbad_rsq=NA,
                      NNallbad_ptoohi=NA,
                      NNallbad_pbias=NA,

                      NNgood_prmse=NA,
                      NNgood_rsq=NA,
                      NNgood_ptoohi=NA,
                      NNgood_pbias=NA,

                      NNall_prmse=NA,
                      NNall_rsq=NA,
                      NNall_ptoohi=NA,
                      NNall_pbias=NA,

                      NNNN_prmse=NA,
                      NNNN_rsq=NA,
                      NNNN_ptoohi=NA,
                      NNNN_pbias=NA

                      )
  }

  return( out )

}
