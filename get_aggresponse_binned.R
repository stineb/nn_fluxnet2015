get_aggresponse_binned <- function( sitename, nam_target="lue_obs_evi", use_fapar=FALSE, makepdf=TRUE, verbose=FALSE ){

  require(dplyr)

  source( "../add_alpha.R")

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

  if (verbose) print(paste("loading", sitename))

  ##------------------------------------------------
  ## load data
  ##------------------------------------------------
  load( paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_fapar, ".Rdata", sep="" )  ) ## gets list 'nn_fluxnet'
  df       <- as.data.frame( nn_fluxnet[[ sitename ]]$nice ) %>% dplyr::select( year_dec, match( fapar_data, names(.) ), fvar, wue_obs, is_drought_byvar, gpp_pmodel, gpp_obs, gpp_obs_gfd, iwue )
  droughts <- nn_fluxnet[[ sitename ]]$droughts        

  filn <- paste( "data/aligned_", sitename, ".Rdata", sep="" )
  if (file.exists(filn)){

    load( filn ) ## loads 'df_dday', defined in 'reshape_align_nn_fluxnet2015.R'

    nintervals_fvar  <- length(fvarbins)-1
    nintervals_fapar <- length(faparbins)-1
    nintervals_iwue  <- length(iwuebins)-1

    ##------------------------------------------------
    ## Aggregate fLUE
    ##------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    fvar <- df_dday %>% group_by( infvarbin ) %>% 
                        summarise( fvar  = median( fvar , na.rm=TRUE ) ) %>% 
                        complete( infvarbin, fill = list( fvar  = NA ) ) %>% 
                        dplyr::select( fvar )

    ## Add to array holding fvar response for each site (sites by row, intervals by column)
    fvar <- unlist( fvar )[1:nintervals_fvar]

    ##------------------------------------------------
    ## Aggregate FAPAR
    ##------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    fapar <- df_dday %>% group_by( infaparbin ) %>% 
                        summarise( fapar  = median( evi , na.rm=TRUE ) ) %>%   ## replace 'evi' with 'fpar' if required
                        complete( infaparbin, fill = list( fapar  = NA ) ) %>% 
                        dplyr::select( fapar )

    ## Add to array holding fapar response for each site (sites by row, intervals by column)
    fapar <- unlist( fapar )[1:nintervals_fapar]

    ## noramlise fapar response relative to median in first bin
    df_dday$evi_norm <- df_dday$evi / fapar[1]
    
    ## save this additional information (evi_norm)
    resave( df_dday, file=paste( "data/aligned_", sitename, ".Rdata", sep="" ) )

    fapar <- fapar / fapar[1]

    ##------------------------------------------------
    ## boxplot fvar vs. days (binned)
    ##------------------------------------------------
    filn <- paste( "/fig_nn_fluxnet2015/aligned_binned/aligned_binned_", sitename, "_", nam_target, char_fapar, ".pdf", sep="")
    if (makepdf) pdf( filn, width=8, height=6 )
    par(las=1)
    fvar_median <- df_dday %>% group_by( dday ) %>% summarise( fvar = median( fvar , na.rm=TRUE ) )
    plot( fvar_median, type='l', col='tomato', ylim = c(0,1.3) )
    abline( h=1.0, lwd=0.5 )
    bp1 <- boxplot( fvar ~ infvarbin, 
                    data    = df_dday, 
                    at      = bincentres_fvar, 
                    main    = sitename, 
                    las     = 1, 
                    outline = FALSE, 
                    na.rm   = TRUE, 
                    col     = add_alpha('tomato', 0.5), 
                    add     = TRUE, 
                    axes    = FALSE, 
                    boxwex  = 5
                  )
    points( bincentres_fvar, fvar, pch=16, col=add_alpha('red', 0.5) )

    fapar_median <- df_dday %>% group_by( dday ) %>% summarise( fapar = median( evi , na.rm=TRUE ) )
    lines( fapar_median, type='l', col='springgreen4' )
    bp1 <- boxplot( evi ~ infaparbin,            ## replace 'evi' with 'fpar' if required
                    data    = df_dday, 
                    at      = bincentres_fapar, 
                    las     = 1,
                    outline = FALSE, 
                    na.rm   = TRUE, 
                    col     = add_alpha('springgreen4', 0.5), 
                    add     = TRUE, 
                    axes    = FALSE, 
                    boxwex  = 5,
                    drop.unused.levels = TRUE
                  )
    abline( h=1.0, lwd=0.5 )  
    points( bincentres_fapar, fapar, pch=16, col=add_alpha('springgreen4', 1) )
    legend( "bottomright", c("fLUE", "EVI"), col=c('tomato', 'springgreen4'), lty=1, bty="n" )
    if (makepdf) dev.off()

  }  else {

    print("no aligned file available.")
    fvar    <- NA
    fapar   <- NA
    # iwue    <- NA
    df_dday <- NA

  }

  out <- list( sitename=sitename, fapar=fapar, fvar=fvar, df_dday=df_dday )  # , iwue=iwue
  return( out )

}
