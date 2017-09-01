plot_aligned_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", bysm=FALSE, use_fapar=FALSE, use_weights=FALSE, makepdf=TRUE, verbose=FALSE, testprofile=FALSE ){

  # ## Debug ----------------
  # sitename   = "FR-Pue"
  # nam_target = "lue_obs_evi"
  # use_weights= FALSE
  # use_fapar  = FALSE
  # makepdf    = TRUE
  # verbose    = TRUE
  # testprofile= TRUE
  # ##--------------------------

  require(dplyr)
  
  source( "add_alpha.R" )

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

  siteinfo <- read.csv( "soilm_data_usability_fluxnet2015.csv", as.is=TRUE )

  ##------------------------------------------------
  ## Load data
  ##------------------------------------------------
  if (verbose) print("loading nn_fVAR file ...")
  if (testprofile) dir <- "data/" else dir <- paste( myhome, "/data/nn_fluxnet/fvar/", sep="" )
  infil <- paste( dir, "nn_fluxnet2015_", sitename, "_", nam_target, char_fapar, ".Rdata", sep="" ) 
  load( infil ) ## gets list 'nn_fluxnet'
  df <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )

  ## Get drought events. Default defined by fLUE, alternative (bysm=TRUE) defined by soil moisture threshold 0.5 (see 'nn_fVAR_fluxnet2015.R')
  if (bysm){
    ## Get soil moisture droughts
    print("get drought events ...")
    df$is_drought_bysoilm <- ( df$soilm_mean < 0.5 )
    droughts <- get_consecutive( 
                                df$is_drought_bysoilm, 
                                leng_threshold = 3, 
                                do_merge       = FALSE
                                )
  } else {
    droughts <- nn_fluxnet[[ sitename ]]$droughts
  }

  load( paste( "data/missing_pri_", nam_target, char_fapar, ".Rdata", sep="") )
  
  load( paste( "data/aligned_", sitename, char_bysm, ".Rdata", sep="" ) ) # loads data_alg_dry, names_alg, fvarbins, faparbins, iwuebins, before, after, bincentres_fvar, bincentres_fapar, bincentres_iwue
  load( paste( "data/df_dday_aggbydday_", sitename, char_bysm, ".Rdata", sep="" ) ) # loads 'df_dday_aggbydday'
  
  if ( is.element( sitename, missing_pri ) )      avl_pri <- FALSE else avl_pri <- TRUE
  if ( any(!is.na(df_dday_aggbydday$dscci_med)) ) avl_pri <- TRUE  else avl_pri <- FALSE

  if ( is.element( "pri", names(df)) ){
    df <- df %>% dplyr::select( year_dec, gpp_obs, var_nn_pot, var_nn_act, ppfd, fvar, soilm_mean, evi, fpar, wue_obs, is_drought_byvar, gpp_pmodel, gpp_obs_gfd, iwue, pri, cci, spri, scci )
  } else {
    df <- df %>% dplyr::select( year_dec, gpp_obs, var_nn_pot, var_nn_act, ppfd, fvar, soilm_mean, evi, fpar, wue_obs, is_drought_byvar, gpp_pmodel, gpp_obs_gfd, iwue )    
    avl_pri <- FALSE
  }

  filn <- paste( "data/aligned_modis_", sitename, ".Rdata", sep="" )
  if (file.exists(filn)) {
    error <- try( load( filn ) )   # loads data_alg_dry_modis, df_dday_modis, df_dday_aggbydday_modis, names_alg_modis, after_modis, before_modis
    if (class(error)=="try-error") avl_modis <- FALSE else avl_modis <- TRUE
  } else {
    avl_modis <- FALSE
  }

  filn <- paste( "data/aligned_mte_", sitename, ".Rdata", sep="" )
  if (file.exists(filn)) {
    error <- try( load( filn ) )   # loads data_alg_dry_mte, df_dday_mte, df_dday_aggbydday_mte, names_alg_mte, after_mte, before_mte
    if (class(error)=="try-error") avl_mte <- FALSE else avl_mte <- TRUE
  } else {
    avl_mte <- FALSE
  }


  ##------------------------------------------------
  ## PLOT MULTIPANEL
  ##------------------------------------------------
  if (verbose) print("plotting muiltipanel ...")
  lue <- TRUE
  panelfiln <- paste( "fig_nn_fluxnet2015/aligned/aligned_potentialgpp_", sitename, "_", nam_target, char_fapar, char_bysm, ".pdf", sep="")

  alpha <- 0.3/(nrow(droughts))
  xvals <- (-before:after)+1
 
  magn <- 2.4
  ncols <- 1
  # if (avl_pri) { nrows <- 4 } else { nrows <- 3 }
  nrows <- 2
  # heights <- rep(1,nrows)*magn
  heights <- c(1,1.4)*magn
  widths  <- rep(1.8,ncols)*magn

  if (makepdf) pdf( panelfiln, width=sum(widths), height=sum(heights), bg="white" )

    panel <- layout(
                    matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                    # matrix( order, nrows, ncols, byrow=TRUE ),
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show( panel )

    ##--------------------------------------------------------
    ## SOIL MOISTURE
    ##--------------------------------------------------------
      if (verbose) print("plot 1/5")
      
      xvals <- (-before:after)+1

      par( las=1, mar=c(0,4,1,4), xpd=FALSE )
      plot( c(-before,after), c(0,1.4), type="n", xlab="day after drought onset", ylab="soil water content (fraction)", axes=FALSE ) # , col.lab="blue"
      # axis( 2, col="blue", col.axis="blue" )
      axis( 2 )
      # axis( 4 )
      abline( h=1.0, col='grey40', lwd=0.5 )
      # mtext( "soil water content (fraction)", side=2, line=3, col="blue", cex=0.65, las=0 )

      rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

      upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$soilm_upp, xout=xvals )$y
      lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$soilm_low, xout=xvals )$y
      mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$soilm_med, xout=xvals )$y
      idxs  <- which( !is.na(upper) & !is.na(lower) )

      # ylim <- range( upper, lower, na.rm=TRUE )
      # par( las=1, mar=c(0,4,0,4), xpd=FALSE )
      # plot( (-before:after), (-before:after), type="n", xlab="day after drought", ylab"soil water content (fraction)", ylim=c(0,1), axes=FALSE )
      # axis( 2 )
      # axis( 4 )

      # boxplot( dplyr::select( dplyr::filter( df_dday, dday==1 ), soilm_mean ), at=1, add=TRUE, outline=FALSE, boxwex=10, axes=FALSE )

      # for ( idx in 1:nrow(droughts) ){
      #   lines( xvals, data_alg_dry[,which(names_alg=="soilm_mean"),idx], col=rgb(0,0,1,0.2) , lwd=1 ) #col=rgb(0,1,0,0.6)
      # }

      ## plot polygon for soil moisture
      polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('blue', 0.3), border=NA )
      lines( xvals, mid, col="blue", lwd=1 )

    ##--------------------------------------------------------
    ## Relative VPD change (in same panel as soil moisture)
    ##--------------------------------------------------------
      par( new=TRUE )
      xvals <- (-before:after)+1

      upper <- 1 / ( approx( df_dday_aggbydday$dday, df_dday_aggbydday$dvpd_upp, xout=xvals )$y ) 
      lower <- 1 / ( approx( df_dday_aggbydday$dday, df_dday_aggbydday$dvpd_low, xout=xvals )$y ) 
      mid   <- 1 / ( approx( df_dday_aggbydday$dday, df_dday_aggbydday$dvpd_med, xout=xvals )$y ) 
      idxs  <- which( !is.na(upper) & !is.na(lower) )

      # ylim <- range( upper, lower, na.rm=TRUE )
      # plot( c(-before,after), ylim, type="n", xlab="", ylab="", axes=FALSE )
      # axis( 4, col="darkgoldenrod4", col.axis="darkgoldenrod4" )
      # mtext( expression( paste( Delta, "VPD", "/VPD"[0] ) ), side=4, line=2.5, las=0, cex=0.7, col="darkgoldenrod4" )

      # for ( idx in 1:nrow(droughts) ){
      #   lines( xvals, data_alg_dry[,which(names_alg=="vpd"),idx], col=rgb(0,0,1,0.2) , lwd=1 ) #col=rgb(0,1,0,0.6)
      # }

      ## plot polygon for soil moisture
      polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('steelblue3', 0.3), border=NA )
      lines( xvals, mid, col=add_alpha('steelblue3', 1), lwd=1 )
      # try( lines( smooth.spline( xvals, mid, spar=0.3 ), col="steelblue3", lwd=1  ) )

      text( -before, 1.2, sitename, cex=1.2, adj=c(0,0), font=2 )
      # text( -before, 1.1, dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid, cex=1.2, adj=c(0,0), font=1 )

      legend( "bottomleft", c("soil water content (fraction)", expression( paste( "VPD"^{-1}, " (relative change)" ) ) ), bty="n", lty=1, lwd=2, col=c("blue", "steelblue3"), cex=1.0 )

    ##--------------------------------------------------------
    ## fLUE
    ##--------------------------------------------------------
      if (verbose) print("plot 2/5")
      par( las=1, mar=c(4,4,0,4), xpd=FALSE )
      plot( c(-before,after), c(0,1.3), type="n", xlab="day after drought onset", ylab="unitless", axes=FALSE )
      axis( 1, xlab="days after drought onset" )
      axis( 2 )
      # axis( 4 )
      abline( h=1.0, col='grey40', lwd=0.5 )
      # title( paste( sitename ) )

      # ## boxplot for levels within bins
      # bp1 <- boxplot( fvar ~ infvarbin, 
      #                 data    = df_dday, 
      #                 at      = bincentres_fvar, 
      #                 outline = FALSE, 
      #                 na.rm   = TRUE, 
      #                 add     = TRUE, 
      #                 axes    = FALSE, 
      #                 boxwex  = 5,
      #                 border  = "grey50",
      #                 lwd     = 0.5
      #               )

      ## GPP over drought aligned
      # for ( idx in 1:nrow(droughts) ){
      #   frac <- data_alg_dry[,which( names_alg=="fvar"),idx]
      #   lines( xvals, frac, col=rgb(0,0,0,0.2) )
      # }
      rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

      ## plot polygon for fvar
      upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_upp, xout=xvals )$y
      lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_low, xout=xvals )$y
      mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_med, xout=xvals )$y
      idxs  <- which( !is.na(upper) & !is.na(lower) )

      polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha("tomato", 0.4), border=NA )
      lines( xvals, mid, col="tomato" )

      # ## create new aligned data frame
      # df_aligned_stat <- data.frame( day=xvals, fvar=mid )

    ##--------------------------------------------------------
    ## FAPAR into same panel as above
    ##--------------------------------------------------------
      if (verbose) print("plot 3/5")
      xvals <- (-before:after)+1
      for ( idx in 1:nrow(droughts) ){
        lines( xvals, data_alg_dry[,which(names_alg==fapar_data),idx], col=add_alpha("springgreen4", 0.2) , lwd=1 ) #col=rgb(0,1,0,0.6)
      }

      # ## boxplot for levels within bins
      # bp1 <- boxplot( df_dday[[ fapar_data ]] ~ infaparbin, 
      #                 data    = df_dday, 
      #                 at      = bincentres_fapar, 
      #                 outline = FALSE, 
      #                 na.rm   = TRUE, 
      #                 add     = TRUE, 
      #                 axes    = FALSE, 
      #                 boxwex  = 10,
      #                 border  = "grey50",
      #                 lwd     = 0.5
      #               )

      ## plot polygon for EVI
      upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$evi_upp, xout=xvals )$y
      lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$evi_low, xout=xvals )$y
      mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$evi_med, xout=xvals )$y

      polygon( c( xvals, rev(xvals)), c( lower, rev(upper)), col=add_alpha('springgreen4', 0.3), border=NA )
      lines( xvals, mid, col="springgreen4", lwd=2 )

      # ## attach to aligned data frame
      # df_aligned_stat[[ fapar_data ]] <- mid

      legend( "left", c("fLUE", "EVI"), bty="n", lty=1, lwd=2, col=c("tomato", "springgreen4"), cex=1.0 )


    # ##--------------------------------------------------------
    # ## PRI and CCI
    # ##--------------------------------------------------------
    # if (avl_pri) {

    #   if (verbose) print("plot 4/5")
    #   par( las=1, mar=c(0,4,0,4), xpd=FALSE )
    #   plot( c(-before,after), c(-0.2,1.5), ylim=range( c( df_dday_aggbydday$dscci_upp, df_dday_aggbydday$dscci_low ), na.rm=TRUE ), type="n", xlab="day after drought onset", ylab=expression( paste( "CCI / CCI"[0] ) ), axes=FALSE )
    #   axis( 2 )
    #   axis( 4 )
    #   abline( h=1.0, col='grey40', lwd=0.5 )
    #   # title( paste( sitename ) )

    #   # ## boxplot for levels within bins
    #   # bp1 <- boxplot( pri ~ inpribin, 
    #   #                 data    = df_dday, 
    #   #                 at      = bincentres_pri, 
    #   #                 outline = FALSE, 
    #   #                 na.rm   = TRUE, 
    #   #                 add     = TRUE, 
    #   #                 axes    = FALSE, 
    #   #                 boxwex  = 5,
    #   #                 border  = "grey50",
    #   #                 lwd     = 0.5
    #   #               )

    #   # # GPP over drought aligned
    #   # for ( idx in 1:nrow(droughts) ){
    #   #   tmp <- data_alg_dry[,which( names_alg=="scci"),idx]
    #   #   tmp <- approx( xvals, tmp, xout=xvals )$y
    #   #   lines( xvals, tmp, col=add_alpha("goldenrod4", 0.4) )
    #   # }

    #   rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

    #   ## plot polygon for cci
    #   upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$dscci_upp, xout=xvals )$y
    #   lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$dscci_low, xout=xvals )$y
    #   mid   <- df_dday_aggbydday$dscci_med
    #   idxs  <- which( !is.na(upper) & !is.na(lower) )

    #   polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha("cadetblue4", 0.4), border=NA )
    #   lines( xvals, mid, col="cadetblue4" )

    #   # ## plot polygon for pri
    #   # upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$dspri_upp, xout=xvals )$y
    #   # lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$dspri_low, xout=xvals )$y
    #   # mid   <- df_dday_aggbydday$dspri_med
    #   # idxs  <- which( !is.na(upper) & !is.na(lower) )

    #   # polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha("cadetblue4", 0.4), border=NA )
    #   # lines( xvals, mid, col="cadetblue4" )

    #   # legend( "bottomleft", c("PRI", "CCI"), bty="n", lty=1, lwd=2, col=c('cadetblue4', "goldenrod4"), cex=1.0, inset=c(0,0.15) )
    #   legend( "bottomleft", c("CCI"), bty="n", lty=1, lwd=2, col=c("cadetblue4"), cex=1.0, inset=c(0,0.15) )

    # }


    # # ##--------------------------------------------------------
    # # ## IWUE
    # # ##--------------------------------------------------------
    # # if (verbose) print("plot 4/5")
    # # ndaysavg <- 10
    # # ref <- mean( data_alg_dry[ max( (before - ndaysavg), 0):before,which(names_alg=="iwue"),], na.rm=TRUE )
    # # xvals <- (-before:after)+1

    # # if (!is.nan(ref)){
    # #   upper <- apply( data_alg_dry[,which(names_alg=="iwue"),] / ref, 1, function(x) quantile(x,0.75, na.rm=TRUE)   )
    # #   upper <- approx( xvals, upper, xout=xvals )$y
    # #   lower <- apply( data_alg_dry[,which(names_alg=="iwue"),] / ref, 1, function(x) quantile(x,0.25, na.rm=TRUE)   )
    # #   lower <- approx( xvals, lower, xout=xvals )$y
    # #   mid   <- apply( data_alg_dry[,which(names_alg=="iwue"),] / ref, 1, function(x) quantile(x,0.5, na.rm=TRUE)   )
    # #   idxs  <- which( !is.na(upper) & !is.na(lower) )
      
    # #   ylim <- range( upper, lower, na.rm=TRUE )
    # #   ylim[1] <- max( 0.0, ylim[1] )

    # #   par( las=1, mar=c(0,4,0,4), xpd=FALSE )
    # #   plot( (-before:after), (-before:after), type="n", xlab="day after drought", ylab=expression( paste( Delta, "IWUE*" ) ), ylim=ylim, axes=FALSE )
    # #   if (verbose) print("y")
    # #   axis( 2 )
    # #   axis( 4 )
    # #   # for ( idx in 1:nrow(droughts) ){
    # #   #   lines( xvals, data_alg_dry[,which(names_alg=="wue_obs"),idx] / ref, col=rgb(0,0,0,0.2) , lwd=1 ) #col=rgb(0,1,0,0.6)
    # #   # }
    # #   rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

    # #   ## boxplot for levels within bins
    # #   bp1 <- boxplot( iwue / ref ~ iniwuebin, 
    # #                   data    = df_dday, 
    # #                   at      = bincentres_iwue, 
    # #                   outline = FALSE, 
    # #                   na.rm   = TRUE, 
    # #                   add     = TRUE, 
    # #                   axes    = FALSE, 
    # #                   boxwex  = 20,
    # #                   border  = "grey50",
    # #                   lwd     = 0.5
    # #                 )

    # #   ## plot polygon for WUE
    # #   polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('blue', 0.3), border=NA )
    # #   lines( xvals, mid, col="blue", lwd=1 )

    # #   # ## attach to aligned data frame
    # #   # df_aligned_stat$iwue <- mid

    # # } else {
    # #   plot( xvals, xvals, type="n", xlab="day after drought", ylab=expression( paste( Delta, "IWUE*" ) ), ylim=c(0,1), axes=FALSE )
    # #   axis( 2 )
    # #   axis( 4 )
    # #   rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )
    # #   # df_aligned_stat$iwue <- rep(NA, nrow(df_aligned_stat))
    # # }


    # ##--------------------------------------------------------
    # ## PLOT ARRANGED BIAS
    # ##--------------------------------------------------------
    # if (verbose) print("plot 5/5")

    # upper <- log( df_dday_aggbydday$bias_pmodel_upp )
    # lower <- log( df_dday_aggbydday$bias_pmodel_low )
    # mid   <- log( df_dday_aggbydday$bias_pmodel_med )
    # idxs  <- which( !is.na(upper) & !is.na(lower) )

    # if (avl_modis){
    #   upper_modis <- log( df_dday_aggbydday_modis$bias_modis_upp )
    #   lower_modis <- log( df_dday_aggbydday_modis$bias_modis_low )
    #   mid_modis   <- log( df_dday_aggbydday_modis$bias_modis_med )
    #   idxs_modis  <- which( !is.na(upper_modis) & !is.na(lower_modis) )  
    # }

    # if (avl_mte){
    #   upper_mte <- log( df_dday_aggbydday_mte$bias_mte_upp )
    #   lower_mte <- log( df_dday_aggbydday_mte$bias_mte_low )
    #   mid_mte   <- log( df_dday_aggbydday_mte$bias_mte_med )
    #   idxs_mte  <- which( !is.na(upper_mte) & !is.na(lower_mte) )
    # } else {
    #   upper_mte <- NA
    #   lower_mte <- NA
    #   mid_mte   <- NA
    # }

    # ylim <- quantile( c( upper, lower, lower_mte, upper_mte ), probs=c(0.02,0.98), na.rm=TRUE )
    # if (is.infinite(ylim[1])) ylim[1] <- -1
    # # ylim <- c(-1,2)
    # par( mar=c(4,4,0,4), xpd=FALSE )
    # plot( c(-before,after), c(0,3), type="n", xlab="day after drought onset", ylab="log(model/observed)", axes=FALSE, ylim=ylim )
    # axis( 2 )
    # axis( 4 )
    # axis( 1, xlab="days after drought onset" )

    # ## MODIS BIAS
    # if (avl_modis){
    #   xvals <- ((-before_modis:after_modis)+1)*8
    #   # for ( idx in 1:nrow(droughts_modis) ){
    #     # frac <- data_alg_dry_modis[,which( names(nice_to_modis)=="bias_modis"),idx]
    #     # lines( xvals, frac, col=rgb(0,0,0,0.2) )
    #   # }
    #   polygon( c( xvals[idxs_modis], rev(xvals[idxs_modis])), c( lower_modis[idxs_modis], rev(upper_modis[idxs_modis])), col=add_alpha("orchid", 0.4), border=NA )
    #   lines( xvals[idxs_modis], mid_modis[idxs_modis], col="orchid" )
    # }

    # ## MTE BIAS
    # if (avl_mte){
    #   xvals <- ((-before_mte:after_mte)+1)*8
    #   polygon( c( xvals[idxs_mte], rev(xvals[idxs_mte])), c( lower_mte[idxs_mte], rev(upper_mte[idxs_mte])), col=add_alpha("goldenrod3", 0.4), border=NA )
    #   lines( xvals[idxs_mte], mid_mte[idxs_mte], col="goldenrod3" )
    # }

    # ## PMODEL BIAS
    # xvals <- (-before:after)+1
    # # for ( idx in 1:nrow(droughts) ){
    #   # lines( xvals, data_alg_dry[,which(names_alg=="bias_pmodel"),idx], col=rgb(0,0,0,0.2) , lwd=1 ) #col=rgb(0,1,0,0.6)
    # # }
    # rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

    # polygon( c( xvals[idxs], rev(xvals[idxs])), c( (lower[idxs]), rev(upper[idxs])), col=add_alpha('royalblue3', 0.3), border=NA )
    # lines( xvals[idxs], mid[idxs], col="royalblue3", lwd=1 )

    # ## Legend
    # if (avl_mte){
    #   legend( "topleft", c("P-model", "MODIS", "FLUXCOM MTE"), bty="n", lty=1, lwd=2, col=c("royalblue3", "orchid", "goldenrod3"), cex=1.0 )
    # } else {
    #   legend( "topleft", c("P-model", "MODIS"), bty="n", lty=1, lwd=2, col=c("royalblue3", "orchid"), cex=1.0 )
    # }
    # abline( h=0.0, lwd=0.5  )

  if (makepdf) dev.off()

}


