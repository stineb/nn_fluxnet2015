lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

smooth_runminmax <- function( x, y ){  
  source( paste( "/cutna_headtail.R", sep="" ) )
  idxs_drop <- cutna_headtail( y )
  if (length(idxs_drop)>0) { x <- x[-idxs_drop] }
  if (length(idxs_drop)>0) { y <- y[-idxs_drop] }
  x[is.nan(x)] <- NA
  x[is.infinite(x)] <- NA
  fld <- approx( x, y, xout=x )$y
  min <- smooth.spline( x, caTools::runmin( fld, k=7 ), spar=0.01 )$y
  max <- smooth.spline( x, caTools::runmax( fld, k=7 ), spar=0.01 )$y
  tmp <- approx( x, zoo::rollmean( fld, k=7, fill=NA ), xout=x )$y
  idxs<- which(!is.na(tmp))
  tmp <- smooth.spline( x[idxs], tmp[idxs], spar=0.01 )$y
  mean <- rep(NA, length(x))
  mean[idxs] <- tmp
  return( list( x=x, min=min, max=max, mean=mean ) )
}

spline_with_gaps <- function( xvals, yvals, nice ){
  idxs <- which( !is.na(yvals) & !is.nan(yvals)  )
  tmp <- smooth.spline( xvals[idxs], yvals[idxs] )
  tmp <- data.frame( year_dec=tmp$x, yvals=tmp$y )
  nice <- nice %>% left_join( tmp )
  return( nice$yvals )       
}


plot_bysite_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", use_fapar=FALSE, use_weights=FALSE, makepdf=TRUE, verbose=FALSE ){

  require( dplyr )
  require( tidyr )
  require( minpack.lm )
  require( LSD )
  require( cluster )
  require( broom )
  require( zoo )

  siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

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

  panelfiln <- paste( "/fig_nn_fluxnet2015/panel_potentialgpp/panel_potentialgpp_", sitename, "_", nam_target, char_wgt, char_fapar, ".pdf", sep="")   
  infil     <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 

  ## this is necessary to avoid plotting into an already existing panel
  # if (!makepdf) dev.off()

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
    if (verbose) print("loading nn_fVAR file ...")
    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    minmax           <- nn_fluxnet[[ sitename ]]$minmax          
    droughts         <- nn_fluxnet[[ sitename ]]$droughts        
    out_evianomalies <- nn_fluxnet[[ sitename ]]$out_evianomalies; fapar_extremes <- out_evianomalies$extremes
    cutoff           <- nn_fluxnet[[ sitename ]]$cutoff
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

  ##------------------------------------------------
  ## Determine data availability
  ##------------------------------------------------
  avl_pri <- with( nice, ifelse( is.element( "pri", names(nice) ), TRUE, FALSE ) )

  ## TEST
  avl_pri <- FALSE

  ##------------------------------------------------
  ## Get MTE-GPP for this site if available
  ##------------------------------------------------
  filn <- paste( "/data/mte_", sitename, ".Rdata", sep="" )
  if ( file.exists(filn) ){
    missing_mte <- FALSE
    load( filn )
  } else {
    missing_mte <- TRUE
  }

  ##------------------------------------------------
  ## Get MODIS-GPP for this site if available
  ##------------------------------------------------
  filn <- paste( "/data/modis_", sitename, ".Rdata", sep="" )
  if ( file.exists(filn) ){
    load( filn )
    avl_modisgpp <- TRUE
  } else {
    avl_modisgpp <- FALSE
  }

  ##------------------------------------------------
  ## Get MODIS-FPAR for this site
  ##------------------------------------------------
  filn <- paste( "/data/fapar_modis_", sitename, ".Rdata", sep="" )
  modis_fpar <- try( read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/sitedata/fapar/", sitename, "/dfapar_fpar_modissubset_", sitename, ".csv", sep="" ), as.is=TRUE ))
  if (class(modis_fpar)!="try-error" &&  sum(!is.na(modis_fpar$data))>2 ){
    modis_fpar <- dplyr::rename( modis_fpar, fapar_modis=data )
    modis_fpar$fapar_modis <- approx( modis_fpar$year_dec, modis_fpar$fapar_modis * 1e-1, xout=modis_fpar$year_dec )$y
    missing_modis_fpar <- FALSE
  } else {
    missing_modis_fpar <- TRUE
  }

  ##------------------------------------------------
  ## Get boxplot for PRI by drought
  ##------------------------------------------------
    if (avl_pri){
      if (makepdf) pdf( paste( "fig_nn_fluxnet2015/boxplot_pri/boxplot_pri_", sitename, ".pdf", sep="" ), width=5, height=6 )
        par(las=1)
        bp <- boxplot( 
          pri ~ is_drought_byvar, 
          # pri / (ppfd * evi) ~ is_drought_byvar, 
          data    = nice, 
          main    = paste( sitename ), 
          col     = c("royalblue3","tomato"), 
          las     = 1, 
          outline = FALSE,
          xlab    = "drought",
          ylab    = "PRI"
          )
      if (makepdf) dev.off()

    }

  ##------------------------------------------------
  ## PLOT MULTIPANEL
  ##------------------------------------------------
    if (verbose) print("plotting muiltipanel ...")
    
    ## OVERRIDE TO AVOID FOURTH PANEL
    # avl_pri <- FALSE

    lue <- TRUE
    nyears <- ceiling(range(nice$year_dec)[2]) - floor(range(nice$year_dec)[1])

    magn <- 3.5
    ncols <- 1
    # order <- c(1,1,2,2)
    if (avl_pri){
      nrows <- 4
    } else {
      nrows <- 3
    }
    heights <- rep(0.5,nrows)*magn
    widths <- 10
    # widths  <- rep( nyears / 5,ncols) * magn

    if (makepdf) pdf( panelfiln, width=sum(widths), height=sum(heights) )

      panel <- layout(
                      matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                      # matrix( order, nrows, ncols, byrow=TRUE ),
                      widths=widths,
                      heights=heights,
                      TRUE
                      )
      # layout.show( panel )

      library(stats)
      source( paste( myhome, "/add_alpha.R", sep="" ) )
      source( paste( myhome, "/loess_range.R", sep="" ) )

      cols <- c( "royalblue3", "tomato" )

      if (nam_target=="lue_obs_evi"||nam_target=="lue_obs_fpar"){
        varnam <- "fLUE"
        descr_outputset <- c(
          "NN, good days (LUE ~ temperature + VPD + PPFD)",
          "NN, all days (LUE ~ temperature + VPD + PPFD + soil moisture)"
          )
        lue <- TRUE
        if (nam_target=="lue_obs_evi"){
          data_fapar <- "evi"
        } else {
          data_fapar <- "fpar"
        }
      } else {
        varnam <- "fGPP"
        descr_outputset <- c(
          "NN, good days (GPP ~ temperature + VPD + PPFD)",
          "NN, all days (GPP ~ temperature + VPD + PPFD + soil moisture)"
          )
        lue <- FALSE
      }

      xlim <- range(nice$year_dec)

      ##----------------------------------------------------------------------------------------
      ## 1st panel: GPP absolute, obs, NN
      ##----------------------------------------------------------------------------------------
      if (verbose) print("plot 1/3")
      par( las=1, mar=c( 0, 4.4, 1, 3 ), xpd=FALSE )

      if (lue){
        ylim <- c( 0, max( nice[[ nam_target ]] * nice$ppfd * nice[[ fapar_data ]], nice$var_nn_pot  * nice$ppfd * nice[[ fapar_data ]], na.rm=TRUE ) )
      } else {
        ylim <- c( 0, max( nice[[ nam_target ]], nice$var_nn_pot, na.rm=TRUE ) )
      }
      plot(  nice$year_dec, nice[[ nam_target ]],
        ylim=ylim, xlim=xlim,
        type="n", 
        xlab="", ylab=expression(paste("GPP (gC m"^{-2}, " d"^{-1}, ")")), axes=FALSE
        )
      # axis( 1, at=unique(floor(nice$year_dec)),       labels=FALSE )
      # axis( 1, at=unique(floor(nice$year_dec))+0.5, , labels=unique(floor(nice$year_dec)), tck=0.0 )
      axis( 2 )
      axis( 4 )
      abline( h=0, lwd=0.5 )
      # box()

      text( xlim[1], 0.9*ylim[2], sitename, cex=1.5, adj=c(0,0), font=2 )
      text( xlim[1]+1.3, 0.9*ylim[2], dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid, cex=1.5, adj=c(0,0), font=1 )

      ## rectangles for droughts
      if (!is.null(droughts) && nrow(droughts)>0 ){
        rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
      }

      ## observed
      if (lue){
        out <- smooth_runminmax( nice$year_dec, nice[[ nam_target ]] * nice$ppfd * nice[[ fapar_data ]] ) ## may have NAs at head or tail in 'nice[[ nam_target ]]'
      } else {
        out <- smooth_runminmax( nice$year_dec, nice[[ nam_target ]] )
      }
      polygon( c( out$x, rev(out$x) ), c( out$min, rev(out$max) ), border=NA, col=add_alpha("black",0.3) )
      lines(  out$x, out$mean, col="grey35", lwd=1 )

      ## potential GPP = modelled without soil moisture limitation
      if (lue){
        out <- smooth_runminmax( nice$year_dec, nice$var_nn_pot * nice$ppfd * nice[[ fapar_data ]] )
      } else {
        out <- smooth_runminmax( nice$var_nn_pot )
      }
      polygon( c( out$x, rev(out$x) ), c( out$min, rev(out$max) ), border=NA, col=add_alpha(cols[1], 0.5) )
      lines(  out$x, out$mean, col=cols[1], lwd=1 )

      ## actual GPP
      if (lue){
        out <- smooth_runminmax( nice$year_dec, nice$var_nn_act * nice$ppfd * nice[[ fapar_data ]] )
      } else {
        out <- smooth_runminmax( nice$year_dec, nice$var_nn_act )
      }
      polygon( c( out$x, rev(out$x) ), c( out$min, rev(out$max) ), border=NA, col=add_alpha(cols[2], 0.5) )
      lines(  out$x, out$mean, col=cols[2], lwd=1 )

      ## Legend
      legend( "topright", c("observed", expression( paste( "NN"[pot])), expression(paste("NN"[act])) ), bty="n", lty=1, lwd=1, col=c("grey35", cols[1], cols[2]), cex=1.2 )

      ##----------------------------------------------------------------------------------------
      ## 2nd panel: fLUE and EVI
      ##----------------------------------------------------------------------------------------
      if (verbose) print("plot 2/3")
      par( mar=c( 0, 4.4, 0, 3 ), xpd=FALSE )
      plot(  nice$year_dec, nice[[ nam_target ]],
        ylim=c(0,1.3), xlim=xlim,
        type="n", 
        xlab="", ylab="unitless", axes=FALSE
        )
      # axis( 1, at=unique(floor(nice$year_dec)),       labels=FALSE )
      # axis( 1, at=unique(floor(nice$year_dec))+0.5, , labels=unique(floor(nice$year_dec)), tck=0.0 )
      axis( 2 )
      axis( 4 )
      abline( h=1, lwd=0.5 )
      # box()

      ## rectangles for droughts
      if (!is.null(droughts) && nrow(droughts)>0 ){
        rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
      }

      ## solid line for fvar_smooth
      lines( nice$year_dec, nice$fvar_smooth, col=cols[2], lwd=1 )

      ## Uncertainty range
      polygon( c( minmax$year_dec, rev(minmax$year_dec) ), c( smooth.spline( minmax$year_dec, minmax$fvar_min_filled, spar=0.01 )$y, rev( smooth.spline( minmax$year_dec, minmax$fvar_max_filled, spar=0.01 )$y ) ), border=NA, col=add_alpha(cols[2],0.5) )

      ## Plot time series of EVI
      lines( nice$year_dec, nice[[ fapar_data ]], col="springgreen3", lwd=2 )
      # if (!missing_modis_fpar){
      #   lines( modis_fpar$year_dec, modis_fpar$fapar_modis*1e1, col='springgreen3', lwd=1 )
      # }

      ## Add rectangle for fAPAR extremes in time series
      if (!is.null(fapar_extremes)){
        # abline( v=fapar_extremes$year_dec_start, col='red' )
        # abline( v=fapar_extremes$year_dec_end, col='red' )
        for (idx in 1:nrow(fapar_extremes)){
          idx_start <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_start[idx]) )
          idx_end   <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_end[idx]) )
          if (idx_start!=1 && idx_end!=1 ){
              rect( nice$year_dec[ idx_start ], ylim[1], nice$year_dec[ idx_end ], ylim[2], col=add_alpha('darkgoldenrod', 0.3), border=NA, lwd=2 )
          }
        }
      }

      ## Legend
      # legend( "bottomright", c("fLUE", "MODIS EVI","MODIS FPAR"), bty="n", lty=1, lwd=c(1,2,1), col=c("tomato", "palegreen3", "springgreen3"), cex=1.2 )
      legend( "bottomright", c("fLUE", "MODIS EVI"), bty="n", lty=1, lwd=c(1,2), col=c("tomato", "springgreen3"), cex=1.2 )


      ##----------------------------------------------------------------------------------------
      ## 3rd panel: GPP-MTE, MODIS GPP, and P-model vs. GPP-obs
      ##----------------------------------------------------------------------------------------
      nice$log_bias_pmodel <- log( nice$bias_pmodel )
      nice$log_bias_pmodel[ which(is.infinite(nice$log_bias_pmodel)) ] <- NA

      if (!missing_mte) { 
        nice_to_mte$log_bias_mte <- log( nice_to_mte$bias_mte )
        nice_to_mte$log_bias_mte[ which(is.infinite(nice_to_mte$log_bias_mte)) ] <- NA 

        # nice_to_mte$log_bias_rf <- log( nice_to_mte$bias_rf )
        # nice_to_mte$log_bias_rf[ which(is.infinite(nice_to_mte$log_bias_rf)) ] <- NA 

        # nice$log_bias_rf <- log( nice$bias_rf )
        # nice$log_bias_rf[ which(is.infinite(nice$log_bias_rf)) ] <- NA 
      }
      if (avl_modisgpp) { 
        nice_to_modis$log_bias_modis  <- log( nice_to_modis$bias_modis )
        nice_to_modis$log_bias_modis[ which(is.infinite(nice_to_modis$log_bias_modis)) ] <- NA 
      }

      if (!missing_mte) { 
        if (avl_modisgpp){
          ylim <- c( quantile( c(nice$log_bias_pmodel, nice_to_modis$log_bias_modis, nice_to_mte$log_bias_mte), 0.01, na.rm=TRUE ), quantile( c(nice$log_bias_pmodel, nice_to_modis$log_bias_modis, nice_to_mte$log_bias_mte), 0.99, na.rm=TRUE ) )
        } else {
          ylim <- c( quantile( c(nice$log_bias_pmodel, nice_to_mte$log_bias_mte), 0.01, na.rm=TRUE ), quantile( c(nice$log_bias_pmodel, nice_to_mte$log_bias_mte), 0.99, na.rm=TRUE ) )        
        }
      } else {
        if (avl_modisgpp){
          ylim <- c( quantile( c(nice$log_bias_pmodel, nice_to_modis$log_bias_modis), 0.01, na.rm=TRUE ), quantile( c(nice$log_bias_pmodel, nice_to_modis$log_bias_modis), 0.99, na.rm=TRUE ) )
        } else {
          ylim <- c( quantile( c(nice$log_bias_pmodel), 0.01, na.rm=TRUE ), quantile( c(nice$log_bias_pmodel), 0.99, na.rm=TRUE ) )          
        }
      }

      if (avl_pri){
        par( mar=c( 0, 4.4, 0, 3 ), xpd=FALSE )
      } else {
        par( mar=c( 4, 4.4, 0, 3 ), xpd=FALSE )
      }
      plot(  nice$year_dec, nice[[ nam_target ]],
        ylim=ylim, xlim=xlim,
        type="n", 
        xlab="year", ylab="log(model/observed)", axes=FALSE
        )
      if (!avl_pri){
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )),     labels=FALSE )
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) ))+0.5, labels=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )), tck=0.0 )
      }
      axis( 2 )
      axis( 4 )
      # box()


      ## rectangles for droughts
      if (!is.null(droughts) && nrow(droughts)>0 ){
        rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
      }

      abline( h=0, lwd=0.5 )

      ## bias of P-model GPP
      bias_pmodel_spl <- spline_with_gaps( nice$year_dec, nice$log_bias_pmodel, nice )
      lines( nice$year_dec, nice$log_bias_pmodel, col=add_alpha('royalblue3', 0.5), lwd=1 )
      lines( nice$year_dec, bias_pmodel_spl, col='royalblue3', lwd=2 )

      ## bias of MODIS GPP
      if (avl_modisgpp){
        lines( nice_to_modis$year_dec, nice_to_modis$log_bias_modis, col=add_alpha('springgreen4',0.5) )
        bias_modis_spl <- spline_with_gaps( nice_to_modis$year_dec, nice_to_modis$log_bias_modis, nice_to_modis )
        lines( nice_to_modis$year_dec, bias_modis_spl, col='springgreen4', lwd=2 )
      }

      ## bias of MTE-GPP
      if (!missing_mte) { 
        mycolnames <- "gpp_mte"
        for (ivar in mycolnames){
          idxs <- which( !is.na(nice_to_mte$log_bias_mte) & !is.nan(nice_to_mte$log_bias_mte) )
          # lines( nice_to_mte$year_dec, (nice_to_mte[[ ivar ]] - nice_to_mte$gpp_obs) / nice_to_mte$gpp_obs, col='tomato')
          if (length(idxs)>10){

            bias_mte_spl  <- spline_with_gaps( nice_to_mte$year_dec, nice_to_mte$log_bias_mte, nice_to_mte )
            # bias_8drf_spl <- spline_with_gaps( nice_to_mte$year_dec, nice_to_mte$log_bias_rf,  nice_to_mte )
            # bias_rf_spl   <- spline_with_gaps( nice$year_dec, nice$log_bias_rf, nice )

            ## 8-days MTE
            lines( nice_to_mte$year_dec, nice_to_mte$log_bias_mte, col=add_alpha('tomato', 0.5))
            lines( nice_to_mte$year_dec, bias_mte_spl, col='tomato', lwd=1.5 )

            # ## 8-days RF
            # lines( nice_to_mte$year_dec, nice_to_mte$log_bias_rf, col=add_alpha('tomato', 0.5))
            # lines( nice_to_mte$year_dec, bias_8drf_spl, col='springgreen4', lwd=1.5 )

            # ## daily RF
            # # lines( nice$year_dec, nice$log_bias_rf, col=add_alpha('tomato', 0.5))
            # lines( nice$year_dec, bias_rf_spl, col='royalblue3', lwd=1.5 )

          } else {
            missing_mte <- TRUE
          }
        }
      }

      if (!missing_mte) { 
        legend( "bottomright", c("P-model", "MODIS", "FLUXCOM MTE"), bty="n", lty=1, lwd=2, col=c("royalblue3", "springgreen4", "tomato"), cex=1.2 )
      } else {
        legend( "bottomright", c("P-model", "MODIS"), bty="n", lty=1, lwd=2, col=c("royalblue3", "springgreen4"), cex=1.2 )
      }

      ## Alternative legend for Tramontana data
      # legend( "bottomright", c("RF daily", "RF 8-days", "MTE 8-days"), bty="n", lty=1, lwd=2, col=c("royalblue3", "springgreen4", "tomato" ), cex=1.2 )

      ##----------------------------------------------------------------------------------------
      ## 4th panel: PRI
      ##----------------------------------------------------------------------------------------
      if (avl_pri){
        if (verbose) print("plot 4/3 - TESTING")
        par( mar=c( 4, 4.4, 0, 3 ), xpd=FALSE )
        with( nice, 
          plot( year_dec, scci, xlim=xlim, type='n', xlab="year", ylab="scaled PRI and CCI", axes=FALSE ) 
          )
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )),     labels=FALSE )
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) ))+0.5, labels=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )), tck=0.0 )
        axis( 2 )
        axis( 4 )

        with( nice, lines( year_dec, scci, col=add_alpha( 'goldenrod4', 0.5 ) ) )
        with( nice, lines( year_dec, spri, col=add_alpha( 'cadetblue4', 0.5 ) ) )

        ## spline curves
        cci_spl <- spline_with_gaps( nice$year_dec, nice$scci, nice )
        pri_spl <- spline_with_gaps( nice$year_dec, nice$spri, nice )

        lines( nice$year_dec, cci_spl, col='goldenrod4', lwd=2 )
        lines( nice$year_dec, pri_spl, col='cadetblue4', lwd=2 )

        ## rectangles for droughts
        if (!is.null(droughts) && nrow(droughts)>0 ){
          rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
        }

        # legend( "bottomright", c("scaled PRI10.12","scaled PRI13d"), bty="n", lty=1, lwd=1, col=c("goldenrod4", "cadetblue4"), cex=1.5 )
        legend( "bottomright", c("PRI", "CCI"), bty="n", lty=1, lwd=2, col=c('cadetblue4', "goldenrod4"), cex=1.5 )
      }

    if (makepdf) dev.off()

}

