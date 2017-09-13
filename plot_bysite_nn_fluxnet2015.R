lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

smooth_runminmax <- function( x, y ){  
  source( "cutna_headtail.R" )
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
  nice <- nice %>% left_join( tmp, by="year_dec" )
  return( nice$yvals )       
}


plot_bysite_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", use_fapar=FALSE, use_weights=FALSE, makepdf=TRUE, verbose=FALSE ){

  source("get_consecutive.R")

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

  panelfiln <- paste( "fig_nn_fluxnet2015/panel_potentialgpp/panel_potentialgpp_", sitename, "_", nam_target, char_wgt, char_fapar, ".pdf", sep="")   

  ## file name of input data file
  dir <- "./data/fvar/"
  infil <- paste( dir, "nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 

  ## this is necessary to avoid plotting into an already existing panel
  # if (!makepdf) dev.off()

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
    if (verbose) print( paste( "loading nn_fVAR file:", infil ) )
    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    minmax           <- nn_fluxnet[[ sitename ]]$minmax          
    droughts         <- nn_fluxnet[[ sitename ]]$droughts        
    out_evianomalies <- nn_fluxnet[[ sitename ]]$out_evianomalies; fapar_extremes <- out_evianomalies$extremes
    cutoff           <- nn_fluxnet[[ sitename ]]$cutoff
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

  ##------------------------------------------------
  ## PLOT MULTIPANEL
  ##------------------------------------------------
    if (verbose) print("plotting muiltipanel ...")

    lue <- TRUE
    nyears <- ceiling(range(nice$year_dec)[2]) - floor(range(nice$year_dec)[1])

    magn <- 3.5
    ncols <- 1

    ## with 2 panels only
    heights <- c(0.5,0.7)*magn
    nrows <- 2

    widths <- 10
    # widths  <- rep( nyears / 5,ncols) * magn

    if (makepdf) pdf( panelfiln, width=sum(widths), height=sum(heights), bg="white" )

      panel <- layout(
                      matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                      # matrix( order, nrows, ncols, byrow=TRUE ),
                      widths=widths,
                      heights=heights,
                      TRUE
                      )
      # layout.show( panel )

      library(stats)
      source( "add_alpha.R" )
      source( "loess_range.R" )

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
      # axis( 4 )
      abline( h=0, lwd=0.5 )
      # box()

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

      # ## removing interpolated NAs from plot again by white rectangle
      # rect( nice$year_dec[na_instances$idx_start], rep(0.01, nrow(na_instances)), nice$year_dec[na_instances$idx_start+na_instances$len-1], rep( 99, nrow(na_instances)), col="white", border=NA )

      ## observed
      if (lue){
        out <- smooth_runminmax( nice$year_dec, nice[[ nam_target ]] * nice$ppfd * nice[[ fapar_data ]] ) ## may have NAs at head or tail in 'nice[[ nam_target ]]'
      } else {
        out <- smooth_runminmax( nice$year_dec, nice[[ nam_target ]] )
      }
      polygon( c( out$x, rev(out$x) ), c( out$min, rev(out$max) ), border=NA, col=add_alpha("black",0.3) )
      lines(  out$x, out$mean, col="grey35", lwd=1 )

      ## removing interpolated NAs from plot again by white rectangle
      na_instances <- get_consecutive( is.na(nice$var_nn_pot), leng_threshold=5, do_merge=FALSE )
      rect( nice$year_dec[na_instances$idx_start], rep(0.01, nrow(na_instances)), nice$year_dec[na_instances$idx_start+na_instances$len-1], rep( 99, nrow(na_instances)), col="white", border=NA )

      ## label of site name and ecosystem type
      text( xlim[1], 0.9*ylim[2], sitename, cex=1.2, adj=c(0,0), font=2 )
      text( xlim[1], 0.75*ylim[2], dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid, cex=1.2, adj=c(0,0), font=1 )

      ## rectangles for droughts
      if (!is.null(droughts) && nrow(droughts)>0 ){
        rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
      }

      ## Legend
      legend( "topright", c("observed", expression( paste( "NN"[pot])), expression(paste("NN"[act])) ), bty="n", lty=1, lwd=2, col=c("grey35", cols[1], cols[2]), cex=0.9 )

      ##----------------------------------------------------------------------------------------
      ## 2nd panel: fLUE and EVI
      ##----------------------------------------------------------------------------------------
      if (verbose) print("plot 2/3")
      if (nrows==2){
        par( mar=c( 4, 4.4, 0, 3 ), xpd=FALSE )
      } else {
        par( mar=c( 0, 4.4, 0, 3 ), xpd=FALSE )
      }
      plot(  nice$year_dec, nice[[ nam_target ]],
        ylim=c(0,1.3), xlim=xlim,
        type="n", 
        xlab="year", ylab="unitless", axes=FALSE
        )
      if (nrows==2){
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )),     labels=FALSE )
        axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) ))+0.5, labels=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )), tck=0.0 )
      }
      axis( 2 )
      # axis( 4 )
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

      # ## Add rectangle for fAPAR extremes in time series
      # if (!is.null(fapar_extremes)){
      #   # abline( v=fapar_extremes$year_dec_start, col='red' )
      #   # abline( v=fapar_extremes$year_dec_end, col='red' )
      #   for (idx in 1:nrow(fapar_extremes)){
      #     idx_start <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_start[idx]) )
      #     idx_end   <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_end[idx]) )
      #     if (idx_start!=1 && idx_end!=1 ){
      #         rect( nice$year_dec[ idx_start ], ylim[1], nice$year_dec[ idx_end ], ylim[2], col=add_alpha('darkgoldenrod', 0.3), border=NA, lwd=2 )
      #     }
      #   }
      # }

      ## Legend
      # legend( "bottomright", c("fLUE", "MODIS EVI","MODIS FPAR"), bty="n", lty=1, lwd=c(1,2,1), col=c("tomato", "palegreen3", "springgreen3"), cex=1.2 )
      legend( "bottomright", c("fLUE", "MODIS EVI"), bty="n", lty=1, lwd=2, col=c("tomato", "springgreen3"), cex=0.9 )

    if (makepdf) dev.off()

}

