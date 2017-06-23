plot_panel_nn <- function( sitename, df, minmax=NULL, nn, dolek, filename, droughts=NULL, markdays=NULL, markdays2=NULL, markdays3=NULL, cutoff=NULL, fapar_extremes=NULL, nam_target=NULL, makepdf=TRUE ){

  # library(NeuralNetTools)
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

  # print(paste("site", sitename))
  ##---------------------------------------------------------
  ## plot daily time series
  ##---------------------------------------------------------
  magn <- 3.5
  if (dolek){
    ncols <- 4
    order <- c(1,1,1,1,2,2,2,2)
  } else {
    ncols <- 2
    order <- c(1,1,2,2)
  }
  nrows <- 2
  heights <- c(1,1)*magn
  widths  <- rep(1.6,ncols)*magn

  if (makepdf) pdf( filename, width=sum(widths), height=sum(heights) )

    panel <- layout(
                    # matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                    matrix( order, nrows, ncols, byrow=TRUE ),
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show( panel )

    ##---------------------------------------------------------
    ## plot time series: observed and modelled GPP
    ##---------------------------------------------------------
    par( las=1, mar=c(4,4.4,3,2) )
    if (lue){
      ylim <- c( 0, max( df[[ nam_target ]], df$var_nn_pot  * df$ppfd * df[[ data_fapar ]], na.rm=TRUE ) )
    } else {
      ylim <- c( 0, max( df[[ nam_target ]], df$var_nn_pot, na.rm=TRUE ) )
    }
    ylab <- expression(paste("GPP (gC m"^{-2}, " d"^{-1}, ")"))
    # if (lue){
    #   ylab <- expression(paste("LUE (gC mol"^{-1}, ")"))
    # } else {
    #   ylab <- expression(paste("GPP (gC m"^{-2}, " d"^{-1}, ")"))
    # }
    plot(  df$year_dec, df[[ nam_target ]],
      ylim=ylim, 
      main=paste( sitename,sep=""),
      type="n", 
      xlab="year", ylab=ylab, axes=FALSE,
      )
    axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
    axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
    axis( 2 )
    box()

    ##----------------------------------------------------------------------------------------
    ## daily values
    ##----------------------------------------------------------------------------------------
    ## observed
    if (lue){
      lines(  df$year_dec, df[[ nam_target ]] * df$ppfd * df[[ data_fapar ]] , col=add_alpha("snow4", 1.0))
      # tmp <- approx( df$year_dec, df[[ nam_target ]] * df$ppfd * df[[ data_fapar ]], xout=df$year_dec)$y
      # idxs <- which( !is.na(tmp) )
      # lines( df$year_dec[idxs], runmed( tmp[idxs], 9), lwd=1 )
    } else {
      lines(  df$year_dec, df[[ nam_target ]], col=add_alpha("snow4", 1.0))      
    }

    ## potential GPP = modelled without soil moisture limitation
    if (lue){
      lines( df$year_dec, df$var_nn_pot * df$ppfd * df[[ data_fapar ]], col=add_alpha(cols[1], 1.0), lwd=1 )
      # tmp <- approx( df$year_dec, df$var_nn_pot * df$ppfd * df[[ data_fapar ]], xout=df$year_dec)$y
      # idxs <- which( !is.na(tmp) )
      # lines( df$year_dec[idxs], runmed( tmp[idxs], 9), col=cols[1], lwd=1 )
    } else {
      lines( df$year_dec, df$var_nn_pot, col=add_alpha(cols[1], 1.0), lwd=1 )
    }
    mtext( descr_outputset[1], col=cols[1], side=3, line=(1-1)*0.75, cex=0.6, adj=0 )

    ## actual GPP
    if (lue){
      lines( df$year_dec, df$var_nn_act * df$ppfd * df[[ data_fapar ]], col=add_alpha(cols[2], 1.0), lwd=1 )
      # tmp <- approx( df$year_dec, df$var_nn_act * df$ppfd * df[[ data_fapar ]], xout=df$year_dec)$y
      # idxs <- which( !is.na(tmp) )
      # lines( df$year_dec[idxs], runmed( tmp[idxs], 9), col=cols[2], lwd=1 )
    } else {
      lines( df$year_dec, df$var_nn_act, col=add_alpha(cols[2], 1.0), lwd=1 )
    }
    mtext( descr_outputset[2], col=cols[2], side=3, line=(2-1)*0.75, cex=0.6, adj=0 )



    ##----------------------------------------------------------------------------------------
    ## Grey bands for droughts
    ##----------------------------------------------------------------------------------------
    if (!is.null(droughts) && nrow(droughts)>0 ){
      rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
    }

    if (!is.null(fapar_extremes)){
      # abline( v=fapar_extremes$year_dec_start, col='red' )
      # abline( v=fapar_extremes$year_dec_end, col='red' )
      for (idx in 1:nrow(fapar_extremes)){
          idx_start <- which.min( abs(df$year_dec-fapar_extremes$year_dec_start[idx]) )
          idx_end   <- which.min( abs(df$year_dec-fapar_extremes$year_dec_end[idx]) )
          if (idx_start!=1 && idx_end!=1 ){
              rect( df$year_dec[ idx_start ], ylim[1], df$year_dec[ idx_end ], ylim[2], col=NA, border='red', lwd=2 )
          }
      }
    }


    ##---------------------------------------------------------
    ## plot time series: fGPP_filled or fLUE_filled (fraction of potential)
    ##---------------------------------------------------------
    ylim <- c( 0, 1.5 )
    par( las=1, mar=c(4,4.4,3,2) )
    plot(  df$year_dec, df[[ nam_target ]],
      ylim=ylim, 
      type="n", 
      xlab="year", ylab=varnam, axes=FALSE,
      )
    axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
    axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
    axis( 2 )
    box()

    ## solid line for fvar_smooth
    # lines( df$year_dec, df$fvar_filled, col=add_alpha(cols[2],0.3), lwd=1 )
    lines( df$year_dec, df$fvar_smooth, col=cols[2], lwd=1 )
    # lines( df$year_dec, df$fvar, col=cols[2], lwd=1 )

    ## Uncertainty range
    polygon( c( minmax$year_dec, rev(minmax$year_dec) ), c( smooth.spline( minmax$year_dec, minmax$fvar_min_filled, spar=0.01 )$y, rev( smooth.spline( minmax$year_dec, minmax$fvar_max_filled, spar=0.01 )$y ) ), border=NA, col=add_alpha(cols[2],0.5) )
    # polygon( c( minmax$year_dec, rev(minmax$year_dec) ), c( minmax$fvar_min_filled, rev( minmax$fvar_max_filled ) ), border=NA, col=add_alpha(cols[2],0.5) )

    # ## fvar based on observational soilmoisture
    # if (any(!is.na(df$fvar_obs))){
    #   lines( df$year_dec, df$fvar_obs, col=cols[2], lwd=1, lty=2 )
    # }

    # tmp <- approx( df$year_dec, df$fvar_filled, xout=df$year_dec)$y
    # idxs <- which( !is.na(tmp) )
    # lines( df$year_dec[idxs], runmed( tmp[idxs], 5), col=cols[2], lwd=1 )

    abline( h=1.0 )
    if (!is.null(cutoff)) { abline( h=cutoff, lty=2 )}
    mtext( varnam, col=cols[2], side=3, line=0, cex=0.6, adj=0 )

    if (!is.null(droughts) && nrow(droughts)>0 ){
      rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
    }

    ## Plot periods of EVI extremes as red boxes
    if (!is.null(fapar_extremes)){
      # abline( v=fapar_extremes$year_dec_start, col='red' )
      # abline( v=fapar_extremes$year_dec_end, col='red' )
      for (idx in 1:nrow(fapar_extremes)){
          idx_start <- which.min( abs(df$year_dec-fapar_extremes$year_dec_start[idx]) )
          idx_end   <- which.min( abs(df$year_dec-fapar_extremes$year_dec_end[idx]) )
          if (idx_start!=1 && idx_end!=1 ){
              rect( df$year_dec[ idx_start ], ylim[1], df$year_dec[ idx_end ], ylim[2], col=add_alpha('darkgoldenrod', 0.3), border=NA, lwd=2 )
          }
      }
    }

    ## Plot time series of EVI
    lines( df$year_dec, df[[ data_fapar ]], col="palegreen3", lwd=2 )


    if (is.element("fet", names(data))){
      ##---------------------------------------------------------
      ## plot time series: observed and modelled ET
      ##---------------------------------------------------------
      par( las=1, mar=c(4,4.4,3,2) )
      plot(  df$year_dec, df$et_obs*1e-6,
        ylim=c( 0, max( df$et_obs*1e-6, df$et_nn_pot*1e-6, na.rm=TRUE ) ), 
        main=paste( sitename,sep=""),
        type="n", 
        xlab="year", ylab=expression(paste("ET (MJ m"^{-2}, " d"^{-1}, ")")), axes=FALSE,
        )
      axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
      axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
      axis( 2 )
      box()

      ## observed
      lines(  df$year_dec, df$et_obs*1e-6, col=add_alpha("snow4", 0.5))

      ## potential et = modelled without temp or soil moisture limitation
      # if (!is.null(df$et_nn_pot)){
        lines( df$year_dec, df$et_nn_pot*1e-6, col=cols[1], lwd=1 )
        mtext( descr_outputset[1], col=cols[1], side=3, line=(1-1)*0.75, cex=0.6, adj=0 )
      # }

      ## NN-predicted et
      lines( df$year_dec, df$et_nn_act*1e-6, col=cols[2], lwd=1 )
      mtext( descr_outputset[2], col=cols[2], side=3, line=(2-1)*0.75, cex=0.6, adj=0 )

      if (!is.null(droughts)&& nrow(droughts)>0 ){
        rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(1,0,0,0.2), border=NA )
      }


      ##---------------------------------------------------------
      ## plot time series: Normalised ET (fraction of potential)
      ##---------------------------------------------------------
      par( las=1, mar=c(4,4.4,3,2) )
      plot(  df$year_dec, df$et_obs,
        ylim=c( 0, 1.5 ), 
        type="n", 
        xlab="year", ylab=expression(paste("act.-et / pot.-et")), axes=FALSE,
        )
      axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
      axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
      axis( 2 )
      box()

      ## NN-predicted et
      lines( df$year_dec, df$fet_filled, col=add_alpha(cols[2],0.3), lwd=1 )
      lines( df$year_dec, df$fet_smooth, col=cols[2], lwd=1 )

      abline( h=1.0 )
      if (!is.null(cutoff)) { abline( h=cutoff, lty=2 )}
      mtext( "act.-ET / pot.-ET", col=cols[2], side=3, line=0, cex=0.6, adj=0 )

      if (!is.null(droughts)&& nrow(droughts)>0 ){
        rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(1,0,0,0.2), border=NA )
      }

    }

    if (is.element("wue_obs", names(data)) && is.element("wue_spline", names(data))){
      ##---------------------------------------------------------
      ## plot time series: observed and modelled wue
      ##---------------------------------------------------------
      par( las=1, mar=c(4,4.4,3,2) )
      plot(  df$year_dec, df$wue_spline,
        ylim=c( 0, max( df$wue_spline, na.rm=TRUE ) ), 
        main=paste( sitename,sep=""),
        type="n", 
        xlab="year", ylab=expression(paste("WUE (gC MJ"^{-1}, ")")), axes=FALSE,
        )
      axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
      axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
      axis( 2 )
      box()

      ## observed
      lines(  df$year_dec, df$wue_obs, col=add_alpha(cols[1],0.3) )

      ## splined observed
      lines(  df$year_dec, df$wue_spline, col=cols[1], lwd=2 )

      # ## potential wue = modelled without temp or soil moisture limitation
      # # if (!is.null(df$wue_nn_pot)){
      #   lines( df$year_dec, df$wue_nn_pot, col=cols[1], lwd=1 )
      #   mtext( descr_outputset[1], col=cols[1], side=3, line=(1-1)*0.75, cex=0.6, adj=0 )
      # # }

      # ## NN-predicted et
      # lines( df$year_dec, df$wue_nn_act, col=cols[2], lwd=1 )
      # mtext( descr_outputset[2], col=cols[2], side=3, line=(2-1)*0.75, cex=0.6, adj=0 )

      if (!is.null(droughts)&& nrow(droughts)>0 ){
        rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(1,0,0,0.2), border=NA )
      }


      if (is.element("fwue", names(data))){
        ##---------------------------------------------------------
        ## plot time series: Normalised wue (fraction of potential)
        ##---------------------------------------------------------
        par( las=1, mar=c(4,4.4,3,2) )
        plot(  df$year_dec, df$fwue,
          ylim=c( 0, 1.5 ), 
          type="n", 
          xlab="year", ylab=expression(paste("fWUE")), axes=FALSE,
          )
        axis( 1, at=unique(floor(df$year_dec)),       labels=FALSE )
        axis( 1, at=unique(floor(df$year_dec))+0.5, , labels=unique(floor(df$year_dec)), tck=0.0 )
        axis( 2 )
        box()

        ## NN-predicted wue
        # lines( df$year_dec, df$fwue_filled, col=add_alpha(cols[2],0.3), lwd=1 )
        lines( df$year_dec, df$fwue_smooth, col=add_alpha(cols[2],0.3), lwd=1 )
        lines( df$year_dec, df$fwue_spline, col=cols[2], lwd=2 )

        abline( h=1.0 )
        if (!is.null(cutoff)) { abline( h=cutoff, lty=2 )}
        mtext( "fWUE", col=cols[2], side=3, line=0, cex=0.6, adj=0 )

        if (!is.null(droughts)&& nrow(droughts)>0 ){
          rect( df$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), df$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(1,0,0,0.2), border=NA )
        }
      }

    }

  if (makepdf) dev.off()

}

