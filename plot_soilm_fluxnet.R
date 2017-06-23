plot_soilm_fluxnet <- function( sitename, ddf, makepdf=FALSE ){

  #---------------------------------------------------------
  # Plot soil moisture calculated from different models and obs.
  #---------------------------------------------------------
  magn <- 6
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  heights <- rep(magn,nrows)
  widths[ncols] <- 0.4*magn
  order <- matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE)

  if (makepdf){ pdf( paste( "fig/soilm/soilm_fluxnet_", sitename, ".pdf", sep="" ), width=sum(widths), height=sum(heights) ) }

    ## observational data is given in volumentric soil water content (SWC), scale to within min and max to best compare against 
    ## SWC as a fraction of field capacity, quantified by the models.
    scale_obs <- TRUE

    panel <- layout(
                    order,
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show(panel)

    par( las=1, mar=c(5,4.5,2.5,0) )

    #---------------------------------------------------------
    # Time series
    #---------------------------------------------------------
    ylim <- c(0,1)
    plot(  ddf$obs$year_dec, ddf$s11$wcont / 220, type='l', ylim=ylim, ylab="soil water content (fraction)", xlab="year", main=sitename )
    lines( ddf$obs$year_dec, ddf$s12$wcont / 220, lty=2 )
    lines( ddf$obs$year_dec, ddf$swc_by_etobs$soilm_from_et / 220, col='red' )
    lines( ddf$obs$year_dec, ddf$swc_by_etobs$soilm_from_et_orthbucket / 220, col='red', lty=2 )
    # lines( ddf$obs$year_dec, soilm_from_et_keenan, col='red', lty=2 )

    if (scale_obs){

      plot_scaled <- TRUE
      if (ncol(ddf$swc_obs)>5){
        maxs       <- apply( ddf$swc_obs[,-(1:4)], 2, max, na.rm=TRUE ) 
        mins       <- apply( ddf$swc_obs[,-(1:4)], 2, min, na.rm=TRUE )
        maxs[]     <- max( maxs ) # use max/min across all depths for nomralisation 
        mins[]     <- min( mins ) # use max/min across all depths for nomralisation
      } else if (ncol(ddf$swc_obs)==5){
        maxs       <- max( ddf$swc_obs[,5], na.rm=TRUE )
        mins       <- min( ddf$swc_obs[,5], na.rm=TRUE )
      } else {
        plot_scaled <- FALSE
      }

      if (plot_scaled){

        ddf$swc_obs_scaled <- cbind( ddf$swc_obs[,(1:4)], as.data.frame( scale(  ddf$swc_obs[,-(1:4)], center = mins, scale = maxs - mins ) ))

        if (ncol(ddf$swc_obs)>4){
          for (icol in 5:ncol(ddf$swc_obs_scaled)){
            lines( ddf$obs$year_dec, ddf$swc_obs_scaled[,icol] , pch=16, col=rgb(0,0,1,0.5) ) 
          }
        }
      }

    } else {

      if (ncol(ddf$swc_obs)>4){
        for (icol in 5:ncol(ddf$swc_obs)){
          lines( ddf$obs$year_dec, ddf$swc_obs[,icol] / 100 , pch=16, col=rgb(0,0,1,0.5) ) 
        }
      }

    }

    par( xpd=TRUE )
    legend( min(ddf$obs$year_dec), -0.2, c("SLPASH", "SWBM", "bucket with ET data", "bucket with ET data and runoff > 0 before full", "measured at different depths"), 
      col=c("black", "black", "red", "red", rgb(0,0,1,0.8) ), lty=c(1,2,1,2,1), bty="n", cex=0.5
      )

    #---------------------------------------------------------
    # Distribution on the right 
    #---------------------------------------------------------
    par( las=1, mar=c(5,3,2.5,2), xpd=TRUE )
    lhist <- 50

    out_s11 <- ddf$s11$wcont / 220

    yhist    <- hist( out_s11, plot=FALSE, breaks=seq( from=0, to=1, length.out=lhist ) ) # note: this uses probability=TRUE
    out_dens <- density( out_s11, bw="SJ", adjust=5, na.rm=TRUE, from=0, to=1 )
    len <- length( yhist$counts )
    plot( out_dens$y, seq( from=0, to=lhist-1, length.out=length(out_dens$y) ), col='red', lwd=2, type='l', ylim=c(0,lhist), xlim=range(yhist$density), axes=FALSE, xlab="" ) # line
    axis( 2, at=seq( from=0, to=lhist, by=10), lab=as.character(seq( from=0, to=lhist, by=10)/lhist) )
    rect( rep( 0, len ), yhist$breaks[1:len]*lhist, yhist$density, yhist$breaks[2:(len+1)]*lhist, col=rgb(0,0,0,0.3) )
    abline( h=median(out_s11)*lhist, lwd=1, lty=1, col='springgreen4' )
   
  if (makepdf) { dev.off() }

}



