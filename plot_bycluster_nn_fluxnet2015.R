library(dplyr)
library(tidyr)
library(LSD)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )

## Manual settings-------
verbose = FALSE
makepdf = TRUE
##-----------------------

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

## Load variables from cluster_aligned_fluxnet2015.R: df_dday_agg, siteinfo_sub, mega, nclust, before, after
# load( "data/data_by_cluster_aligned.Rdata" )
# load( "data/data_by_cluster_fvar_vs_soilm.Rdata" )

## Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( "data/nice_agg_lue_obs_evi.Rdata" )   # loads 'nice_agg'
load( "data/nice_resh_lue_obs_evi.Rdata" )  # loads 'nice_resh'

load( "data/overview_data_fluxnet2015_L5.Rdata" )  # loads 'overview', written by cluster_fvar_vs_soilm.R or cluster_aligned_fluxnet2015.R

## Add cluster information to nice_agg
nice_agg  <- nice_agg %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster ) )
nice_resh <- nice_resh%>% left_join( dplyr::select( overview, mysitename, finalcluster ), by="mysitename" )

## Load aggregated data from aligned MODIS and MTE data
load( "data/data_aligned_agg.Rdata" )
df_dday_agg <- df_dday_agg %>% left_join( dplyr::select( overview, mysitename, finalcluster ), by="mysitename" )

before <- 30
after  <- 100

##------------------------------------------------
## Plot VPD vs. soil moisture
##------------------------------------------------
  magn <- 4
  ncols <- 2
  nrows <- 1
  heights <- 1*magn
  widths <- c(1.2, 0.8)*magn

  plotfiln <- paste( "./fig_nn_fluxnet2015/vpd_vs_soilm/vpd_vs_soilm_lue_obs_evi_ALL.pdf", sep="" )
  print(paste("plotting", plotfiln))
  if (makepdf) pdf( plotfiln, width=sum(widths), height=sum(heights) )      
    par( las=1, mar=c(4,4,1,1) )      
    panel <- layout(
                    matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                    # matrix( order, nrows, ncols, byrow=TRUE ),
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show( panel )

    with( 
          dplyr::filter( nice_resh, finalcluster %in% c(1,2) ), 
          heatscatter( 
                      soilm, 
                      vpd*1e-3, 
                      main="",
                      ylab="VPD (kPa)",
                      xlab="soil water content (fraction)",
                      xlim=c(0,1)
                      )
          )

    # ## draw the legend
    # legend( "topright", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", inset=c(0.05,0.25), cex=0.8 )

    boxplot( vpd*1e-3 ~ moist, data=dplyr::filter( nice_resh, finalcluster %in% c(1,2) ), outline=FALSE, col="grey70", xlab="moist", ylab="VPD (kPa)" )

    # vioplot( dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist )$vpd*1e-3, dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist )$vpd*1e-3, col="grey70" )

  if (makepdf) dev.off()


##------------------------------------------------
## fLUE vs. soilmoisture by cluster
##------------------------------------------------
  magn <- 3
  ncols <- 2
  nrows <- 3
  heights <- c(1,1,1.1)*magn
  widths  <- c(1.1,1)*magn*1.2

  plotfiln <- paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_agg_ALL.pdf", sep="" ) 
  print(paste("plotting", plotfiln))
  pdf( plotfiln, width=sum(widths), height=sum(heights) )

    panel <- layout(
                  matrix( c(1:(nrows*ncols)), nrows, ncols, byrow=TRUE ),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
    # layout.show( panel )

    ## histogram of values in bin
    soilm_lim <- 0.1
    lhist <- 100
    yhist <- with( dplyr::filter( nice_agg, finalcluster %in% c(1,2,3) & soilm_mean < soilm_lim & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    len <- length( yhist$counts )

    par( las=1, mar=c(2,4,1,1), xaxs="i", yaxs="i" )
    plot( c(0, max( yhist$counts ) ), c( 0, 1.4 ), type="n", xlim=c(-max(yhist$counts),0), axes=FALSE, xlab="", ylab="fLUE", cex.lab=1.2 )
    rect( -yhist$counts, yhist$breaks[1:len], rep( 0, len ), yhist$breaks[2:(len+1)], col=rgb(0,0,0,0.15) )
    axis( 2, lwd=1.5, cex.lab=2 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )

    yhist_1 <- with( dplyr::filter( nice_agg, finalcluster %in% c(1) & soilm_mean < soilm_lim & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_1$counts, yhist_1$breaks[1:len], rep( 0, len ), yhist_1$breaks[2:(len+1)], col=add_alpha("tomato", 0.5) )

    yhist_2 <- with( dplyr::filter( nice_agg, finalcluster %in% c(2) & soilm_mean < soilm_lim & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_2$counts, yhist_2$breaks[1:len], rep( 0, len ), yhist_2$breaks[2:(len+1)], col=add_alpha("royalblue3", 0.5) )

    yhist_3 <- with( dplyr::filter( nice_agg, finalcluster %in% c(3) & soilm_mean < soilm_lim & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_3$counts, yhist_3$breaks[1:len], rep( 0, len ), yhist_3$breaks[2:(len+1)], col=add_alpha("springgreen3", 0.5) )

    text( -320, 1.25, "a)", adj=c(0,0), font=2, cex=1.2 )

    legend( "bottomleft", legend=c("cDD", "cGR", "cLS"), fill=c( add_alpha("tomato", 0.6),  add_alpha("royalblue3", 0.6),  add_alpha("springgreen3", 0.6)), bty="n", inset=c(0,0) )


    ## clusters 1 2 3 pooled
    par( las=1, mar=c(2,2,1,1), xaxs="i", yaxs="i" )
    with( 
          dplyr::filter( nice_agg, finalcluster %in% c(1,2,3) ),
          heatscatter( 
                      soilm_mean, 
                      fvar,
                      ylim=c(0,1.4),
                      xlim=c(0,1),
                      main="",
                      xlab="soil water content (fraction)",
                      ylab="fLUE", 
                      cex.lab=1.2
                      ) 

        )
    box( lwd=1.5 )
    axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 2, lwd=1.5, cex.lab=2 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    text( 0.05, 1.25, "b)  cDD, cGR, and cLS pooled", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )
    abline( v=soilm_lim, lty=2 )


    ## cluster 1: big change in fAPAR
    par( las=1, mar=c(2,4,1,1), xaxs="i", yaxs="i" )
    with( 
          dplyr::filter( nice_agg, finalcluster==1 ),
          heatscatter( 
                      soilm_mean, 
                      fvar,
                      ylim=c(0,1.4),
                      xlim=c(0,1),
                      main="",
                      xlab="",
                      ylab="fLUE", 
                      cex.lab=1.2
                      ) 

        )
    box( lwd=1.5 )
    axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 2, lwd=1.5, cex.lab=2 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    text( 0.05, 1.25, "c)  cDD", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )

    ## cluster 2: small change in fAPAR
    par( mar=c(2,2,1,1) )
    with( 
          dplyr::filter( nice_agg, finalcluster==2 ),
          heatscatter( 
                      soilm_mean, 
                      fvar,
                      ylim=c(0,1.4),
                      xlim=c(0,1),
                      main="",
                      xlab="",
                      ylab=""
                      ) 

        )
    box( lwd=1.5 )
    axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 2, lwd=1.5 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5 ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    text( 0.05, 1.25, "d)  cGR", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )

    ## cluster 3: low LUE reduction with soil moisture (groundwater access?)
    par( mar=c(4,4,1,1) )
    with( 
          dplyr::filter( nice_agg, finalcluster==3 ),
          heatscatter( 
                      soilm_mean, 
                      fvar,
                      ylim=c(0,1.4),
                      xlim=c(0,1),
                      main="",
                      xlab="soil water content (fraction)",
                      ylab="fLUE",
                      cex.lab=1.2
                      ) 

        )
    box( lwd=1.5 )
    axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 2, lwd=1.5 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    text( 0.05, 1.25, "e)  cLS", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )

    ## cluster 4: no data in low soil moisture bin
    par( mar=c(4,2,1,1) )
    with( 
          dplyr::filter( nice_agg, finalcluster==4 ),
          heatscatter( 
                      soilm_mean, 
                      fvar,
                      ylim=c(0,1.4),
                      xlim=c(0,1),
                      main="",
                      xlab="soil water content (fraction)",
                      ylab="",
                      cex.lab=1.2
                      ) 

        )
    box( lwd=1.5 )
    axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 2, lwd=1.5 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
    axis( 4, lwd=1.5 ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
    text( 0.05, 1.25, "f)  cNA", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )

    legend( "bottomleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )

  dev.off()


##------------------------------------------------
## Aligned plots for all sites combined into clusters (3 rows)
##------------------------------------------------
  nclust <- 2
  magn <- 2
  ncols <- nclust
  nrows <- 3
  heights <- c(1,1,1.2)*magn
  widths  <- c(1.45,1.63)*magn

  plotfiln <- "fig_nn_fluxnet2015/aligned_cluster/aligned_cluster_3rows.pdf"
  print(paste("printing", plotfiln ))
  pdf( plotfiln, width=sum(widths), height=sum(heights) )

    panel <- layout(
                  matrix( c(1:(nrows*ncols)), nrows, ncols, byrow=FALSE ),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
    # layout.show( panel )

    for (iclust in 1:nclust){
      ##--------------------------------------------------------
      ## fLUE
      ##--------------------------------------------------------
      if (iclust==1) par( las=1, mar=c(0,4,1,1), xpd=FALSE, xaxs="i", yaxs="r" )
      if (iclust==2) par( las=1, mar=c(0,4,1,3), xpd=FALSE, xaxs="i", yaxs="r" )

      plot( c(-before,after), c(0,1.3), type="n", xlab="day after drought onset", ylab="fLUE", axes=FALSE )
      axis( 2 )
      if (iclust==1) axis( 4, labels=FALSE )
      if (iclust==2) axis( 4 )
      abline( h=1.0, col='grey40', lwd=0.5 )

      # # lines for each event
      # for (sitename in unique(df_dday_agg$mysitename)){
      #   if ( siteinfo_sub$alignedcluster[ which( siteinfo_sub$mysitename==sitename ) ]==iclust && siteinfo_sub$fvar_min[ which( siteinfo_sub$mysitename==sitename ) ]<0.9 ){
      #     sub <- dplyr::filter( df_dday_agg, mysitename==sitename )
      #     instances <- unlist( unique( ( dplyr::select( sub, inst ) ) ) )
      #     for (iinst in instances){
      #       subsub <- dplyr::filter( sub, inst==iinst )
      #       # rect( 0, -99, max( subsub$dday ), 99, col=rgb(0,0,0,0.002), border=NA )
      #       with( subsub, lines( dday, fvar, col=add_alpha( 'tomato', 0.1 ) ) )
      #     }
      #   }
      # }

      rect( 0, -99, after, 99, col=rgb(0,0,0,0.2), border=NA )

      ## filter data for this cluster
      sub <- df_dday_agg %>% dplyr::filter( finalcluster==iclust )

      ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      median <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = median( fvar , na.rm=TRUE ) ) %>% 
                        complete( dday )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.75, na.rm=TRUE ) ) %>% 
                        complete( dday )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.25, na.rm=TRUE ) ) %>% 
                        complete( dday )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.90, na.rm=TRUE ) ) %>% 
                        complete( dday )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.10, na.rm=TRUE ) ) %>% 
                        complete( dday )

      polygon( c( median$dday, rev(median$dday) ), c( lower$fvar,  rev(upper$fvar) ),  col=add_alpha("tomato", 0.3), border=NA )
      polygon( c( median$dday, rev(median$dday) ), c( llower$fvar, rev(uupper$fvar) ), col=add_alpha("tomato", 0.3), border=NA )
      lines( median, col='tomato2', lwd=2 )

      if (iclust==1) {
        text( 40, 1.2, "cDD", font=2, cex=1.2 )
      } else {
        text( 40, 1.2, "cGR", font=2, cex=1.2 )
      }


      ##--------------------------------------------------------
      ## fAPAR
      ##--------------------------------------------------------
      if (iclust==1) par( las=1, mar=c(0,4,0,1), xpd=FALSE, xaxs="i", yaxs="r" )
      if (iclust==2) par( las=1, mar=c(0,4,0,3), xpd=FALSE, xaxs="i", yaxs="r" )

      plot( c(-before,after), c(0,1.3), type="n", ylab=expression( paste( "EVI / EVI"[0] ) ), axes=FALSE, xlab="days after drought onset" )
      axis( 2 )
      if (iclust==1) axis( 4, labels=FALSE )
      if (iclust==2) axis( 4 )
      # axis( 1, xlab="days after drought onset" )

      # ## lines for each event
      # for (sitename in unique(df_dday_agg$mysitename)){
      #   if ( siteinfo_sub$alignedcluster[ which( siteinfo_sub$mysitename==sitename ) ]==iclust ){
      #     sub <- dplyr::filter( df_dday_agg, mysitename==sitename )
      #     instances <- unlist( unique( ( dplyr::select( sub, inst ) ) ) )
      #     for (iinst in instances){
      #       subsub <- dplyr::filter( sub, inst==iinst )
      #       with( subsub, lines( dday, evi_norm, col=add_alpha( 'springgreen3', 0.1 ) ) )
      #     }
      #   }
      # }

      rect( 0, -99, after, 99, col=rgb(0,0,0,0.2), border=NA )

      ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      median <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = median( evi_norm , na.rm=TRUE ) ) %>% 
                        complete( dday )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.75, na.rm=TRUE ) ) %>% 
                        complete( dday )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.25, na.rm=TRUE ) ) %>% 
                        complete( dday )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.90, na.rm=TRUE ) ) %>% 
                        complete( dday )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.10, na.rm=TRUE ) ) %>% 
                        complete( dday )

      polygon( c( median$dday, rev(median$dday) ), c( lower$evi_norm, rev(upper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
      polygon( c( median$dday, rev(median$dday) ), c( llower$evi_norm, rev(uupper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
      lines( median, col='springgreen4', lwd=2 )


      ##--------------------------------------------------------
      ## Soil moisture and VPD
      ##--------------------------------------------------------
      if (iclust==1) par( las=1, mar=c(4,4,0,1), xpd=FALSE, xaxs="i", yaxs="r" )
      if (iclust==2) par( las=1, mar=c(4,4,0,3), xpd=FALSE, xaxs="i", yaxs="r" )
      plot( c(-before,after), c(0,1.4), type="n", xlab="day after drought onset", ylab="relative change", axes=FALSE )
      axis( 2 )
      axis( 1, xlab="days after drought onset" )
      if (iclust==1) axis( 4, labels=FALSE )
      if (iclust==2) axis( 4 )
      axis( 1, xlab="days after drought onset" )

      # ##??lines for each event
      # for (sitename in unique(df_dday_agg$mysitename)){
      #   if ( siteinfo_sub$alignedcluster[ which( siteinfo_sub$mysitename==sitename ) ]==iclust ){
      #     sub <- dplyr::filter( df_dday_agg, mysitename==sitename )
      #     instances <- unlist( unique( ( select( sub, inst ) ) ) )
      #     for (iinst in instances){
      #       subsub <- dplyr::filter( sub, inst==iinst )
      #       with( subsub, lines( dday, bias_pmodel, col=add_alpha( 'royalblue3', 0.1 ) ) )
      #     }
      #   }
      # }

      rect( 0, -99, after, 99, col=rgb( 0, 0, 0, 0.2 ), border=NA )

      ## VPD change
      ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      median <- sub %>% group_by( dday ) %>% 
                        summarise( dvpd = median( 1/dvpd , na.rm=TRUE ) ) %>%
                        complete( dday )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( dvpd = quantile( 1/dvpd, 0.75, na.rm=TRUE ) ) %>%
                        complete( dday )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( dvpd = quantile( 1/dvpd, 0.90, na.rm=TRUE ) ) %>%
                        complete( dday )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( dvpd = quantile( 1/dvpd, 0.25, na.rm=TRUE ) ) %>%
                        complete( dday )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( dvpd = quantile( 1/dvpd, 0.10, na.rm=TRUE ) ) %>%
                        complete( dday )


      polygon( c( median$dday, rev(median$dday) ), c( lower$dvpd, rev(upper$dvpd) ), col=add_alpha("steelblue3", 0.4), border=NA )
      #polygon( c( median$dday, rev(median$dday) ), c( llower$bias_pmodel, rev(uupper$bias_pmodel) ), col=add_alpha("steelblue3", 0.4), border=NA )
      lines( median, col='steelblue3', lwd=2 )


      ## SOIL MOISTURE
      ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      median <- sub %>% group_by( dday ) %>% 
                        summarise( soilm_norm = median( soilm_norm , na.rm=TRUE ) ) %>%
                        complete( dday )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( soilm_norm = quantile( soilm_norm, 0.75, na.rm=TRUE ) ) %>%
                        complete( dday )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( soilm_norm = quantile( soilm_norm, 0.90, na.rm=TRUE ) ) %>%
                        complete( dday )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( soilm_norm = quantile( soilm_norm, 0.25, na.rm=TRUE ) ) %>%
                        complete( dday )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( soilm_norm = quantile( soilm_norm, 0.10, na.rm=TRUE ) ) %>%
                        complete( dday )

      polygon( c( median$dday, rev(median$dday) ), c( lower$soilm_norm, rev(upper$soilm_norm) ), col=add_alpha("royalblue3", 0.4), border=NA )
      #polygon( c( median$dday, rev(median$dday) ), c( llower$bias_pmodel, rev(uupper$bias_pmodel) ), col=add_alpha("royalblue3", 0.4), border=NA )
      lines( median, col='royalblue4', lwd=2 )

      legend( "topright", c("soil moisture (relative change)", expression( paste("VPD"^{-1}, " (relative change)")) ), bty="n", lty=1, lwd=2, col=c("royalblue4", "steelblue3"), cex=1.0 )


    }

  dev.off()

