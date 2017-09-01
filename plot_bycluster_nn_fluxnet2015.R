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
load( "data/nice_agg_lue_obs_evi.Rdata" )       # loads 'nice_agg'

load( "data/overview_data_fluxnet2015_L5.Rdata" )  # loads 'overview', written by cluster_fvar_vs_soilm.R or cluster_aligned_fluxnet2015.R

## Add cluster information to nice_agg
nice_agg          <- nice_agg          %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster ) )


## Load aggregated data from aligned MODIS and MTE data
load( "data/data_aligned_agg.Rdata" )
df_dday_agg       <- df_dday_agg       %>% left_join( dplyr::select( overview, mysitename, finalcluster ), by="mysitename" )
# df_dday_modis_agg <- df_dday_modis_agg %>% left_join( dplyr::select( overview, mysitename, finalcluster  ), by="mysitename" )
# df_dday_mte_agg   <- df_dday_mte_agg   %>% left_join( dplyr::select( overview, mysitename, finalcluster  ), by="mysitename" )


##------------------------------------------------
## Plot VPD vs. soil moisture
##------------------------------------------------
  magn <- 4
  ncols <- 2
  nrows <- 1
  heights <- 1*magn
  widths <- c(1.2, 0.8)*magn

  if (makepdf) pdf( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/fig_nn_fluxnet2015/vpd_vs_soilm/vpd_vs_soilm_lue_obs_evi_ALL.pdf", sep="" ), width=sum(widths), height=sum(heights) )      
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
  nice_agg$ratio_obs_mod <- remove_outliers( nice_agg$ratio_obs_mod, coef=5 )

  magn <- 3
  ncols <- 2
  nrows <- 3
  heights <- c(1,1,1.1)*magn
  widths  <- c(1.1,1)*magn*1.2

  pdf( paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_agg_ALL.pdf", sep="" ), width=sum(widths), height=sum(heights) )

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
## fLUE vs. soilmoisture by cluster
##------------------------------------------------
  nice_agg$ratio_obs_mod <- remove_outliers( nice_agg$ratio_obs_mod, coef=5 )

  magn <- 3
  ncols <- 2
  nrows <- 1
  heights <- 1.5*magn
  widths  <- c(0.8,2)*magn

  pdf( paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_agg_123.pdf", sep="" ), width=sum(widths), height=sum(heights) )
  # pdf( paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_agg_123.pdf", sep="" ), width=8, height=6 )

    panel <- layout(
                  matrix( c(1:(nrows*ncols)), nrows, ncols, byrow=TRUE ),
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
    # layout.show( panel )

    ## histogram of values in bin
    lhist <- 100
    yhist <- with( dplyr::filter( nice_agg, finalcluster %in% c(1,2,3) & soilm_mean < 0.1 & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    len <- length( yhist$counts )

    par( las=1, mar=c(4,1,1,1), xaxs="i", yaxs="i" )
    plot( c(0, max( yhist$counts ) ), c( 0, 1.4 ), type="n", xlim=c(-max(yhist$counts),0), axes=FALSE, xlab="", ylab="" )
    rect( -yhist$counts, yhist$breaks[1:len], rep( 0, len ), yhist$breaks[2:(len+1)], col=rgb(0,0,0,0.15) )

    yhist_1 <- with( dplyr::filter( nice_agg, finalcluster %in% c(1) & soilm_mean < 0.1 & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_1$counts, yhist_1$breaks[1:len], rep( 0, len ), yhist_1$breaks[2:(len+1)], col=add_alpha("tomato", 0.5) )

    yhist_2 <- with( dplyr::filter( nice_agg, finalcluster %in% c(2) & soilm_mean < 0.1 & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_2$counts, yhist_2$breaks[1:len], rep( 0, len ), yhist_2$breaks[2:(len+1)], col=add_alpha("royalblue3", 0.5) )

    yhist_3 <- with( dplyr::filter( nice_agg, finalcluster %in% c(3) & soilm_mean < 0.1 & fvar > 0 & fvar < 1.4 ), hist( fvar, breaks=seq( from=0, to=1.4, length.out=lhist ), plot=FALSE ) )
    rect( -yhist_3$counts, yhist_3$breaks[1:len], rep( 0, len ), yhist_3$breaks[2:(len+1)], col=add_alpha("springgreen3", 0.5) )

    legend( "bottomleft", legend=c("cDD", "cGR", "cLS"), fill=c( add_alpha("tomato", 0.6),  add_alpha("royalblue3", 0.6),  add_alpha("springgreen3", 0.6)), bty="n", inset=c(0,0) )


    ## clusters 1 2 3 pooled
    par( las=1, mar=c(4,4,1,1), xaxs="i", yaxs="i" )
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
    text( 0.05, 1.25, "cDD, cGR, and cLS pooled", adj=c(0,0), font=2, cex=1.2 )
    abline( h=1.0, lwd=0.5 )

    legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )

  dev.off()


# ##------------------------------------------------
# ## Aligned plots for all sites combined into clusters (2 rows)
# ##------------------------------------------------
#   magn <- 2
#   ncols <- 2
#   nrows <- 3
#   heights <- c(1,1.2)*magn
#   widths  <- c(1.45,1.63)*magn

#   pdf( paste( "fig_nn_fluxnet2015/aligned_cluster/aligned_cluster.pdf", sep="" ), width=sum(widths), height=sum(heights) )

#     panel <- layout(
#                   matrix( c(1:(nrows*ncols)), nrows, ncols, byrow=FALSE ),
#                   widths=widths,
#                   heights=heights,
#                   TRUE
#                   )
#     # layout.show( panel )

#     for (iclust in 1:2){
#       ##--------------------------------------------------------
#       ## fLUE
#       ##--------------------------------------------------------
#       if (verbose) print(paste(iclust, "/", 2))
      
#       if (iclust==1) par( las=1, mar=c(0,4,1,1), xpd=FALSE, xaxs="i", yaxs="r" )
#       if (iclust==2) par( las=1, mar=c(0,4,1,3), xpd=FALSE, xaxs="i", yaxs="r" )

#       plot( c(-before,after), c(0,1.3), type="n", xlab="day after drought onset", ylab="fLUE", axes=FALSE )
#       axis( 2 )
#       if (iclust==1) axis( 4, labels=FALSE )
#       if (iclust==2) axis( 4 )
#       abline( h=1.0, col='grey40', lwd=0.5 )

#       # # lines for each event
#       # for (sitename in unique(df_dday_agg$mysitename)){
#       #   if ( siteinfo_sub$alignedcluster[ which( siteinfo_sub$mysitename==sitename ) ]==iclust && siteinfo_sub$fvar_min[ which( siteinfo_sub$mysitename==sitename ) ]<0.9 ){
#       #     sub <- dplyr::filter( df_dday_agg, mysitename==sitename )
#       #     instances <- unlist( unique( ( dplyr::select( sub, inst ) ) ) )
#       #     for (iinst in instances){
#       #       subsub <- dplyr::filter( sub, inst==iinst )
#       #       # rect( 0, -99, max( subsub$dday ), 99, col=rgb(0,0,0,0.002), border=NA )
#       #       with( subsub, lines( dday, fvar, col=add_alpha( 'tomato', 0.1 ) ) )
#       #     }
#       #   }
#       # }

#       rect( 0, -99, after, 99, col=rgb(0,0,0,0.2), border=NA )

#       ## filter data for this cluster
#       sub <- df_dday_agg %>% dplyr::filter( finalcluster==iclust )

#       ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
#       median <- sub %>% group_by( dday ) %>% 
#                         summarise( fvar = median( fvar , na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       upper  <- sub %>% group_by( dday ) %>% 
#                         summarise( fvar = quantile( fvar, 0.75, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       lower  <- sub %>% group_by( dday ) %>% 
#                         summarise( fvar = quantile( fvar, 0.25, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       uupper <- sub %>% group_by( dday ) %>% 
#                         summarise( fvar = quantile( fvar, 0.90, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       llower <- sub %>% group_by( dday ) %>% 
#                         summarise( fvar = quantile( fvar, 0.10, na.rm=TRUE ) ) %>% 
#                         complete( dday )

#       polygon( c( median$dday, rev(median$dday) ), c( lower$fvar,  rev(upper$fvar) ),  col=add_alpha("tomato", 0.3), border=NA )
#       polygon( c( median$dday, rev(median$dday) ), c( llower$fvar, rev(uupper$fvar) ), col=add_alpha("tomato", 0.3), border=NA )
#       lines( median, col='tomato2', lwd=2 )

#       text( 40, 1.2, paste("Cluster", iclust), font=2 )

#       ##--------------------------------------------------------
#       ## fAPAR
#       ##--------------------------------------------------------
#       if (iclust==1) par( las=1, mar=c(4,4,0,1), xpd=FALSE, xaxs="i", yaxs="r" )
#       if (iclust==2) par( las=1, mar=c(4,4,0,3), xpd=FALSE, xaxs="i", yaxs="r" )

#       plot( c(-before,after), c(0,1.3), type="n", ylab=expression( paste( "EVI / EVI"[0] ) ), axes=FALSE, xlab="days after drought onset" )
#       axis( 2 )
#       axis( 1, xlab="days after drought onset" )
#       if (iclust==1) axis( 4, labels=FALSE )
#       if (iclust==2) axis( 4 )
#       # axis( 1, xlab="days after drought onset" )

#       # ## lines for each event
#       # for (sitename in unique(df_dday_agg$mysitename)){
#       #   if ( siteinfo_sub$alignedcluster[ which( siteinfo_sub$mysitename==sitename ) ]==iclust ){
#       #     sub <- dplyr::filter( df_dday_agg, mysitename==sitename )
#       #     instances <- unlist( unique( ( dplyr::select( sub, inst ) ) ) )
#       #     for (iinst in instances){
#       #       subsub <- dplyr::filter( sub, inst==iinst )
#       #       with( subsub, lines( dday, evi_norm, col=add_alpha( 'springgreen3', 0.1 ) ) )
#       #     }
#       #   }
#       # }

#       rect( 0, -99, after, 99, col=rgb(0,0,0,0.2), border=NA )

#       ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
#       median <- sub %>% group_by( dday ) %>% 
#                         summarise( evi_norm = median( evi_norm , na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       upper  <- sub %>% group_by( dday ) %>% 
#                         summarise( evi_norm = quantile( evi_norm, 0.75, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       lower  <- sub %>% group_by( dday ) %>% 
#                         summarise( evi_norm = quantile( evi_norm, 0.25, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       uupper <- sub %>% group_by( dday ) %>% 
#                         summarise( evi_norm = quantile( evi_norm, 0.90, na.rm=TRUE ) ) %>% 
#                         complete( dday )
#       llower <- sub %>% group_by( dday ) %>% 
#                         summarise( evi_norm = quantile( evi_norm, 0.10, na.rm=TRUE ) ) %>% 
#                         complete( dday )

#       polygon( c( median$dday, rev(median$dday) ), c( lower$evi_norm, rev(upper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
#       polygon( c( median$dday, rev(median$dday) ), c( llower$evi_norm, rev(uupper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
#       lines( median, col='springgreen4', lwd=2 )

#     }

#   dev.off()


##------------------------------------------------
## Plot bias of NNact, NNpot, and NNvpd during dry and moist days
## aggregating data from all sites and done with all soil moisture datasets
##------------------------------------------------
  tmp <- rbind( 
                data.frame( obs=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist )$lue_obs_evi, 
                            vpd=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist )$var_nn_vpd, 
                            act=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist )$var_nn_act,
                            pot=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist )$var_nn_pot,
                            moist=rep( TRUE, nrow( dplyr::filter( nice_resh, finalcluster %in% c(1,2) & moist ) ) )
                            ),
                data.frame( obs=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist )$lue_obs_evi, 
                            vpd=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist )$var_nn_vpd, 
                            act=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist )$var_nn_act,
                            pot=dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist )$var_nn_pot,
                            moist=rep( FALSE, nrow( dplyr::filter( nice_resh, finalcluster %in% c(1,2) & !moist ) ) )
                            )
                )
  stats_vpdtest <- data.frame( 
                                rmse_nn_act_moist = Metrics::rmse( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, moist,  finalcluster %in% c(1,2) )$var_nn_act ) %>% format( digits=2 ), 
                                rmse_nn_act_dry   = Metrics::rmse( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_act ) %>% format( digits=2 ), 
                                rmse_nn_vpd_moist = Metrics::rmse( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, moist,  finalcluster %in% c(1,2) )$var_nn_vpd ) %>% format( digits=2 ), 
                                rmse_nn_vpd_dry   = Metrics::rmse( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_vpd ) %>% format( digits=2 ), 
                                rmse_nn_pot_moist = Metrics::rmse( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, moist,  finalcluster %in% c(1,2) )$var_nn_pot ) %>% format( digits=2 ), 
                                rmse_nn_pot_dry   = Metrics::rmse( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi,  dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_pot ) %>% format( digits=2 ), 

                                rsq_nn_act_moist = summary( lm( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$var_nn_act ~ dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),
                                rsq_nn_act_dry   = summary( lm( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_act ~ dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),
                                rsq_nn_vpd_moist = summary( lm( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$var_nn_vpd ~ dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),
                                rsq_nn_vpd_dry   = summary( lm( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_vpd ~ dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),
                                rsq_nn_pot_moist = summary( lm( dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$var_nn_pot ~ dplyr::filter( nice_resh,  moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),
                                rsq_nn_pot_dry   = summary( lm( dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$var_nn_pot ~ dplyr::filter( nice_resh, !moist, finalcluster %in% c(1,2) )$lue_obs_evi ) )$adj.r.squared %>% format( digits=2 ),

                                pval_nn_act_moist = with( dplyr::filter( tmp, moist ),  t.test( obs, act, paired=TRUE) )$p.value  %>% format( digits=2 ),
                                pval_nn_act_dry   = with( dplyr::filter( tmp, !moist ), t.test( obs, act, paired=TRUE) )$p.value  %>% format( digits=2 ),
                                pval_nn_vpd_moist = with( dplyr::filter( tmp, moist ),  t.test( obs, vpd, paired=TRUE) )$p.value  %>% format( digits=2 ),
                                pval_nn_vpd_dry   = with( dplyr::filter( tmp, !moist ), t.test( obs, vpd, paired=TRUE) )$p.value  %>% format( digits=2 ),
                                pval_nn_pot_moist = with( dplyr::filter( tmp, moist ),  t.test( obs, pot, paired=TRUE) )$p.value  %>% format( digits=2 ),
                                pval_nn_pot_dry   = with( dplyr::filter( tmp, !moist ), t.test( obs, pot, paired=TRUE) )$p.value  %>% format( digits=2 )
                                )

  if (makepdf) pdf( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/fig_nn_fluxnet2015/ratio_vpdtest/ratio_vpdtest_lue_obs_evi_ALL.pdf", sep="" ), width=8, height=6 )      
    
    par( las=1, xaxs="i" )

    ylim <- c(-1.2,2.6)

    ## NNpot
    bp1 <- boxplot( 
      log( var_nn_pot / lue_obs_evi ) ~ moist, 
      data =  dplyr::filter( nice_resh, finalcluster %in% c(1,2) ), 
      col   =c("grey50"), 
      ylim = ylim, 
      # xlim=c(min(soilm_threshold), max(soilm_threshold)+0.025), 
      las = 1, 
      xlab = "", 
      ylab = "log of modeled/observed", 
      yaxs = "i",
      at   = c(1,5) ,
      xlim = c(0,8),
      outline=FALSE,
      axes=FALSE
      )
    # text( 0.1, 2.7-0.5-0.03, expression(paste("R"^2)), adj=c(0,0))
    # text( 0.1, 2.5-0.5-0.03, expression(paste("RMSE")), adj=c(0,0))
    text( 0.1, 2.7-0.5-0.03, "p-value", adj=c(0,0) )

    text( 1,   2.9-0.5, expression(paste("LUE"[pot])))
    text( 1+4, 2.9-0.5, expression(paste("LUE"[pot])))

    # text( 1,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_dry   ) )
    # text( 1+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_moist ) )
    text( 1,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_dry   ) )
    text( 1+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_moist   ) )

    # text( 1,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_dry   ) )
    # text( 1+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_moist ) )

    ## NNvpd
    bp2 <- boxplot( 
      log( var_nn_vpd / lue_obs_evi ) ~ moist, 
      data  = dplyr::filter( nice_resh, finalcluster %in% c(1,2) ), 
      col   =c("grey70"), 
      xlab  ="", 
      ylab  ="", 
      yaxs  ="i",
      at   = c(1,5)+1 ,
      # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
      add   = TRUE,
      axes  = FALSE,
      outline=FALSE
      )
    text( 2,   2.9-0.5, expression(paste("LUE"[VPD])))
    text( 2+4, 2.9-0.5, expression(paste("LUE"[VPD])))

    # text( 2,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_dry   ) )
    # text( 2+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_moist ) )
    text( 2,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_dry   ) )
    text( 2+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_moist   ) )

    # text( 2,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_dry   ) )
    # text( 2+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_moist ) )

    ## NNact
    bp3 <- boxplot( 
      log( var_nn_act / lue_obs_evi ) ~ moist, 
      data  = dplyr::filter( nice_resh, finalcluster %in% c(1,2) ), 
      col  = c("springgreen"), 
      xlab  ="", 
      ylab  ="", 
      yaxs  ="i",
      at   = c(1,5)+2 ,
      # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
      add   = TRUE,
      axes  = FALSE,
      outline=FALSE
      ) 
    text( 3,   2.9-0.5, expression(paste("LUE"[act]))) 
    text( 3+4, 2.9-0.5, expression(paste("LUE"[act]))) 

    # text( 3,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_dry   ) )
    # text( 3+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_moist ) )
    text( 3,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_dry   ) )
    text( 3+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_moist   ) )

    # text( 3,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_dry   ) )
    # text( 3+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_moist ) )

    axis( 1, at=c(2,6), labels=c("dry","moist"), tick=FALSE, cex.axis=1.5 )
    axis( 2 )
    box()        
    abline( h=0.0, col='red' )
    rect( 0, -2, 4, 4, border=NA, col=rgb(0,0,0,0.2) )

  if (makepdf) dev.off()



##------------------------------------------------
## Plot bias of extreme combinations:
## 1. high soil moisture, high VPD
##------------------------------------------------
  df_test_hihi <- c()
  for (sitename in unique( nice_resh$mysitename ) ){
    tmp  <- dplyr::filter( nice_resh, mysitename==sitename )
    tmp2 <- dplyr::filter( tmp, vpd > quantile( tmp$vpd, 0.75 ) & soilm > 0.75 )
    df_test_hihi <- rbind( df_test_hihi, tmp2 )
  }
  df_test_hihi <- df_test_hihi %>% mutate( hihi=TRUE )

  df_test_lolo <- c()
  for (sitename in unique( nice_resh$mysitename ) ){
    tmp  <- dplyr::filter( nice_resh, mysitename==sitename )
    tmp2 <- dplyr::filter( tmp, vpd < quantile( tmp$vpd, 0.25 ) & soilm < 0.25 )
    df_test_lolo <- rbind( df_test_lolo, tmp2 )
  }
  df_test_lolo <- df_test_lolo %>% mutate( hihi=FALSE )

  df_test <- rbind( df_test_hihi, df_test_lolo )


  # vpd_trh <- quantile( nice_agg$vpd, 0.75 )
  # soilm_trh <- 0.75
  # df_test_hihi <- dplyr::filter( nice_agg, vpd > vpd_trh, soilm_mean > soilm_trh )
  # df_test_hihi <- df_test_hihi %>% mutate( hihi=TRUE )

  # vpd_trh <- quantile( nice_agg$vpd, 0.25 )
  # soilm_trh <- 0.25
  # df_test_lolo <- dplyr::filter( nice_agg, vpd < vpd_trh, soilm_mean < soilm_trh )
  # df_test_lolo <- df_test_lolo %>% mutate( hihi=FALSE )

  # df_test <- rbind( df_test_hihi, df_test_lolo )

  if (makepdf) pdf( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/fig_nn_fluxnet2015/ratio_vpdtest/ratio_vpdtest2_lue_obs_evi_ALL.pdf", sep="" ), width=8, height=6 )      
    
    par( las=1, xaxs="i" )

    ylim <- c(-1.2,2.6)

    ## NNpot
    bp1 <- boxplot( 
      log( var_nn_pot / lue_obs_evi ) ~ hihi, 
      data  = dplyr::filter( df_test, finalcluster %in% c(1,2) ),
      col   = c("grey50"), 
      ylim  = ylim, 
      # xlim=c(min(soilm_threshold), max(soilm_threshold)+0.025), 
      las = 1, 
      xlab = "", 
      ylab = "log of modeled/observed", 
      yaxs = "i",
      at   = c(1,5) ,
      xlim = c(0,8),
      outline=FALSE,
      axes=FALSE
      )
    # text( 0.1, 2.7-0.5-0.03, expression(paste("R"^2)), adj=c(0,0))
    # text( 0.1, 2.5-0.5-0.03, expression(paste("RMSE")), adj=c(0,0))
    # text( 0.1, 2.7-0.5-0.03, "p-value", adj=c(0,0) )

    text( 1,   2.9-0.5, expression(paste("LUE"[pot])))
    text( 1+4, 2.9-0.5, expression(paste("LUE"[pot])))

    # # text( 1,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_dry   ) )
    # # text( 1+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_moist ) )
    # text( 1,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_dry   ) )
    # text( 1+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_moist   ) )

    # # text( 1,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_dry   ) )
    # text( 1+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_moist ) )

    ## NNvpd
    bp2 <- boxplot( 
      log( var_nn_vpd / lue_obs_evi ) ~ hihi, 
      data  = dplyr::filter( df_test, finalcluster %in% c(1,2) ), 
      col   =c("grey70"), 
      xlab  ="", 
      ylab  ="", 
      yaxs  ="i",
      at   = c(1,5)+1 ,
      # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
      add   = TRUE,
      axes  = FALSE,
      outline=FALSE
      )
    text( 2,   2.9-0.5, expression(paste("LUE"[VPD])))
    text( 2+4, 2.9-0.5, expression(paste("LUE"[VPD])))

    # # text( 2,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_dry   ) )
    # # text( 2+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_moist ) )
    # text( 2,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_dry   ) )
    # text( 2+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_moist   ) )

    # # text( 2,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_dry   ) )
    # # text( 2+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_moist ) )

    ## NNact
    bp3 <- boxplot( 
      log( var_nn_act / lue_obs_evi ) ~ hihi, 
      data  = dplyr::filter( df_test, finalcluster %in% c(1,2) ), 
      col  = c("springgreen"), 
      xlab  ="", 
      ylab  ="", 
      yaxs  ="i",
      at   = c(1,5)+2 ,
      # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
      add   = TRUE,
      axes  = FALSE,
      outline=FALSE
      ) 
    text( 3,   2.9-0.5, expression(paste("LUE"[act]))) 
    text( 3+4, 2.9-0.5, expression(paste("LUE"[act]))) 

    # # text( 3,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_dry   ) )
    # # text( 3+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_moist ) )
    # text( 3,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_dry   ) )
    # text( 3+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_moist   ) )

    # # text( 3,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_dry   ) )
    # # text( 3+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_moist ) )

    axis( 1, at=c(2,6), labels=c("hihi","lolo"), tick=FALSE, cex.axis=1.5 )
    axis( 2 )
    box()        
    abline( h=0.0, col='red' )
    rect( 0, -2, 4, 4, border=NA, col=rgb(0,0,0,0.2) )

  if (makepdf) dev.off()


##------------------------------------------------
## Aligned plots for all sites combined into clusters (3 rows)
##------------------------------------------------
  nclust <- 2
  magn <- 2
  ncols <- nclust
  nrows <- 3
  heights <- c(1,1,1.2)*magn
  widths  <- c(1.45,1.63)*magn

  pdf( paste( "fig_nn_fluxnet2015/aligned_cluster/aligned_cluster_3rows.pdf", sep="" ), width=sum(widths), height=sum(heights) )

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
      print(paste(iclust, "/", nclust))
      
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


      # ## P-model bias
      # ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      # median <- sub %>% group_by( dday ) %>% 
      #                   summarise( soilm_mean = median( soilm_mean , na.rm=TRUE ) ) %>%
      #                   complete( dday )
      # upper  <- sub %>% group_by( dday ) %>% 
      #                   summarise( soilm_mean = quantile( soilm_mean, 0.75, na.rm=TRUE ) ) %>%
      #                   complete( dday )
      # uupper <- sub %>% group_by( dday ) %>% 
      #                   summarise( soilm_mean = quantile( soilm_mean, 0.90, na.rm=TRUE ) ) %>%
      #                   complete( dday )
      # lower  <- sub %>% group_by( dday ) %>% 
      #                   summarise( soilm_mean = quantile( soilm_mean, 0.25, na.rm=TRUE ) ) %>%
      #                   complete( dday )
      # llower <- sub %>% group_by( dday ) %>% 
      #                   summarise( soilm_mean = quantile( soilm_mean, 0.10, na.rm=TRUE ) ) %>%
      #                   complete( dday )

      # abline( h=1.0, col='grey40', lwd=0.5 )

      # polygon( c( median$dday, rev(median$dday) ), c( lower$soilm_mean, rev(upper$soilm_mean) ), col=add_alpha("royalblue3", 0.4), border=NA )
      # #polygon( c( median$dday, rev(median$dday) ), c( llower$bias_pmodel, rev(uupper$bias_pmodel) ), col=add_alpha("royalblue3", 0.4), border=NA )
      # lines( median, col='royalblue4', lwd=2 )

      # ## MODIS BIAS
      # ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      # median <- df_dday_modis_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                                 group_by( dday ) %>% 
      #                                 summarise( bias_modis = median( bias_modis , na.rm=TRUE ) )
      # upper  <- df_dday_modis_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                                 group_by( dday ) %>% 
      #                                 summarise( bias_modis = quantile( bias_modis, 0.75, na.rm=TRUE ) )
      # uupper <- df_dday_modis_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                                 group_by( dday ) %>% 
      #                                 summarise( bias_modis = quantile( bias_modis, 0.90, na.rm=TRUE ) )
      # lower  <- df_dday_modis_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                                 group_by( dday ) %>% 
      #                                 summarise( bias_modis = quantile( bias_modis, 0.25, na.rm=TRUE ) )
      # llower <- df_dday_modis_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                                 group_by( dday ) %>% 
      #                                 summarise( bias_modis = quantile( bias_modis, 0.10, na.rm=TRUE ) )


      # polygon( c( median$dday, rev(median$dday) ), c( lower$bias_modis, rev(upper$bias_modis) ), col=add_alpha("orchid3", 0.4), border=NA )
      # #polygon( c( median$dday, rev(median$dday) ), c( llower$bias_modis, rev(uupper$bias_modis) ), col=add_alpha("orchid3", 0.4), border=NA )
      # lines( median, col='orchid4', lwd=2 )


      # ## MTE BIAS
      # ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      # median <- df_dday_mte_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                               group_by( dday ) %>% 
      #                               summarise( bias_mte = median( bias_mte , na.rm=TRUE ) )
      # upper  <- df_dday_mte_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                               group_by( dday ) %>% 
      #                               summarise( bias_mte = quantile( bias_mte, 0.75, na.rm=TRUE ) )
      # uupper <- df_dday_mte_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                               group_by( dday ) %>% 
      #                               summarise( bias_mte = quantile( bias_mte, 0.90, na.rm=TRUE ) )
      # lower  <- df_dday_mte_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                               group_by( dday ) %>% 
      #                               summarise( bias_mte = quantile( bias_mte, 0.25, na.rm=TRUE ) )
      # llower <- df_dday_mte_agg %>% dplyr::filter( alignedcluster==iclust ) %>% 
      #                               group_by( dday ) %>% 
      #                               summarise( bias_mte = quantile( bias_mte, 0.10, na.rm=TRUE ) )


      # polygon( c( median$dday, rev(median$dday) ), c( lower$bias_mte, rev(upper$bias_mte) ), col=add_alpha("goldenrod3", 0.4), border=NA )
      # #polygon( c( median$dday, rev(median$dday) ), c( llower$bias_mte, rev(uupper$bias_mte) ), col=add_alpha("goldenrod3", 0.4), border=NA )
      # lines( median, col='goldenrod4', lwd=2 )


    }

  dev.off()




# ##------------------------------------------------
# ## BOXPLOT GPP BIAS BY MODEL
# ##------------------------------------------------
#   print("plot aggregated gpp bias...")

#   pdf( paste( "fig_nn_fluxnet2015/boxplot_gppbias/boxplot_gppbias_ALL.pdf", sep="" ), width=7, height=6 )

#     ylim <- c(0,5)
#     xlim <- c(0,6)

#     plot( xlim, ylim, axes=FALSE, type='n', xlab="", ylab="GPP mod. / GPP obs.", main="ALL" )
#     rect( c(0,4), rep(-99,2), c(2,6), rep(99,2), col=rgb(0,0,0,0.2), border=NA )

#     par( las=1, mar=c(4.5,4.5,2,2) )
#     bp1 <- boxplot( 
#       bias_pmodel ~ is_drought_byvar, 
#       data = filter( nice_agg, quadfitcluster==1 | quadfitcluster==2 ),
#       col  = c( "royalblue3", "tomato" ), 
#       las  = 1, 
#       xlab = "",
#       ylab = "",
#       main = "",
#       at   = c(5,6)-0.5,
#       yaxs ="i",
#       add  = TRUE,
#       axes = FALSE,
#       outline = FALSE
#       # outpch=16,
#       # outcol=rgb(0,0,0,0.3)
#       # boxwex = 0.2,
#       # outline=FALSE
#       )
#     axis( 2 )
#     linmod <- lm( bias_pmodel ~ is_drought_byvar, data = filter( nice_agg, quadfitcluster==1 | quadfitcluster==2 ) )
#     pval   <- lmp(linmod)
#     text( 4.5, ylim[2]*0.95, "P-model", adj=c(0,0), cex=1.0 )
#     #text( 4.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( filter( nice_agg, quadfitcluster==1 | quadfitcluster==2 ) ) ) ), adj=c(0,0), cex=0.75 )

#     par( las=1, mar=c(4.5,4.5,2,2) )
#     bp1 <- boxplot( 
#       bias_modis ~ is_drought_byvar, 
#       data = filter( nice_to_modis_agg, quadfitcluster==1 | quadfitcluster==2 ),
#       col  = c( "royalblue3", "tomato" ), 
#       las  = 1, 
#       xlab = "",
#       ylab = "",
#       at   = c(3,4)-0.5,
#       yaxs ="i",
#       add  = TRUE,
#       # axes = FALSE,
#       xlab = "",
#       ylab = "",
#       axes = FALSE,
#       outline = FALSE
#       # outpch=16,
#       # outcol=rgb(0,0,0,0.3)
#       # boxwex = 0.2,
#       # outline=FALSE
#       )
#     linmod <- lm( bias_modis ~ is_drought_byvar, data = filter( nice_to_modis_agg, quadfitcluster==1 | quadfitcluster==2 ) )
#     pval   <- lmp(linmod)
#     text( 2.5, ylim[2]*0.95, "MODIS", adj=c(0,0), cex=1.0 )
#     #text( 2.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( gpp_agg_modis ) ) ), adj=c(0,0), cex=0.75 )

#     par( las=1, mar=c(4.5,4.5,2,2) )
#     if (nrow(nice_to_mte_agg)>0){
#       bp1 <- boxplot( 
#         bias_mte ~ is_drought_byvar, 
#         data = filter( nice_to_mte_agg, quadfitcluster==1 | quadfitcluster==2 ),
#         col  = c( "royalblue3", "tomato" ), 
#         las  = 1, 
#         xlab = "",
#         ylab = "",
#         at   = c(1,2)-0.5,
#         yaxs ="i",
#         add  = TRUE,
#         # axes = FALSE,
#         xlab = "",
#         ylab = "",
#         axes = FALSE,
#         outline = FALSE
#         # outpch=16,
#         # outcol=rgb(0,0,0,0.3)
#         # boxwex = 0.2,
#         # outline=FALSE
#         )
#       linmod <- lm( bias_mte ~ is_drought_byvar, data = filter( nice_to_mte_agg, quadfitcluster==1 | quadfitcluster==2 ) )
#       pval   <- lmp(linmod)
#       text( 0.5, ylim[2]*0.95, "MTE", adj=c(0,0), cex=1.0 )
#       #text( 0.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( gpp_agg_mte ) ) ), adj=c(0,0), cex=0.75 )
#     }

#     abline( h=1, lwd=0.5 )

#     legend( "bottomright", c("non-drought", "drought"), bty="n", fill=c("royalblue3", "tomato"))

#   dev.off()
#   print("done.")



# ##------------------------------------------------
# ## PLOT AGGREGATED GPP BIAS BY CLUSTER
# ##------------------------------------------------
#   print("plot aggregated gpp bias by cluster...")
#   for (icluster in 1:nclusters){

#     ## filter by cluster
#     gpp_agg_sub       <- dplyr::filter( gpp_agg, quadfitcluster==icluster )
#     gpp_agg_sub_modis <- dplyr::filter( gpp_agg_modis, quadfitcluster==icluster )
#     gpp_agg_sub_mte   <- dplyr::filter( gpp_agg_mte, quadfitcluster==icluster )


#     pdf( paste( "fig_nn_fluxnet2015/boxplot_gppbias/boxplot_gppbias_CLUSTER", icluster,".pdf", sep="" ), width=7, height=6 )

#       ylim <- c(-2.5,5)
#       xlim <- c(0,6)

#       plot( xlim, ylim, axes=FALSE, type='n', xlab="", ylab="log of relative bias", main=paste("cluster", icluster) )
#       rect( c(0,4), rep(-99,2), c(2,6), rep(99,2), col=rgb(0,0,0,0.2), border=NA )

#       par( las=1, mar=c(4.5,4.5,2,2) )
#       bp1 <- boxplot( 
#         bias_pmodel ~ is_drought_byvar, 
#         data = gpp_agg_sub,
#         col  = c( "royalblue3", "tomato" ), 
#         las  = 1, 
#         xlab = "",
#         ylab = "",
#         main = "",
#         at   = c(5,6)-0.5,
#         yaxs ="i",
#         add  = TRUE,
#         axes = FALSE,
#         outline = FALSE
#         # outpch=16,
#         # outcol=rgb(0,0,0,0.3)
#         # boxwex = 0.2,
#         # outline=FALSE
#         )
#       axis( 2 )
#       linmod <- lm( bias_pmodel ~ is_drought_byvar, data = gpp_agg_sub )
#       pval   <- lmp(linmod)
#       text( 4.5, ylim[2]*0.95, "P-model", adj=c(0,0), cex=1.0 )
#       text( 4.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( gpp_agg_sub ) ) ), adj=c(0,0), cex=0.75 )

#       par( las=1, mar=c(4.5,4.5,2,2) )
#       if (nrow(aggregate(bias_modis ~ is_drought_byvar, data = gpp_agg_sub_modis, FUN=median, na.rm=T ))==2){
#         bp1 <- boxplot( 
#           bias_modis ~ is_drought_byvar, 
#           data = gpp_agg_sub_modis,
#           col  = c( "royalblue3", "tomato" ), 
#           las  = 1, 
#           xlab = "",
#           ylab = "",
#           at   = c(3,4)-0.5,
#           yaxs ="i",
#           add  = TRUE,
#           # axes = FALSE,
#           xlab = "",
#           ylab = "",
#           axes = FALSE,
#           outline = FALSE
#           # outpch=16,
#           # outcol=rgb(0,0,0,0.3)
#           # boxwex = 0.2,
#           # outline=FALSE
#           )
#         linmod <- lm( bias_modis ~ is_drought_byvar, data = gpp_agg_sub_modis )
#         pval   <- lmp(linmod)
#         text( 2.5, ylim[2]*0.95, "MODIS", adj=c(0,0), cex=1.0 )
#         text( 2.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( gpp_agg_sub_modis ) ) ), adj=c(0,0), cex=0.75 )
#       } 

#       par( las=1, mar=c(4.5,4.5,2,2) )
#       if (nrow(gpp_agg_sub_mte)>0){
#         if (nrow(aggregate(bias_mte ~ is_drought_byvar, data = gpp_agg_sub_mte, FUN=median, na.rm=T ))==2){
#           bp1 <- boxplot( 
#             bias_mte ~ is_drought_byvar, 
#             data = gpp_agg_sub_mte,
#             col  = c( "royalblue3", "tomato" ), 
#             las  = 1, 
#             xlab = "",
#             ylab = "",
#             at   = c(1,2)-0.5,
#             yaxs ="i",
#             add  = TRUE,
#             # axes = FALSE,
#             xlab = "",
#             ylab = "",
#             axes = FALSE,
#             outline = FALSE
#             # outpch=16,
#             # outcol=rgb(0,0,0,0.3)
#             # boxwex = 0.2,
#             # outline=FALSE
#             )
#           linmod <- lm( bias_mte ~ is_drought_byvar, data = gpp_agg_sub_mte )
#           pval   <- lmp(linmod)
#           text( 0.5, ylim[2]*0.95, "MTE", adj=c(0,0), cex=1.0 )
#           text( 0.1, ylim[2]*0.85, paste( "P(>F) =", format( pval, digits = 2 ), ", N =", as.character( nrow( gpp_agg_sub_mte ) ) ), adj=c(0,0), cex=0.75 )
#         }
#       }

#       abline( h=0, lwd=0.5 )

#       legend( "bottomright", c("non-drought", "drought"), bty="n", fill=c("royalblue3", "tomato"))

#     dev.off()
#     print("done.")

#   }

