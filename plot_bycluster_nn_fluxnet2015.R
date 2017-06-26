library(dplyr)
library(LSD)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )

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
# nice_to_modis_agg <- nice_to_modis_agg %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster  ) )
# nice_to_mte_agg   <- nice_to_mte_agg   %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster  ) )

## Load aggregated data from aligned MODIS and MTE data
load( "data/data_aligned_agg.Rdata" )
df_dday_agg       <- df_dday_agg       %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster ) )
# df_dday_modis_agg <- df_dday_modis_agg %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster  ) )
# df_dday_mte_agg   <- df_dday_mte_agg   %>% left_join( dplyr::select( overview, mysitename, alignedcluster, quadfitcluster, finalcluster  ) )

verbose <- FALSE

before <- -30
after <- 100

##------------------------------------------------
## fLUE vs. soilmoisture
##------------------------------------------------
nice_agg$ratio_obs_mod <- remove_outliers( nice_agg$ratio_obs_mod, coef=5 )

pdf( paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_agg_ALL.pdf", sep="" ), width=8, height=7 )

  ## cluster 1: big change in fAPAR
  par( las=1, mfrow=c(2,2), mar=c(4,4,3,1) )
  with( 
        dplyr::filter( nice_agg, alignedcluster==1 & quadfitcluster!=2 & finalcluster!=4 ),
        heatscatter( 
                    soilm_mean, 
                    fvar,
                    ylim=c(0,1.3),
                    main="",
                    xlab="soil water content (fraction)",
                    ylab="fLUE"
                    ) 

      )
  mtext( "Cluster 1", line=1, adj=c(0.5) )
  abline( h=1.0, lwd=0.5 )

  ## cluster 2: small change in fAPAR
  with( 
        dplyr::filter( nice_agg, alignedcluster==2 & quadfitcluster!=2 & finalcluster!=4 ),
        heatscatter( 
                    soilm_mean, 
                    fvar,
                    ylim=c(0,1.3),
                    xlim=c(0,1),
                    main="",
                    xlab="soil water content (fraction)",
                    ylab="fLUE"
                    ) 

      )
  mtext( "Cluster 2", line=1, adj=c(0.5) )
  abline( h=1.0, lwd=0.5 )

  ## cluster 3: low LUE reduction with soil moisture (groundwater access?)
  with( 
        dplyr::filter( nice_agg, quadfitcluster==2 ),
        heatscatter( 
                    soilm_mean, 
                    fvar,
                    ylim=c(0,1.3),
                    main="",
                    xlab="soil water content (fraction)",
                    ylab="fLUE"
                    ) 

      )
  mtext( "Cluster 3", line=1, adj=c(0.5) )
  abline( h=1.0, lwd=0.5 )

  ## cluster 4: no data in low soil moisture bin
  with( 
        dplyr::filter( nice_agg, finalcluster==4 ),
        heatscatter( 
                    soilm_mean, 
                    fvar,
                    ylim=c(0,1.3),
                    xlim=c(0,1),
                    main="",
                    xlab="soil water content (fraction)",
                    ylab="fLUE"
                    ) 

      )
  mtext( "Cluster 4", line=1, adj=c(0.5) )
  abline( h=1.0, lwd=0.5 )

  legend( "bottomleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )

dev.off()


##------------------------------------------------
## Aligned plots for all sites combined into clusters (2 rows)
##------------------------------------------------
  nclust <- 2
  magn <- 2
  ncols <- nclust
  nrows <- 2
  heights <- c(1,1.2)*magn
  widths  <- c(1.45,1.63)*magn

  pdf( paste( "fig_nn_fluxnet2015/aligned_cluster/aligned_cluster.pdf", sep="" ), width=sum(widths), height=sum(heights) )

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
      if (verbose) print(paste(iclust, "/", nclust))
      
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
      if (iclust==1 | iclust==2){
        sub <- df_dday_agg %>% dplyr::filter( alignedcluster==iclust & quadfitcluster!=2 & finalcluster!=4 )
      } else if (iclust==3) {
        sub <- df_dday_agg %>% dplyr::filter( quadfitcluster==2 )
      } 

      ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
      median <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = median( fvar , na.rm=TRUE ) )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.75, na.rm=TRUE ) )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.25, na.rm=TRUE ) )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.90, na.rm=TRUE ) )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( fvar = quantile( fvar, 0.10, na.rm=TRUE ) )

      polygon( c( median$dday, rev(median$dday) ), c( lower$fvar, rev(upper$fvar) ), col=add_alpha("tomato", 0.3), border=NA )
      polygon( c( median$dday, rev(median$dday) ), c( llower$fvar, rev(uupper$fvar) ), col=add_alpha("tomato", 0.3), border=NA )
      lines( median, col='tomato2', lwd=2 )

      text( 40, 1.2, paste("Cluster", iclust), font=2 )

      ##--------------------------------------------------------
      ## fAPAR
      ##--------------------------------------------------------
      if (iclust==1) par( las=1, mar=c(4,4,0,1), xpd=FALSE, xaxs="i", yaxs="r" )
      if (iclust==2) par( las=1, mar=c(4,4,0,3), xpd=FALSE, xaxs="i", yaxs="r" )

      plot( c(-before,after), c(0,1.3), type="n", ylab=expression( paste( "EVI / EVI"[0] ) ), axes=FALSE, xlab="days after drought onset" )
      axis( 2 )
      axis( 1, xlab="days after drought onset" )
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
                        summarise( evi_norm = median( evi_norm , na.rm=TRUE ) )
      upper  <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.75, na.rm=TRUE ) )
      lower  <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.25, na.rm=TRUE ) )
      uupper <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.90, na.rm=TRUE ) )
      llower <- sub %>% group_by( dday ) %>% 
                        summarise( evi_norm = quantile( evi_norm, 0.10, na.rm=TRUE ) )

      polygon( c( median$dday, rev(median$dday) ), c( lower$evi_norm, rev(upper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
      polygon( c( median$dday, rev(median$dday) ), c( llower$evi_norm, rev(uupper$evi_norm) ), col=add_alpha("springgreen3", 0.3), border=NA )
      lines( median, col='springgreen4', lwd=2 )

    }

  dev.off()
