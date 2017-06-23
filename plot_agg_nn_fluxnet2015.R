rm(list=ls(all=TRUE))

library(dplyr)

source( paste( myhome, "/analyse_modobs.R", sep="" ) )
source( paste( myhome, "/remove_outliers.R", sep="" ) )

## Manual settings ----------------
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
nrep       = 5
dotrain    = FALSE
overwrite_modis = TRUE
overwrite_mte = TRUE
##--------------------------

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

load( paste("data/missing_pri_", nam_target, char_fapar, ".Rdata", sep="") )
siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( paste( "data/nice_agg_", char_wgt, nam_target, ".Rdata", sep="" ) )       # loads 'nice_agg'
load( paste( "data/nice_modis_agg_", char_wgt, nam_target, ".Rdata", sep="" ) ) # loads 'nice_to_modis_agg'
load( paste( "data/nice_mte_agg_", char_wgt, nam_target, ".Rdata", sep="" ) )   # loads 'nice_to_mte_agg'

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_modis_agg', 'df_dday_mte_agg', 

## add vegetation type info to nice_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )


##--------------------------------------
## MOD VS OBS OF ALL DATA AGGREGATED
##--------------------------------------

  ##--------------------------------------
  ## NN evaluation
  ##--------------------------------------
    ## GPP
    pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_gpp_rct_ALL_FROMNICE", char_wgt, ".pdf", sep=""), width=8, height=8 )
      par( mfrow=c(2,2) )
      stats_tmp <- analyse_modobs( 
                                  nice_agg$gpp_nn_act, 
                                  nice_agg$gpp_obs, 
                                  plot.title=expression( paste("NN"[pot])), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65      
                                  )
      mtext( "a", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_obs, 
                                  plot.title=expression( paste("NN"[pot], "   moist days")), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "b", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, is_drought_byvar )$gpp_obs, 
                                  plot.title=expression( paste("NN"[pot], "   dry days")), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "c", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_act, 
                                  plot.title=expression( paste("NN"[pot], " vs. NN"[act], "  moist days")), 
                                  plot.xlab=expression(paste("actual GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("potential GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "d", line = 1, font = 2, adj=c(0,0) )

    dev.off()

    ## LUE
    pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_lue_rct_ALL_FROMNICE", char_wgt, ".pdf", sep=""), width=8, height=8 )
      par( mfrow=c(2,2) )
      stats_tmp <- analyse_modobs( 
                                  nice_agg$var_nn_act, 
                                  nice_agg$lue_obs_evi, 
                                  plot.title=expression( paste("NN"[pot])), 
                                  plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "a", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$var_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$lue_obs_evi, 
                                  plot.title=expression( paste("NN"[pot], "   moist days")), 
                                  plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "b", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, is_drought_byvar )$var_nn_pot, 
                                  dplyr::filter( nice_agg, is_drought_byvar )$lue_obs_evi, 
                                  plot.title=expression( paste("NN"[pot], "   dry days")), 
                                  plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "c", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$var_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$var_nn_act, 
                                  plot.title=expression( paste("NN"[pot], " vs. NN"[act], "   moist days")), 
                                  plot.xlab=expression(paste("actual LUE (gC mol"^{-1}, ")")),
                                  plot.ylab=expression(paste("potential LUE (gC mol"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "d", line = 1, font = 2, adj=c(0,0) )

    dev.off()


##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------
  df_dday_agg$ratio_obs_mod             <- remove_outliers( df_dday_agg$ratio_obs_mod, coef=5 )
  df_dday_modis_agg$ratio_obs_mod_modis <- remove_outliers( df_dday_modis_agg$ratio_obs_mod_modis, coef=5 )
  df_dday_mte_agg$ratio_obs_mod_mte     <- remove_outliers( df_dday_mte_agg$ratio_obs_mod_mte, coef=5 )

  magn <- 3
  ncols <- 4
  nrows <- 2
  widths <- c(magn, 0.3*magn, magn, 0.3*magn )
  heights <- 1.2*c(magn,0.8*magn)
  order <- matrix(c(1,1,2,2,3,4,5,6),nrows,ncols,byrow=TRUE)

  pdf( "fig_nn_fluxnet2015/bias_vs_fvar/bias_vs_fvar_ALL.pdf", width=sum(widths), height=sum(heights) )

    panel <- layout(
                    order,
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show(panel)

    #---------------------------------------------------------
    # P-model
    #---------------------------------------------------------
    ## point cloud
    par( las=1, mar=c(5,4.5,3,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,2.5)
    with( 
          dplyr::filter( df_dday_agg, ratio_obs_mod<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod, 
                      xlab="fLUE",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main=""
                    )
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "P-model", line=1, adj=0.5 )

    ## use only data during droughts for stats
    sub <- dplyr::filter( df_dday_agg, is_drought_byvar==1 ) 
    stats <- analyse_modobs( sub$ratio_obs_mod, sub$fvar, do.plot=FALSE )

    # write stats into plot
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )
    
    ## draw the legend
    legend( "topleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", inset=c(0.05,0.25), cex=0.8 )


    # Distribution 
    par( las=1, mar=c(5,0,3,5), xpd=FALSE )

    boxplot( filter( nice_agg, !is_drought_byvar )$ratio_obs_mod, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50', xlim=c(0,5), at=0.5 )
    lines( c(-2,1), c(1,1), lwd=0.5, lty=2 )
    # abline( h=1.0, lwd=0.5, lty=2 )


    #---------------------------------------------------------
    # MODIS
    #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,4)
    with( 
          dplyr::filter( df_dday_modis_agg, ratio_obs_mod_modis<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_modis, 
                      xlab="fLUE",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      cexplot=1.2,
                      main=""
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "MODIS", line=1, adj=0.5 )

    ## use only data during droughts for stats
    sub <- dplyr::filter( df_dday_modis_agg, is_drought_byvar==1 ) 
    stats <- analyse_modobs( sub$ratio_obs_mod_modis, sub$fvar, do.plot=FALSE )

    ## write stats into plot
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    # text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( 1.0 / filter( nice_to_modis_agg, !is_drought_byvar )$bias_modis, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )


    #---------------------------------------------------------
    # MTE
    #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    with( 
          dplyr::filter( df_dday_mte_agg, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with heatscatter()
          plot( 
                fvar, 
                ratio_obs_mod_mte, 
                xlab="fLUE",
                ylab="GPP observed / GPP modelled",
                xlim=c(0,1.2),
                ylim=ylim,
                cex=1.2,
                pch=16,
                col=add_alpha("black", 0.3),
                main=""
              ) 

        )
    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "FLUXCOM MTE", line=1, adj=0.5 )


    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( 1.0 / filter( nice_to_mte_agg, !is_drought_byvar )$bias_mte, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )

  dev.off()



##------------------------------------------------
## CCI and PRI
##------------------------------------------------
  df_dday_agg$dscci <- remove_outliers( df_dday_agg$dscci, coef=5 )
  df_dday_agg$dspri <- remove_outliers( df_dday_agg$dspri, coef=5 )

  ##------------------------------------------------
  ## EVI, CCI and PRI vs. GPP
  ##------------------------------------------------
  pdf( "fig_nn_fluxnet2015/cci_vs_gpp/evi_vs_gpp_ALL.pdf" )
    xlim <- c(0,25)
    ylim <- c(0,1)
    # sub <- dplyr::filter( nice_agg, classid %in% c("GRA") )
    # sub <- dplyr::filter( nice_agg, mysitename=="FI-Hyy" )
    sub <- nice_agg
    par(las=1)
    heatscatter( sub$gpp_obs, sub$evi, ylim=ylim, xlim=xlim, ylab="EVI", xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), main="" )
    linmod <- lm( evi ~ gpp_obs, data=sub )
    abline( linmod, col="red", lty=2 )
    rsq <- summary( linmod )$adj.r.squared
    numb <- sum(!is.na(sub$gpp_obs) & !is.na(sub$cci))
    text( 20, 0.2, bquote( italic(R)^2 == .(format( rsq, digits = 2) ) ),  adj=0.0, cex=1.0 )
    text( 20, 0.25, paste( "N =", format( numb, digits = 1 ) ), adj=0.0, cex=1.0 )
    legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=16, col=colorRampPalette( c("gray80", "navy", "red", "yellow"))(5), bty="n",  cex=0.8 )
  dev.off()

  if ("cci" %in% names(nice_agg)){
    pdf( "fig_nn_fluxnet2015/cci_vs_gpp/cci_vs_gpp_ALL.pdf" )
      xlim <- c(0,25)
      ylim <- c(-0.4, 0.4)
      # sub <- dplyr::filter( nice_agg, classid %in% c("GRA") )
      # sub <- dplyr::filter( nice_agg, mysitename=="FI-Hyy" )
      sub <- nice_agg
      par(las=1)
      heatscatter( sub$gpp_obs, sub$cci, ylim=ylim, xlim=xlim, ylab="CCI", xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), main="" )
      linmod <- lm( cci ~ gpp_obs, data=sub )
      abline( linmod, col="red", lty=2 )
      rsq <- summary( linmod )$adj.r.squared
      numb <- sum(!is.na(sub$gpp_obs) & !is.na(sub$cci))
      text( 20, -0.2, bquote( italic(R)^2 == .(format( rsq, digits = 2) ) ),  adj=0.0, cex=1.0 )
      text( 20, -0.16, paste( "N =", format( numb, digits = 1 ) ), adj=0.0, cex=1.0 )
      legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=16, col=colorRampPalette( c("gray80", "navy", "red", "yellow"))(5), bty="n",  cex=0.8 )
    dev.off()
  }

  if ("pri" %in% names(nice_agg)){
    pdf( "fig_nn_fluxnet2015/cci_vs_gpp/pri_vs_gpp_ALL.pdf" )
      xlim <- c(0,25)
      ylim <- c(-0.4, 0.1)
      sub <- nice_agg
      par(las=1)
      heatscatter( sub$gpp_obs, sub$pri, ylim=ylim, xlim=xlim, ylab="PRI", xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), main="" )
      linmod <- lm( cci ~ gpp_obs, data=sub )
      abline( linmod, col="red", lty=2 )
      rsq <- summary( linmod )$adj.r.squared
      numb <- sum(!is.na(sub$gpp_obs) & !is.na(sub$cci))
      text( 20, -0.2, bquote( italic(R)^2 == .(format( rsq, digits = 2) ) ),  adj=0.0, cex=1.0 )
      text( 20, -0.16, paste( "N =", format( numb, digits = 1 ) ), adj=0.0, cex=1.0 )
      legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=16, col=colorRampPalette( c("gray80", "navy", "red", "yellow"))(5), bty="n",  cex=0.8 )
    dev.off()
  }

  # ##------------------------------------------------
  # ## CCI and PRI vs. LUE - NO CORRELATION WITH R2 OF 0.024
  # ##------------------------------------------------
  # if ("cci" %in% names(nice_agg)){
  #   xlim <- c(0,3.5)
  #   ylim <- c(-0.4, 0.4)
  #   # sub <- dplyr::filter( nice_agg, classid %in% c("GRA") )
  #   sub <- nice_agg

  #   pdf( "fig_nn_fluxnet2015/cci_vs_lue/cci_vs_lue_ALL.pdf" )
  #     par(las=1)
  #     heatscatter( sub$lue_obs_evi, sub$cci, ylab="CCI", xlab="LUE", main="", ylim=ylim, xlim=xlim )
  #     linmod <- lm( cci ~ lue_obs_evi, data=sub )
  #     abline( linmod, col="red", lty=2 )
  #     rsq <- summary( linmod )$adj.r.squared
  #     numb <- sum(!is.na(sub$lue_obs_evi) & !is.na(sub$cci))
  #     text( 3.0, -0.2, bquote( italic(R)^2 == .(format( rsq, digits = 2) ) ),  adj=0.0, cex=1.0 )
  #     text( 3.0, -0.16, paste( "N =", format( numb, digits = 1 ) ), adj=0.0, cex=1.0 )
  #   dev.off()
  # }


  ##------------------------------------------------
  ## CCI and PRI
  ##------------------------------------------------
  if ("dscci" %in% names(df_dday_agg) && "dspri" %in% names(df_dday_agg)){

    pdf("fig_nn_fluxnet2015/cci_vs_fvar/pri_vs_fvar_ALL.pdf")

      par(las=1)
      xlim <- c(0,1.2)
      ylim <- c(0,1.5)

      with( 
            df_dday_agg,
            plot( 
                  fvar, 
                  dscci, 
                  xlab="fLUE",
                  ylab="scaled and normalised CCI (PRI)",
                  xlim=xlim,
                  ylim=ylim,
                  pch=16,
                  col=add_alpha('goldenrod4', 0.3),
                  cex=0.7,
                  main="All"
                )
          )

     with( 
            df_dday_agg,
            points( 
                    fvar,
                    dspri, 
                    pch=16,
                    col=add_alpha('cadetblue4', 0.3),
                    cex=0.7
                  )
          )

     with( 
            df_dday_aggbydday_agg,
            points( 
                    fvar_med, 
                    dscci_med, 
                    pch=16,
                    col=add_alpha('goldenrod4', 0.7),
                    cex=1.3
                  )
          )

     with( 
            df_dday_aggbydday_agg,
            points( 
                    fvar_med, 
                    dspri_med, 
                    pch=16,
                    col=add_alpha('cadetblue4', 0.7),
                    cex=1.3
                  )
          )

      abline( h=1.0, lwd=0.5, lty=2 )
      abline( v=1.0, lwd=0.5, lty=2 )
      lines( c(-99,99), c(-99,99), col='red' )

      ## linear fit
      ## CCI
      linmod <- lm( dscci ~ fvar, data=df_dday_agg )
      abline( linmod, col="goldenrod4", lty=1 )

      ## PRI
      linmod <- lm( dspri ~ fvar, data=df_dday_agg )
      abline( linmod, col="cadetblue4", lty=1 )

      ## use only data during droughts for stats
      sub <- dplyr::filter( df_dday_agg, is_drought_byvar==1 ) 
      stats_cci <- analyse_modobs( sub$dscci, sub$fvar, do.plot=FALSE )
      stats_pri <- analyse_modobs( sub$dspri, sub$fvar, do.plot=FALSE )

      # write stats into plot
      x0 <- 0.03*xlim[2]
      y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
      # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
      text( x0, y0+0.15, "CCI:",  adj=0.0, cex=1.2, col='goldenrod4' )
      text( x0, y0+0.05, "PRI:",  adj=0.0, cex=1.2, col='cadetblue4' )

      text( x0+0.1, y0+0.15, bquote( italic(R)^2 == .(format( stats_cci$rsq, digits = 2) ) ),  adj=0.0, cex=1.2, col='goldenrod4' )
      text( x0+0.1, y0+0.05, bquote( italic(R)^2 == .(format( stats_pri$rsq, digits = 2) ) ),  adj=0.0, cex=1.2, col='cadetblue4' )
      # text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )

      # mtext( paste( "N =", format( stats$N, digits = 1 ) ), side=3, line=0, adj=0.0, cex=0.8 )
      # mtext( bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ), side=3, line=1, adj=0.0, cex=0.8 )

    dev.off()

  }

  if ("dscci" %in% names(df_dday_agg) && "dspri" %in% names(df_dday_agg)){
    for (sitename in unique(df_dday_agg$mysitename)){

      if ( !is.element( sitename, missing_pri ) ) {

        pdf( paste( "fig_nn_fluxnet2015/cci_vs_fvar/pri_vs_fvar_", sitename, ".pdf", sep="" ))

          par(las=1)
          xlim <- c(0,1.2)
          ylim <- c(0,1.5)
          with( 
                filter( df_dday_agg, mysitename==sitename ),
                plot( 
                      fvar, 
                      dscci, 
                      xlab="fLUE",
                      ylab="scaled and normalised CCI (PRI)",
                      xlim=xlim,
                      ylim=ylim,
                      pch=16,
                      col=add_alpha('goldenrod4', 0.3),
                      cex=0.8,
                      main=paste( sitename, " ", dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid )
                    )
              )
          with( 
                filter( df_dday_aggbydday_agg, mysitename==sitename ),
                points( 
                        fvar_med, 
                        dscci_med,
                        pch=16,
                        col=add_alpha('goldenrod4', 0.7),
                        cex=1.3
                      )
              )

          with( 
                filter( df_dday_agg, mysitename==sitename ),
                points( 
                        fvar, 
                        dspri, 
                        pch=16,
                        col=add_alpha('cadetblue4', 0.3),
                        cex=0.8,
                        main=paste( sitename, " ", dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid )
                      )
              )
          with( 
                filter( df_dday_aggbydday_agg, mysitename==sitename ),
                points( 
                        fvar_med, 
                        dspri_med,
                        pch=16,
                        col=add_alpha('cadetblue4', 0.7),
                        cex=1.3
                      )
              )

          abline( h=1.0, lwd=0.5, lty=2 )
          abline( v=1.0, lwd=0.5, lty=2 )
          lines( c(-99,99), c(-99,99), col='red', lty=2 )

          ## linear fit
          ## CCI
          linmod <- lm( dscci ~ fvar, data=filter( df_dday_agg,  mysitename==sitename ) )
          abline( linmod, col="goldenrod4", lty=1 )

          ## PRI
          linmod <- lm( dspri ~ fvar, data=filter( df_dday_agg,  mysitename==sitename ) )
          abline( linmod, col="cadetblue4", lty=1 )

          ## use only data during droughts for stats
          sub <- dplyr::filter( df_dday_agg, is_drought_byvar==1 & mysitename==sitename ) 
          stats_cci <- analyse_modobs( sub$dscci, sub$fvar, do.plot=FALSE )
          stats_pri <- analyse_modobs( sub$dspri, sub$fvar, do.plot=FALSE )

          # write stats into plot
          x0 <- 0.03*xlim[2]
          y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
          # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
          text( x0, y0+0.15, "PRI:",  adj=0.0, cex=1.2, col='cadetblue4' )
          text( x0, y0+0.05, "CCI:",  adj=0.0, cex=1.2, col='goldenrod4' )

          text( x0+0.1, y0+0.15, bquote( italic(R)^2 == .(format( stats_pri$rsq, digits = 2) ) ),  adj=0.0, cex=1.2, col='cadetblue4' )
          text( x0+0.1, y0+0.05, bquote( italic(R)^2 == .(format( stats_cci$rsq, digits = 2) ) ),  adj=0.0, cex=1.2, col='goldenrod4' )

        dev.off()

      }
    }
  }

##------------------------------------------------
## EVI vs. FPAR
##------------------------------------------------
pdf("fig_nn_fluxnet2015/fpar_vs_evi.pdf")
  par(las=1)
  with( nice_agg,
        heatscatter( evi, fpar, main="", xlab="EVI", ylab="FPAR", xlim=c(0,1), ylim=c(0,1) )
        )
  abline( c(0,0), c(1,1), col="red" )
  legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=16, col=colorRampPalette( c("gray80", "navy", "red", "yellow"))(5), bty="n",  cex=0.8 )
dev.off()



