library(dplyr)

source( paste( workingdir, "/analyse_modobs.R", sep="" ) )
source( paste( workingdir, "/remove_outliers.R", sep="" ) )

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

filn <- paste( "data/nice_modis_agg_", char_wgt, nam_target, ".Rdata", sep="" )
if (file.exists(filn)){
  avl_modis <- TRUE
  load( filn ) # loads 'nice_to_modis_agg'
} else {
  avl_modis <- FALSE
}

filn <- paste( "data/nice_mte_agg_", char_wgt, nam_target, ".Rdata", sep="" )
if (file.exists(filn)){
  avl_mte <- TRUE
  load( filn ) # loads 'nice_to_mte_agg'
} else {
  avl_mte <- FALSE
}

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_modis_agg', 'df_dday_mte_agg', 

## add vegetation type info to nice_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )


##--------------------------------------
## ALPHA by fLUE drought
##--------------------------------------
pdf( paste( "fig_nn_fluxnet2015/boxplot_fluedrought_vs_alpha.pdf", sep=""), width=6, height=4 )

  par( las=1 )

  nice_agg <- nice_agg %>% mutate( alpha=aet_pmodel/pet_pmodel )

  ## alpha
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha) )$alpha )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha)  )$alpha )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, !is.na(alpha) & !is.infinite(alpha) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays  <- dplyr::filter( nice_agg, !is.na(alpha) & !is.infinite(alpha) ) %>% nrow()

  with( dplyr::filter( nice_agg, finalcluster %in% c(1,2) ), boxplot( alpha ~ is_drought_byvar, outline=FALSE, xlab="fLUE drought", ylab="AET/PET", col="grey70" ) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )

dev.off()


##--------------------------------------
## SPEI by fLUE drought
##--------------------------------------
pdf( paste( "fig_nn_fluxnet2015/boxplot_fluedrought_vs_alpha_spi_spei_1mo.pdf", sep=""), width=6, height=4, bg="white" )

  par( las=1, mfrow=c(1,3) )

  ## spi 1
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & finalcluster %in% c(1,2) & !is.na(spi1) & !is.infinite(spi1) )$spi1 )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & finalcluster %in% c(1,2) & !is.na(spi1) & !is.infinite(spi1)  )$spi1 )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(spi1) & !is.infinite(spi1) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays  <- dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(spi1) & !is.infinite(spi1) ) %>% nrow()

  with( dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(spi1) & !is.infinite(spi1) ), boxplot( spi1 ~ is_drought_byvar, xlab="fLUE drought", ylab="SPI, 1 mo.", col="grey70", outline=FALSE) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )
  mtext( "a)", adj=0, line=1, font=2 )


  ## spei 1
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & finalcluster %in% c(1,2) & !is.na(spei1) & !is.infinite(spei1) )$spei1 )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & finalcluster %in% c(1,2) & !is.na(spei1) & !is.infinite(spei1)  )$spei1 )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(spei1) & !is.infinite(spei1) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays <- dplyr::filter(  nice_agg, finalcluster %in% c(1,2) & !is.na(spei1) & !is.infinite(spei1) ) %>% nrow()

  with( dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(spei1) & !is.infinite(spei1) ), boxplot( spei1 ~ is_drought_byvar, xlab="fLUE drought", ylab="SPEI, 1 mo.", col="grey70", outline=FALSE) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )
  mtext( "b)", adj=0, line=1, font=2 )


  ## AET/PET
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha) )$alpha )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha)  )$alpha )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays  <- dplyr::filter( nice_agg, finalcluster %in% c(1,2) & !is.na(alpha) & !is.infinite(alpha) ) %>% nrow()

  with( dplyr::filter( nice_agg, finalcluster %in% c(1,2) ), boxplot( alpha ~ is_drought_byvar, outline=FALSE, xlab="fLUE drought", ylab="AET/PET", col="grey70" ) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )  
  mtext( "c)", adj=0, line=1, font=2 )

dev.off()


pdf( paste( "fig_nn_fluxnet2015/boxplot_fluedrought_vs_spi_spei_3mo.pdf", sep=""), width=6, height=4 )

  par( las=1, mfrow=c(1,2) )

  ## spei 1
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & !is.na(spei3) & !is.infinite(spei3) )$spei3 )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & !is.na(spei3) & !is.infinite(spei3)  )$spei3 )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, !is.na(spei3) & !is.infinite(spei3) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays <- dplyr::filter( nice_agg, !is.na(spei3) & !is.infinite(spei3) ) %>% nrow()

  with( dplyr::filter( nice_agg, !is.na(spei3) & !is.infinite(spei3) ), boxplot( spei3 ~ is_drought_byvar, xlab="fLUE drought", ylab="SPEI, 3 mo.", col="grey70", outline=FALSE) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )

  ## spi 1
  var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & !is.na(spi3) & !is.infinite(spi3) )$spi3 )
  var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & !is.na(spi3) & !is.infinite(spi3)  )$spi3 )
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
  wtest <- wilcox.test( var1, var2, na.action=na.omit )
  nsites <- dplyr::filter( nice_agg, !is.na(spi3) & !is.infinite(spi3) ) %>% dplyr::select( mysitename ) %>% unique() %>% nrow()
  ndays <- dplyr::filter( nice_agg, !is.na(spi3) & !is.infinite(spi3) ) %>% nrow()

  with( dplyr::filter( nice_agg, !is.na(spi3) & !is.infinite(spi3) ), boxplot( spi3 ~ is_drought_byvar, xlab="fLUE drought", ylab="SPI, 3 mo.", col="grey70", outline=FALSE) )
  mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  mtext( paste( "N =", as.character( ndays ) ), adj=1, line=1, cex=0.8 )
  mtext( paste( "sites =", as.character( nsites ) ), adj=1, line=2, cex=0.8 )

dev.off()

# ## soil moisture anomaly
# with( nice_agg, boxplot( soilm_mean_anom ~ is_drought_byvar, xlab="fLUE drought", ylab="soil moisture anomaly by DOY", outline=FALSE) )

# var1 <- unique( dplyr::filter( nice_agg, is_drought_byvar  & !is.na(soilm_mean_anom) )$soilm_mean_anom )
# var2 <- unique( dplyr::filter( nice_agg, !is_drought_byvar & !is.na(soilm_mean_anom)  )$soilm_mean_anom )

# ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )
# wtest <- wilcox.test( var1, var2, na.action=na.omit )

# with( nice_agg, heatscatter( fvar, spei3, pch=16, xlab="fLUE", ylab="SPEI, 3 mo."))


##--------------------------------------
## fLUE (rolling mean) by EVI extreme
##--------------------------------------
  var1 <- dplyr::filter( nice_agg, is_fapar_extreme  & finalcluster %in% c(1,2,3,4) )$fvar_rollmean
  var2 <- dplyr::filter( nice_agg, !is_fapar_extreme & finalcluster %in% c(1,2,3,4)  )$fvar_rollmean
  ttest <- t.test( var1, var2, paired=FALSE, na.action=na.omit )

  pdf( paste( "fig_nn_fluxnet2015/boxplot_fvar_by_eviextreme.pdf", sep=""), width=4, height=4 )
    par( las=1, mar=c(4,4,1,1) )
    boxplot( 
      fvar_rollmean ~ is_fapar_extreme, 
      data=dplyr::filter( nice_agg, finalcluster %in% c(1,2,3,4)), 
      ylim=c(0.75,1.2), 
      col="grey70",
      xlab="EVI extreme", ylab="fLUE (365 d moving average)"
      )
    # mtext( paste( "p =", format( ttest$p.value, digits=2 ) ), adj=1, cex=0.8 )
  dev.off()

##--------------------------------------
## MOD VS OBS OF ALL DATA AGGREGATED
##--------------------------------------

  ##--------------------------------------
  ## NN evaluation: LUE and GPP vs. obs.
  ## aggregated over soil moisture datasets and using fLUE drought identification for split
  ##--------------------------------------

    ## GPP
    print("plotting mod vs obs for GPP in nice_agg ...")
    pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_gpp_rct_ALL_FROMNICE", char_wgt, ".pdf", sep=""), width=8, height=8 )
      par( mfrow=c(2,2) )
      stats_tmp <- analyse_modobs( 
                                  nice_agg$gpp_nn_act, 
                                  nice_agg$gpp_obs, 
                                  plot.title=expression( paste("NN"[act], "  all days")), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65      
                                  )
      mtext( "a)", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_obs, 
                                  plot.title=expression( paste("NN"[pot], "   moist days")), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "b)", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, is_drought_byvar )$gpp_obs, 
                                  plot.title=expression( paste("NN"[pot], "   dry days")), 
                                  plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "c)", line = 1, font = 2, adj=c(0,0) )

      stats_tmp <- analyse_modobs( 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_pot, 
                                  dplyr::filter( nice_agg, !is_drought_byvar )$gpp_nn_act, 
                                  plot.title=expression( paste("NN"[pot], " vs. NN"[act], "  moist days")), 
                                  plot.xlab=expression(paste("actual GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  plot.ylab=expression(paste("potential GPP (gC m"^{-2}, " d"^{-1}, ")")),
                                  lab.xpos=0.65
                                  )
      mtext( "d)", line = 1, font = 2, adj=c(0,0) )

    dev.off()

    ## LUE
    print("plotting mod vs obs for LUE in nice_agg ...")
    pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_lue_rct_ALL_FROMNICE", char_wgt, ".pdf", sep=""), width=8, height=8 )
      par( mfrow=c(2,2) )
      stats_tmp <- analyse_modobs( 
                                  nice_agg$var_nn_act, 
                                  nice_agg$lue_obs_evi, 
                                  plot.title=expression( paste("NN"[act], "  all days")), 
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


  # ##--------------------------------------
  # ## EXPANDED BY SOIL MOISTURE DATASET
  # ##--------------------------------------
  #  ## GPP
  #   print("plotting mod vs obs for GPP in nice_resh ...")
  #   pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_gpp_rct_ALL_FROMNICE_RESH", char_wgt, ".pdf", sep=""), width=8, height=8 )
  #     par( mfrow=c(2,2) )
  #     with( nice_resh,
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_act * iabs, 
  #                                   lue_obs_evi * iabs, 
  #                                   plot.title=expression( paste("NN"[pot])), 
  #                                   plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   lab.xpos=0.65      
  #                                   )
        
  #       )
  #     mtext( "d", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot * iabs, 
  #                                   lue_obs_evi * iabs, 
  #                                   plot.title=expression( paste("NN"[pot], "   moist days")), 
  #                                   plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "c", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, !moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot * iabs, 
  #                                   lue_obs_evi * iabs, 
  #                                   plot.title=expression( paste("NN"[pot], "   dry days")), 
  #                                   plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "b", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot * iabs, 
  #                                   var_nn_act * iabs, 
  #                                   plot.title=expression( paste("NN"[pot], " vs. NN"[act], "  moist days")), 
  #                                   plot.xlab=expression(paste("actual GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("potential GPP (gC m"^{-2}, " d"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "a", line = 1, font = 2, adj=c(0,0) )

  #   dev.off()

  #   ## LUE
  #   print("plotting mod vs obs for LUE in nice_resh ...")
  #   pdf( paste( "fig_nn_fluxnet2015/modobs/modobs_lue_rct_ALL_FROMNICE_RESH", char_wgt, ".pdf", sep=""), width=8, height=8 )
  #     par( mfrow=c(2,2) )
  #     with( nice_resh,
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_act, 
  #                                   lue_obs_evi, 
  #                                   plot.title=expression( paste("NN"[pot])), 
  #                                   plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
  #                                   lab.xpos=0.65      
  #                                   )
        
  #       )
  #     mtext( "d", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot, 
  #                                   lue_obs_evi, 
  #                                   plot.title=expression( paste("NN"[pot], "   moist days")), 
  #                                   plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "c", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, !moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot, 
  #                                   lue_obs_evi, 
  #                                   plot.title=expression( paste("NN"[pot], "   dry days")), 
  #                                   plot.xlab=expression(paste("observed LUE (gC mol"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("predicted LUE (gC mol"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "b", line = 1, font = 2, adj=c(0,0) )

  #     with( dplyr::filter( nice_resh, moist ),
  #       stats_tmp <- analyse_modobs( 
  #                                   var_nn_pot, 
  #                                   var_nn_act, 
  #                                   plot.title=expression( paste("NN"[pot], " vs. NN"[act], "  moist days")), 
  #                                   plot.xlab=expression(paste("actual LUE (gC mol"^{-1}, ")")),
  #                                   plot.ylab=expression(paste("potential LUE (gC mol"^{-1}, ")")),
  #                                   lab.xpos=0.65
  #                                   )
        
  #       )
  #     mtext( "a", line = 1, font = 2, adj=c(0,0) )

  #   dev.off()


##------------------------------------------------
## EVI vs. FPAR
##------------------------------------------------
  pdf("fig_nn_fluxnet2015/fpar_vs_evi.pdf")
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
