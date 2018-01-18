library(dplyr)
library(Hmisc)

source( "predict_nn.R" )

load( "data/overview_data_fluxnet2015_L5.Rdata" ) # loads 'overview'
do.sites <- filter( overview, finalcluster %in% c(1,2) & !(mysitename=="AU-Wom") ) %>% select( mysitename )
do.sites <- do.sites$mysitename

overwrite = TRUE
makepdf = TRUE

df_vpdtest_agg <- c()
df_soilmtest_agg <- c()
hilo_agg <- c()

for (sitename in do.sites){

  ##------------------------------------------------
  ## Get data
  ##------------------------------------------------
  filn <- paste( "data/aligned_", sitename, ".Rdata", sep="" )

  if (file.exists(filn)){

    ## aligned data
    load( filn )    # loading 'df_dday'

    ##------------------------------------------------
    ## Get NN models
    ##------------------------------------------------
    outfiln <- paste0( "./data/df_", sitename, ".Rdata" )

    if (file.exists(outfiln) & !overwrite){

      load( outfiln )
   
    } else {

      outdir <- "./data/profile/"
      filn <- paste( outdir, "profile_light_lue_obs_evi_nn_", sitename, ".Rdata", sep="" )

      deleteagain <- FALSE
      if (!file.exists(filn)){
        print( "downloading profile file ..." )
        system( paste("rsync -avz bstocker@login.cx1.hpc.ic.ac.uk:/home/bstocker/data/nn_fluxnet/profile/profile_light_lue_obs_evi_nn_", sitename, ".Rdata"," ./data/profile/", sep="" ) )
        # deleteagain <- TRUE
        print("done.")
      }

      load( filn )

      isoilm_data <- "soilm_swbm"
      ipackage <- "nnet"
      trh <- profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$best_soilm_trh[1]
      nn_act <- profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$nn_act
      nn_vpd <- profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$nn_vpd
      nn_pot <- profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]][[ paste0( "smtrh_", trh ) ]]$nn_pot

      ##------------------------------------------------
      ## Sensitivity: Create data frames to use as NN input
      ##------------------------------------------------
      predictors_act <- nn_act$coefnames
      predictors_vpd <- nn_vpd$coefnames
      predictors_pot <- nn_pot$coefnames
      med <- df_dday  %>% group_by( infvarbin ) %>%
                          summarise( ppfd       = median( ppfd, na.rm=TRUE ),
                                     temp       = median( temp, na.rm=TRUE ),
                                     vpd        = median( vpd,  na.rm=TRUE ),
                                     soilm_swbm = median( soilm_mean, na.rm=TRUE )
                                   )

      ibin <- which( med$infvarbin=="(-20,0]" )
      nsteps <- 100
      df_vpdtest <-  data.frame( mysitename = rep( sitename, nsteps ),
                                      ppfd       = rep( med$ppfd[ibin], nsteps ),
                                      temp       = rep( med$temp[ibin], nsteps ),
                                      vpd        = seq(0, 3000, 3000/(nsteps-1)),
                                      # soilm_swbm = rep( 1, nsteps ),
                                      soilm_swbm = rep( med$soilm_swbm[ibin], nsteps ),
                                      ndatapoints= length( profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$year_dec )
                                    )

      df_soilmtest <-  data.frame( mysitename = rep( sitename, nsteps ),
                                        ppfd       = rep( med$ppfd[ibin], nsteps ),
                                        temp       = rep( med$temp[ibin], nsteps ),
                                        vpd        = rep( med$vpd[ibin], nsteps),
                                        soilm_swbm = seq(1, 0, -1/(nsteps-1)),
                                        ndatapoints= length( profile_nn_light[[ sitename ]][[ isoilm_data ]][[ ipackage ]]$year_dec )
                                      )

      ##------------------------------------------------
      ## Predict values for sensitivity to VPD 
      ##------------------------------------------------
      ## NN_pot
      df_vpdtest$lue_pot <- predict_nn( 
                                    data       = df_vpdtest, 
                                    weights    = rep( 1, nsteps ),
                                    predictors = predictors_pot, 
                                    nam_target = "", 
                                    nn         = nn_pot, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals
      df_vpdtest$nlue_pot <- df_vpdtest$lue_pot / df_vpdtest$lue_pot[1]

      ## NN_vpd
      df_vpdtest$lue_vpd <- predict_nn( 
                                    data       = df_vpdtest, 
                                    weights    = rep( 1, nsteps ),
                                    predictors = predictors_vpd, 
                                    nam_target = "", 
                                    nn         = nn_vpd, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals
      df_vpdtest$nlue_vpd <- df_vpdtest$lue_vpd / df_vpdtest$lue_vpd[1]

      ## NN_act
      df_vpdtest$lue_act <- predict_nn( 
                                    data       = df_vpdtest, 
                                    weights    = rep( 1, nsteps ),  
                                    predictors = predictors_act, 
                                    nam_target = "", 
                                    nn         = nn_act, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals
      df_vpdtest$nlue_act <- df_vpdtest$lue_act / df_vpdtest$lue_act[1]

      ##------------------------------------------------
      ## Predict values for sensitivity to soil moisture 
      ##------------------------------------------------
      ## NN_act
      df_soilmtest$lue_act <- predict_nn( 
                                    data       = df_soilmtest, 
                                    weights    = rep( 1, nsteps ),  
                                    predictors = predictors_act, 
                                    nam_target = "", 
                                    nn         = nn_act, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals  

      ## save intermediate step
      save( df_vpdtest, df_soilmtest, file=outfiln )   


      ##------------------------------------------------
      ## Predict hihi and lolo values
      ##------------------------------------------------
      ## all data
      filn <- paste( "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_lue_obs_evi.Rdata", sep="" )
      load( filn ) ## gets list 'nn_fluxnet'
      nice <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )  

      nice <- nice %>% dplyr::select( ppfd, temp, vpd, soilm_swbm, lue_obs_evi )

      ## create input data frames hihi and lolo
      hihi <- nice %>% dplyr::filter( soilm_swbm > 0.75 ) 
      vpd_trh <- quantile( hihi$vpd, 0.9 )
      hihi <- hihi %>% dplyr::filter( vpd > vpd_trh )

      lolo <- nice %>% dplyr::filter( soilm_swbm < 0.25 ) 
      vpd_trh <- quantile( lolo$vpd, 0.1 )
      lolo <- lolo %>% dplyr::filter( vpd < vpd_trh )      

      ## predict hihi values
      ## NN_pot
      hihi$lue_pot <- predict_nn( 
                                    data       = hihi, 
                                    weights    = rep( 1, nrow(hihi) ),
                                    predictors = predictors_pot, 
                                    nam_target = "", 
                                    nn         = nn_pot, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      ## NN_vpd
      hihi$lue_vpd <- predict_nn( 
                                    data       = hihi, 
                                    weights    = rep( 1, nrow(hihi) ),
                                    predictors = predictors_vpd, 
                                    nam_target = "", 
                                    nn         = nn_vpd, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      ## NN_act
      hihi$lue_act <- predict_nn( 
                                    data       = hihi, 
                                    weights    = rep( 1, nrow(hihi) ),  
                                    predictors = predictors_act, 
                                    nam_target = "", 
                                    nn         = nn_act, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      hihi <- hihi %>% mutate( hihi=TRUE )

      ## predict lolo values
      ## NN_pot
      lolo$lue_pot <- predict_nn( 
                                    data       = lolo, 
                                    weights    = rep( 1, nrow(lolo) ),
                                    predictors = predictors_pot, 
                                    nam_target = "", 
                                    nn         = nn_pot, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      ## NN_vpd
      lolo$lue_vpd <- predict_nn( 
                                    data       = lolo, 
                                    weights    = rep( 1, nrow(lolo) ),
                                    predictors = predictors_vpd, 
                                    nam_target = "", 
                                    nn         = nn_vpd, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      ## NN_act
      lolo$lue_act <- predict_nn( 
                                    data       = lolo, 
                                    weights    = rep( 1, nrow(lolo) ),  
                                    predictors = predictors_act, 
                                    nam_target = "", 
                                    nn         = nn_act, 
                                    do_predict = TRUE, 
                                    package    = ipackage
                                    )$vals

      lolo <- lolo %>% mutate( hihi=FALSE )
      hilo <- rbind( hihi, lolo )

      hilo_agg <- rbind( hilo_agg, hilo )

      # boxplot( log( lue_vpd / lue_obs_evi ) ~ hihi , data=hilo, main="NN_vpd", outline=FALSE ); abline( h=0, col='red' )
      # boxplot( log( lue_act / lue_obs_evi ) ~ hihi , data=hilo, main="NN_act", outline=FALSE ); abline( h=0, col='red' )
      # boxplot( log( lue_pot / lue_obs_evi ) ~ hihi , data=hilo, main="NN_pot", outline=FALSE ); abline( h=0, col='red' )
      
      ##------------------------------------------------
      ## Plot hihi and lolo
      ##------------------------------------------------
      if (nrow(hihi)>0 && nrow(lolo)>0){

        if (makepdf) pdf( paste( "./fig_nn_fluxnet2015/ratio_vpdtest/ratio_vpdtest2_lue_obs_evi_", sitename, ".pdf", sep="" ), width=8, height=6 )      
          
          par( las=1, xaxs="i" )

          ylim <- c(-1.2,2.6)

          ## NNpot
          bp1 <- boxplot( 
            log( lue_pot / lue_obs_evi ) ~ hihi, 
            data  = hilo,
            col   = c("grey50"), 
            ylim  = ylim, 
            # xlim=c(min(soilm_threshold), max(soilm_threshold)+0.025), 
            las = 1, 
            xlab = "", 
            ylab = "log of modeled/observed", 
            yaxs = "i",
            at   = c(5,1) ,
            xlim = c(0,8),
            outline=FALSE,
            axes=FALSE,
            main=sitename
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
            log( lue_vpd / lue_obs_evi ) ~ hihi, 
            data  = hilo, 
            col   =c("grey70"), 
            xlab  ="", 
            ylab  ="", 
            yaxs  ="i",
            at   = c(5,1)+1 ,
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
            log( lue_act / lue_obs_evi ) ~ hihi, 
            data  = hilo, 
            col  = c("springgreen"), 
            xlab  ="", 
            ylab  ="", 
            yaxs  ="i",
            at   = c(5,1)+2 ,
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

          axis( 1, at=c(2,6), labels=c("high VPD, high soil moisture", "low VPD, low soil moisture"), tick=FALSE, cex.axis=0.8 )
          # axis( 1 )
          axis( 2 )
          box()        
          abline( h=0.0, col='red' )
          rect( 0, -2, 4, 4, border=NA, col=rgb(0,0,0,0.2) )

        if (makepdf) dev.off()

      }

    }

    ## add to dataframe holding all sites' data
    df_vpdtest_agg   <- rbind( df_vpdtest_agg, df_vpdtest )
    df_soilmtest_agg <- rbind( df_vpdtest_agg, df_vpdtest )

    # if (deleteagain){
    #   system( paste("rm /alphadata01/bstocker/data/nn_fluxnet/profile/profile_lue_obs_evi_nn_", sitename, ".Rdata", sep="" ) )
    # }


    ##------------------------------------------------
    ## Plot values for sensitivity analysis
    ##------------------------------------------------
    pdf( paste0( "fig_nn_fluxnet2015/vpd_sensitivity/vpd_sensitivity_", sitename, ".pdf" ) )
      par( las=1 )
      with( df_vpdtest, plot(  vpd, nlue_vpd, pch=16, type='l', col='tomato', xlab="VPD (Pa)", main=sitename, ylab=expression( paste( "LUE/LUE"[0] ) ) ) )
      with( df_vpdtest, lines( vpd, nlue_pot, pch=16, type='l', col='royalblue3' ) )
      with( df_vpdtest, lines( vpd, nlue_act, pch=16, type='l', col='black' ) )
      legend( "bottomleft", c("NN_act", "NN_pot", "NN_vpd"), col=c("black", "royalblue3", "tomato"), lty=1, bty="n" )
    dev.off()

    # soilm_lev <- seq(0, 1, by=1/(10-1))
    # for (itest in 1:100){
    #   tmp <- inp_med %>% mutate( soilm_swbm = rep( soilm_lev[itest], nsteps ) )
    #   tmp$lue_act <- predict_nn( 
    #                               data       = tmp, 
    #                               weights    = rep( 1, nsteps ),  
    #                               predictors = c( predictors, isoilm_data ), 
    #                               nam_target = "", 
    #                               nn         = nn_act, 
    #                               do_predict = TRUE, 
    #                               package    = ipackage
    #                               )$vals
    #   with( tmp, lines( vpd, lue_act/lue_act[1], pch=16, type='l', col=add_alpha('black', 0.1 ) ) )
    # }

    # with( df_soilmtest, plot(  soilm_swbm, lue_act/lue_act[1], pch=16, type='l', col='tomato', xlab="fractional soil water (unitless)", ylab=expression( paste( "LUE/LUE"[0] ) ) ) )

  }

}

## Save aggregated data
save( df_vpdtest_agg, df_soilmtest_agg, file="./data/df_vpdtest_agg.Rdata" )
save( hilo_agg, file="./data/hilo_agg.Rdata" )

load( "./data/df_vpdtest_agg.Rdata" )
load( "./data/hilo_agg.Rdata" )

## aggregated plots
df_vpdtest_agg_aggbyvpd <- df_vpdtest_agg %>% group_by( vpd ) %>% 
                                                        summarise( 
                                                                    lue_pot_mean = wtd.mean( nlue_pot, weights=ndatapoints ),
                                                                    lue_vpd_mean = wtd.mean( nlue_vpd, weights=ndatapoints ),
                                                                    lue_act_mean = wtd.mean( nlue_act, weights=ndatapoints ),

                                                                    lue_pot_median = wtd.quantile( nlue_pot, weights=ndatapoints, probs=0.5 ),
                                                                    lue_vpd_median = wtd.quantile( nlue_vpd, weights=ndatapoints, probs=0.5 ),
                                                                    lue_act_median = wtd.quantile( nlue_act, weights=ndatapoints, probs=0.5 ),
                                                                    
                                                                    lue_pot_q40 = wtd.quantile( nlue_pot, weights=ndatapoints, probs=0.40 ),
                                                                    lue_vpd_q40 = wtd.quantile( nlue_vpd, weights=ndatapoints, probs=0.40 ),
                                                                    lue_act_q40 = wtd.quantile( nlue_act, weights=ndatapoints, probs=0.40 ),
                                                                    
                                                                    lue_pot_q60 = wtd.quantile( nlue_pot, weights=ndatapoints, probs=0.60 ),
                                                                    lue_vpd_q60 = wtd.quantile( nlue_vpd, weights=ndatapoints, probs=0.60 ),
                                                                    lue_act_q60 = wtd.quantile( nlue_act, weights=ndatapoints, probs=0.60 )
                                                                    )
par( las=1 )
with( df_vpdtest_agg_aggbyvpd, plot(  vpd, lue_vpd_mean, pch=16, type='l', col='tomato', xlab="VPD (Pa)", main="aggregated", ylab=expression( paste( "LUE/LUE"[0] ) ) ) )
with( df_vpdtest_agg_aggbyvpd, lines( vpd, lue_pot_mean, pch=16, type='l', col='royalblue3' ) )
with( df_vpdtest_agg_aggbyvpd, lines( vpd, lue_act_mean, pch=16, type='l', col='black' ) )

## Plot sensitivity, aggregated data

magn <- 4
ncols <- 2
nrows <- 1
heights <- c(1)*magn
widths  <- c(1.2,1.2)*magn

plotfiln <- "./fig_nn_fluxnet2015/vpd_soilmoisture_test.pdf"
if (makepdf) print(paste("plotting figure", plotfiln))
if (makepdf) pdf( plotfiln, width=sum(widths), height=sum(heights) )
  
  panel <- layout(
                matrix( c(1:(nrows*ncols)), nrows, ncols, byrow=TRUE ),
                widths=widths,
                heights=heights,
                TRUE
                )
  # layout.show( panel )

  ## 1st panel: sensitivity curve
  par( las=1, mar=c(4,4,2,1) )

  with( df_vpdtest_agg_aggbyvpd, plot(  vpd, lue_vpd_median, pch=16, type='l', col='tomato', xlab="VPD (Pa)", main="", ylab=expression( paste( "LUE/LUE"[0] ) ) ) )
  with( df_vpdtest_agg_aggbyvpd, lines( vpd, lue_pot_median, pch=16, type='l', col='royalblue3' ) )
  with( df_vpdtest_agg_aggbyvpd, lines( vpd, lue_act_median, pch=16, type='l', col='black' ) )

  with( df_vpdtest_agg_aggbyvpd, polygon( c(vpd, rev(vpd) ), c( lue_vpd_q40, rev(lue_vpd_q60) ), border=NA, col=add_alpha('tomato', 0.5) ) )
  with( df_vpdtest_agg_aggbyvpd, polygon( c(vpd, rev(vpd) ), c( lue_pot_q40, rev(lue_pot_q60) ), border=NA, col=add_alpha('royalblue3', 0.5) ) )
  with( df_vpdtest_agg_aggbyvpd, polygon( c(vpd, rev(vpd) ), c( lue_act_q40, rev(lue_act_q60) ), border=NA, col=add_alpha('black', 0.5) ) )

  legend( "bottomleft", c( expression( paste("NN"[act])),  expression( paste("NN"[pot])),  expression( paste("NN"[vpd]))), col=c("black", "royalblue3", "tomato"), lty=1, bty="n", lwd=2 )
  mtext( "a)", font=2, line=0.5, adj=0 )

  ## 2nd panel: boxplot
  par( las=1, xaxs="i" )

  ylim <- c(-1.2,2.6)

  ## NNpot
  bp1 <- boxplot( 
    log( lue_pot / lue_obs_evi ) ~ hihi, 
    data  = hilo_agg,
    col   = c("grey50"), 
    ylim  = ylim, 
    # xlim=c(min(soilm_threshold), max(soilm_threshold)+0.025), 
    las = 1, 
    xlab = "", 
    ylab = "log of modeled/observed", 
    yaxs = "i",
    at   = c(5,1) ,
    xlim = c(0,8),
    outline=FALSE,
    axes=FALSE,
    main=""
    )
  # text( 0.1, 2.7-0.5-0.03, expression(paste("R"^2)), adj=c(0,0))
  # text( 0.1, 2.5-0.5-0.03, expression(paste("RMSE")), adj=c(0,0))
  # text( 0.1, 2.7-0.5-0.03, "p-value", adj=c(0,0) )

  text( 1,   2.9-0.5, expression(paste("LUE"[pot])), cex=0.8 )
  text( 1+4, 2.9-0.5, expression(paste("LUE"[pot])), cex=0.8 )
  mtext( "b)", font=2, line=0.5, adj=0 )

  # # text( 1,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_dry   ) )
  # # text( 1+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_pot_moist ) )
  # text( 1,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_dry   ) )
  # text( 1+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_pot_moist   ) )

  # # text( 1,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_dry   ) )
  # text( 1+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_pot_moist ) )

  ## NNvpd
  bp2 <- boxplot( 
    log( lue_vpd / lue_obs_evi ) ~ hihi, 
    data  = hilo_agg, 
    col   =c("grey70"), 
    xlab  ="", 
    ylab  ="", 
    yaxs  ="i",
    at   = c(5,1)+1 ,
    # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
    add   = TRUE,
    axes  = FALSE,
    outline=FALSE
    )
  text( 2,   2.9-0.5, expression(paste("LUE"[VPD])), cex=0.8 )
  text( 2+4, 2.9-0.5, expression(paste("LUE"[VPD])), cex=0.8 )

  # # text( 2,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_dry   ) )
  # # text( 2+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_vpd_moist ) )
  # text( 2,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_dry   ) )
  # text( 2+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_vpd_moist   ) )

  # # text( 2,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_dry   ) )
  # # text( 2+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_vpd_moist ) )

  ## NNact
  bp3 <- boxplot( 
    log( lue_act / lue_obs_evi ) ~ hihi, 
    data  = hilo_agg, 
    col  = c("springgreen"), 
    xlab  ="", 
    ylab  ="", 
    yaxs  ="i",
    at   = c(5,1)+2 ,
    # xlim  =range(soilm_thrsh_avl)+c(-0.02,0.02),
    add   = TRUE,
    axes  = FALSE,
    outline=FALSE
    ) 
  text( 3,   2.9-0.5, expression(paste("LUE"[act])), cex=0.8 ) 
  text( 3+4, 2.9-0.5, expression(paste("LUE"[act])), cex=0.8 ) 

  # # text( 3,   2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_dry   ) )
  # # text( 3+4, 2.7-0.5, as.character(stats_vpdtest$rsq_nn_act_moist ) )
  # text( 3,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_dry   ) )
  # text( 3+4,   2.7-0.5, as.character(stats_vpdtest$pval_nn_act_moist   ) )

  # # text( 3,   2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_dry   ) )
  # # text( 3+4, 2.5-0.5, as.character(stats_vpdtest$rmse_nn_act_moist ) )

  axis( 1, at=c(2,6), labels=c("high VPD, high soil moisture", "low VPD, low soil moisture"), tick=FALSE, cex.axis=0.8 )
  # axis( 1 )
  axis( 2 )
  box()        
  abline( h=0.0, col='red' )
  rect( 0, -2, 4, 4, border=NA, col=rgb(0,0,0,0.2) )

if (makepdf) dev.off()

