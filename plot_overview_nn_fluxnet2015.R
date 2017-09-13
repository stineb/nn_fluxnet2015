library(gplots)
library(graphics)
library(dplyr)

source( "analyse_modobs.R" )
source( "get_spi_spei.R" )

load( "data/overview_data_fluxnet2015_L5.Rdata" ) # loads 'overview'

load( "data/wtd_fluxnet2015.Rdata" )  # loads 'df_wtd'
load( "data/greve_fluxnet2015.Rdata" )  # df_greve
load( "data/ai_fluxnet2015.Rdata" )  # df_ai
load( "data/soilparams_fluxnet2015.Rdata" )  # df_soil
load( "data/wtd_degraaf_fluxnet2015.Rdata" )  # loads 'df_wtd_degraaf'
load( "data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

overview <- overview %>% left_join( df_wtd,         by="mysitename" )
overview <- overview %>% left_join( df_greve,       by="mysitename" )
overview <- overview %>% left_join( df_ai,          by="mysitename" )
overview <- overview %>% left_join( df_alpha,       by="mysitename" )
overview <- overview %>% left_join( df_soil,        by="mysitename" )
overview <- overview %>% left_join( df_wtd_degraaf, by="mysitename" ) %>% mutate( wtd_degraaf=wtd_degraaf/10 )


##------------------------------------------------
## Function definition for plotting the colorful table
##------------------------------------------------
plot_table <- function( table, colors, marginColor, col.names=colnames(table), row.names=rownames(table), main="", text.cex=1.0 ){

  plot(c(-1,ncol(table)),c(0,nrow(table)+1), type="n", xaxt="n", yaxt="n", xlab="",ylab="",main=main, bty="n")

  ## margin top
  for (c in 1:ncol(table)) {
    rect( c-1, nrow(table), c, nrow(table) + 1,    col=marginColor )
    text( c-.5,nrow(table) +.5,col.names[c], cex=text.cex )
  }

  ## margin left
  for (r in 1:nrow(table)) {
    rect(-1, r-1, 0, r, col=marginColor )
    text(-.5, r-.5,row.names[nrow(table) - r + 1], cex=text.cex )
  }

  ## all cells
  for (r in 1:nrow(table)){
    for (c in 1:ncol(table)) {
      rect( c-1,  r-1, c, r, col=as.character( colors[ nrow(table) - r + 1, c ] ) )
      text( c-.5, r-.5, format( table[ nrow(table) - r + 1, c ], digits=2), cex=text.cex )
    }
  }

}

col.names <- c( expression( paste( Delta, "GPP (%)")), 
                expression( paste( Delta, "GPP"[dr], " (%)")), 
                "dr (%)", 
                expression( paste( "fLUE"[0])), 
                expression( paste( "fLUE"[1])), 
                expression( paste( Delta, "EVI (%)")), 
                "AI", 
                "AET/PET", 
                expression( paste("WTD"[FMM])), 
                expression( paste("WTD"[DG])), 
                expression( paste("drain"[HWSD])), 
                expression( paste("AWC"[HWSD])), 
                expression( paste( "veg"[IGBP])) 
                )

sub <- dplyr::select( overview, mysitename, fgpp_total, fgpp_drought, fdroughtdays, fvar_min, fvar_max, evi_norm3, ai, alpha, wtd, wtd_degraaf, drain, awc, classid, alignedcluster, quadfitcluster, finalcluster )  # NNall_rsq, NNgood_rsq
rownames(sub) <- sub$mysitename
sub <- sub %>% dplyr::select( -mysitename )

sub$evi_norm3 <- 100 * (1.0 - sub$evi_norm3 )

##------------------------------------------------
## Overview plot
##------------------------------------------------
  ## Color data frame
  col <- sub
  ncols <- 5
  colors <- colorRampPalette( c("red", "orange", "springgreen3") )( ncols )

  lev <- c( 0, 1, ncols )
  out.mycolorbar <- mycolorbar( c("red","orange","springgreen3"), lev, orient="v", plot=FALSE )

  col$fgpp_total   <- cut( sub$fgpp_total,   breaks=c( -100, seq( from=0, to=40, length.out=(ncols-1)), 60  ), labels=rev(colors) )
  col$fgpp_drought <- cut( sub$fgpp_drought, breaks=c( -100, seq( from=0, to=50, length.out=(ncols-1)), 250 ), labels=rev(colors) )
  col$fdroughtdays <- cut( sub$fdroughtdays, breaks=c( seq( from=0,   to=80.0, length.out=(ncols)), 100), labels=rev(colors) )
  col$fvar_min     <- cut( sub$fvar_min    , breaks=c(seq( from=0, to=0.95, length.out=(ncols)), 2), labels=colors )
  col$fvar_max     <- cut( sub$fvar_max    , breaks=c( seq( from=0, to=0.95, length.out=(ncols)), 9999), labels=colors )
  # col$NNall_rsq    <- cut( sub$NNall_rsq   , breaks=seq( from=range( sub$NNall_rsq   , na.rm=TRUE )[1], to=range( sub$NNall_rsq   , na.rm=TRUE )[2], length.out=(ncols+1)), labels=colors )
  # col$NNgood_rsq   <- cut( sub$NNgood_rsq  , breaks=seq( from=range( sub$NNgood_rsq  , na.rm=TRUE )[1], to=range( sub$NNgood_rsq  , na.rm=TRUE )[2], length.out=(ncols+1)), labels=colors )
  col$evi_norm3       <- cut( sub$evi_norm3      , breaks=c( -100, seq( from=0, to=70, length.out=(ncols))), labels=rev(colors) )
  col$wtd          <- cut( sub$wtd         , breaks=c( seq( from=-0.01, to=30, length.out=(ncols)), 300), labels=rev(colors) )
  col$wtd_degraaf  <- cut( sub$wtd_degraaf , breaks=c( seq( from=-0.01, to=30, length.out=(ncols)), 300), labels=rev(colors) )
  col$ai           <- cut( sub$ai          , breaks=c( 0, 0.05, 0.25, 0.5, 0.75, 1, 2.7 ), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
  col$alpha        <- cut( sub$alpha       , breaks=c( 0, 0.05, 0.25, 0.5, 0.75, 1, 2.7 ), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
  col$drain        <- cut( sub$drain       , breaks=seq(0,7), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 7 ) )
  col$awc          <- cut( as.numeric( as.character(sub$awc)) , breaks=seq( 0, 150, by=25), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
  # col$greve        <- cut( sub$greve       , breaks=c( -0.01, 0.5, 1, 1.5, 2), labels=c("royalblue3", "wheat3", "wheat3", "tomato") )

  ##??Add IGBP colors 
  library(plotKML)
  data( worldgrids_pal )
  df_cols <- data.frame( classid=unique(sub$classid) )
  df_cols$col <- rep( NA, nrow(df_cols) )

  df_cols$col[ which(df_cols$classid=="MF") ]  <- col2hex("springgreen4")
  df_cols$col[ which(df_cols$classid=="ENF") ] <- col2hex("dodgerblue4")
  df_cols$col[ which(df_cols$classid=="GRA") ] <- col2hex("khaki3")
  df_cols$col[ which(df_cols$classid=="WSA") ] <- worldgrids_pal$IGBP["Woody savannas " ]
  df_cols$col[ which(df_cols$classid=="SAV") ] <- worldgrids_pal$IGBP["Savannas " ]
  df_cols$col[ which(df_cols$classid=="EBF") ] <- col2hex("forestgreen")
  df_cols$col[ which(df_cols$classid=="WET") ] <- col2hex("violetred")
  df_cols$col[ which(df_cols$classid=="DBF") ] <- col2hex("springgreen3")
  df_cols$col[ which(df_cols$classid=="OSH") ] <- col2hex("white")
  df_cols$col[ which(df_cols$classid=="CRO") ] <- worldgrids_pal$IGBP["Croplands " ]
  df_cols$col[ which(df_cols$classid=="CSH") ] <- worldgrids_pal$IGBP["Closed shrublands " ]

  # df_cols$col[ which(df_cols$classid=="MF") ]  <- worldgrids_pal$IGBP["Mixed forest " ]
  # df_cols$col[ which(df_cols$classid=="ENF") ] <- worldgrids_pal$IGBP["Evergreen Needleleaf forest " ]
  # df_cols$col[ which(df_cols$classid=="GRA") ] <- worldgrids_pal$IGBP["Grasslands " ]
  # df_cols$col[ which(df_cols$classid=="WSA") ] <- worldgrids_pal$IGBP["Woody savannas " ]
  # df_cols$col[ which(df_cols$classid=="SAV") ] <- worldgrids_pal$IGBP["Savannas " ]
  # df_cols$col[ which(df_cols$classid=="EBF") ] <- worldgrids_pal$IGBP["Evergreen Broadleaf forest " ]
  # df_cols$col[ which(df_cols$classid=="WET") ] <- worldgrids_pal$IGBP["Permanent wetlands " ] 
  # df_cols$col[ which(df_cols$classid=="DBF") ] <- worldgrids_pal$IGBP["Deciduous Broadleaf forest " ]
  # df_cols$col[ which(df_cols$classid=="OSH") ] <- worldgrids_pal$IGBP["Open shrublands " ]
  # df_cols$col[ which(df_cols$classid=="CRO") ] <- worldgrids_pal$IGBP["Croplands " ]
  # df_cols$col[ which(df_cols$classid=="CSH") ] <- worldgrids_pal$IGBP["Closed shrublands " ]

  tmp <- sub %>% left_join( df_cols, by="classid" )
  col$classid <- tmp$col

  pdf( "fig_nn_fluxnet2015/overview_cluster1.pdf", width=12, height=4 )
  # pdf( "fig_nn_fluxnet2015/overview_ALL.pdf", width=12, height=20 )
  par( mar=c(0.1, 0.1, 3, 0.1 ) )
  # par( mfrow=c(4,1) )
  idxs <- which( sub$finalcluster==1 ) 
  plot_table( sub[idxs,1:(ncol(sub)-3)], col[idxs,1:(ncol(col)-3)], "gray", main="cDD", text.cex=0.8, col.names=col.names )
  dev.off()

  pdf( "fig_nn_fluxnet2015/overview_cluster2.pdf", width=12, height=5.5 )
  par( mar=c(0.1, 0.1, 3, 0.1 ) )
  idxs <- which( sub$finalcluster==2 )
  plot_table( sub[idxs,1:(ncol(sub)-3)], col[idxs,1:(ncol(col)-3)], "gray", main="cGR", text.cex=0.8, col.names=col.names )
  dev.off()

  pdf( "fig_nn_fluxnet2015/overview_cluster3.pdf", width=12, height=7 )
  par( mar=c(0.1, 0.1, 3, 0.1 ) )
  idxs <- which( sub$finalcluster==3 )
  plot_table( sub[idxs,1:(ncol(sub)-3)], col[idxs,1:(ncol(col)-3)], "gray", main="cLS", text.cex=0.8, col.names=col.names )
  dev.off()

  pdf( "fig_nn_fluxnet2015/overview_cluster4.pdf", width=12, height=6.5 )
  par( mar=c(0.1, 0.1, 3, 0.1 ) )
  idxs <- which( sub$finalcluster==4 )
  plot_table( sub[idxs,1:(ncol(sub)-3)], col[idxs,1:(ncol(col)-3)], "gray", main="cNA", text.cex=0.8, col.names=col.names )
  dev.off()


##------------------------------------------------
## Some evaluations of the explanatory power of WTD:
##------------------------------------------------
  ## Agreement between WTD from Fan and DeGraaf dataset: very poor
  ## R2 of 0.22
  ##------------------------------------------------
  print( "Correlation of water table depth site data from Fan & Miguez-Macho (2013) compared with De Graaf et al. (2015):")
  out <- analyse_modobs( overview$wtd_degraaf*1e-1, overview$wtd, do.plot=FALSE )
  print( paste( "adjusted R2:", out$rsq ) )

  ## Correlation of fLUE_0 with WTD
  ##------------------------------------------------
  print( "Correlation of fLUE_0 with water table depth, site data from Fan & Miguez-Macho (2013)")
  out <- analyse_modobs( overview$y_x0, overview$wtd, do.plot=FALSE )
  print( paste( "adjusted R2:", out$rsq ) )

  print( "Correlation of fLUE_0 with water table depth, site data from De Graaf et al. (2015)")
  out <- analyse_modobs( overview$y_x0, overview$wtd_degraaf, do.plot=FALSE )
  print( paste( "adjusted R2:", out$rsq ) )


save( overview, file="data/overview_data_fluxnet2015_L6.Rdata" ) # loads 'overview'


