library(gplots)
library(graphics)
library(dplyr)

load( "data/overview_data_fluxnet2015_L5.Rdata" ) # loads 'overview'

load( "data/wtd_fluxnet2015.Rdata" )  # loads 'df_wtd'
load( "data/greve_fluxnet2015.Rdata" )  # df_greve
load( "data/ai_fluxnet2015.Rdata" )  # df_ai
load( "data/soilparams_fluxnet2015.Rdata" )  # df_soil
load( "data/wtd_degraaf_fluxnet2015.Rdata" )  # loads 'df_wtd_degraaf'
load( "data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

overview <- overview %>% left_join( df_wtd, by="mysitename" )
overview <- overview %>% left_join( df_greve, by="mysitename" )
overview <- overview %>% left_join( df_ai, by="mysitename" )
overview <- overview %>% left_join( df_alpha, by="mysitename" )
overview <- overview %>% left_join( df_soil, by="mysitename" )
overview <- overview %>% left_join( df_wtd_degraaf, by="mysitename" ) %>% mutate( wtd_degraaf=wtd_degraaf/10 )

# ## determine final clusters
# overview$finalcluster <- with( overview, ifelse( alignedcluster==4, 4, ifelse( quadfitcluster==2, 3, ifelse( alignedcluster==1, 1, ifelse( alignedcluster==2, 2, NA ) ) ) ) )

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

sub <- dplyr::select( overview, mysitename, fgpp_total, fgpp_drought, fdroughtdays, fvar_min, fvar_max, fapar3, ai, alpha, wtd, wtd_degraaf, drain, awc, classid, alignedcluster, quadfitcluster, finalcluster )  # NNall_rsq, NNgood_rsq
rownames(sub) <- sub$mysitename
sub <- sub %>% dplyr::select( -mysitename )

sub$fapar3 <- 100 * (1.0 - sub$fapar3 )

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
col$fdroughtdays <- cut( sub$fdroughtdays, breaks=seq( from=0,   to=100.0, length.out=(ncols+1)), labels=rev(colors) )
col$fvar_min     <- cut( sub$fvar_min    , breaks=c(seq( from=0, to=0.95, length.out=(ncols)), 2), labels=colors )
col$fvar_max     <- cut( sub$fvar_max    , breaks=c(seq( from=0, to=0.95, length.out=(ncols)), 2), labels=colors )
# col$NNall_rsq    <- cut( sub$NNall_rsq   , breaks=seq( from=range( sub$NNall_rsq   , na.rm=TRUE )[1], to=range( sub$NNall_rsq   , na.rm=TRUE )[2], length.out=(ncols+1)), labels=colors )
# col$NNgood_rsq   <- cut( sub$NNgood_rsq  , breaks=seq( from=range( sub$NNgood_rsq  , na.rm=TRUE )[1], to=range( sub$NNgood_rsq  , na.rm=TRUE )[2], length.out=(ncols+1)), labels=colors )
col$fapar3       <- cut( sub$fapar3      , breaks=c( -100, seq( from=0, to=70, length.out=(ncols))), labels=rev(colors) )
col$wtd          <- cut( sub$wtd         , breaks=c( seq( from=-0.01, to=70, length.out=(ncols)), 300), labels=rev(colors) )
col$wtd_degraaf  <- cut( sub$wtd_degraaf , breaks=c( seq( from=-0.01, to=70, length.out=(ncols)), 300), labels=rev(colors) )
col$ai           <- cut( sub$ai          , breaks=c( 0, 0.05, 0.25, 0.5, 0.75, 1, 2.7 ), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
col$alpha        <- cut( sub$alpha       , breaks=c( 0, 0.05, 0.25, 0.5, 0.75, 1, 2.7 ), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
col$drain        <- cut( sub$drain       , breaks=seq(0,7), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 7 ) )
col$awc          <- cut( as.numeric( as.character(sub$awc)) , breaks=seq( 0, 150, by=25), labels=colorRampPalette( c("red", "orange", "springgreen3") )( 6 ) )
# col$greve        <- cut( sub$greve       , breaks=c( -0.01, 0.5, 1, 1.5, 2), labels=c("royalblue3", "wheat3", "wheat3", "tomato") )

## Add IGBP colors 
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
plot_table( sub[idxs,1:(ncol(sub)-4)], col[idxs,1:(ncol(col)-4)], "gray", main="Cluster 1", text.cex=0.8, col.names=col.names )
dev.off()

pdf( "fig_nn_fluxnet2015/overview_cluster2.pdf", width=12, height=5.5 )
par( mar=c(0.1, 0.1, 3, 0.1 ) )
idxs <- which( sub$finalcluster==2 )
plot_table( sub[idxs,1:(ncol(sub)-4)], col[idxs,1:(ncol(col)-4)], "gray", main="Cluster 2", text.cex=0.8, col.names=col.names )
dev.off()

pdf( "fig_nn_fluxnet2015/overview_cluster3.pdf", width=12, height=6.3 )
par( mar=c(0.1, 0.1, 3, 0.1 ) )
idxs <- which( sub$finalcluster==3 )
plot_table( sub[idxs,1:(ncol(sub)-4)], col[idxs,1:(ncol(col)-4)], "gray", main="Cluster 3", text.cex=0.8, col.names=col.names )
dev.off()

pdf( "fig_nn_fluxnet2015/overview_cluster4.pdf", width=12, height=8.5 )
par( mar=c(0.1, 0.1, 3, 0.1 ) )
idxs <- which( sub$finalcluster==4 )
plot_table( sub[idxs,1:(ncol(sub)-4)], col[idxs,1:(ncol(col)-4)], "gray", main="Cluster 4", text.cex=0.8, col.names=col.names )
dev.off()


##------------------------------------------------
## boxplots by cluster
##------------------------------------------------
pdf( "fig_nn_fluxnet2015/boxplot_AI_by_cluster.pdf" )
boxplot( ai ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="aridity index", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4"), col="grey70", las=1 )
dev.off()


# boxplot( wtd ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="WTD (Fan et al.)", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4") )
# boxplot( wtd_degraaf ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="WTD (De Graaf et al.)", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4") )
# boxplot( drain ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="drainage class", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4") )

# overview$awc <- as.numeric( as.character( overview$awc ) )
# boxplot( awc ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="available water capacity", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4") )

# agg <- overview %>% group_by( finalcluster ) %>% summarise( drain=mean( drain ), awc=mean( awc ), ai=mean( ai ) )


##------------------------------------------------
## boxplot of WTD and alpha (AET/PET) by cluster
##------------------------------------------------
pdf( "fig_nn_fluxnet2015/boxplot_alpha_by_cluster.pdf" )
boxplot( alpha ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="annual mean AET/PET", names=c("cluster 1", "cluster 2", "cluster 3", "cluster 4"), col="grey70", las=1 )
dev.off()


##------------------------------------------------
## Agreement between WTD from Fan and DeGraaf dataset: very poor
## R2 of 0.22
##------------------------------------------------
print( "Correlation of water table depth site data from Fan & Miguez-Macho (2013) compared with De Graaf et al. (2015):")
out <- analyse_modobs( overview$wtd_degraaf*1e-1, overview$wtd, do.plot=FALSE )
print( paste( "adjusted R2:", out$rsq ) )


# ##------------------------------------------------
# ## regression tree
# ##------------------------------------------------
# library( rpart )
# # fit <- rpart( y_x0 ~ fapar3 + ai + drain + awc + wtd + classid, method="anova", data=dplyr::filter( overview, !is.na(finalcluster) ) )
# fit <- rpart( y_x0 ~ fapar3 + ai + classid, method="anova", data=dplyr::filter( overview, !is.na(finalcluster) ) )
# plot( fit, uniform=TRUE, main="Regression Tree")
# text( fit, use.n=TRUE, all=TRUE, cex=0.8 )

##------------------------------------------------
## correlation of fapar change and LUE reduction
##------------------------------------------------
growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

pdf( "fig_nn_fluxnet2015/corr_fvar_vs_dfapar.pdf", width=6, height=5.5 )
  par( las=1 )
  with( dplyr::filter( overview, finalcluster %in% c(1,2) ), 
        plot( fapar3, y_x0, pch=16, col='tomato', type='n', xlim=c(-0.2,1.2), ylim=c(-0.2,1.2), 
              axes=FALSE, xlab=expression(paste( "1 - EVI / EVI"[0] )), ylab=expression(paste( "fLUE"[0] ))
              )
    )
  xaxslab <- seq(1.2,-0.2,by=-0.2)
  xaxslab[7] <- 0.0 
  axis( 1, labels=xaxslab, at=seq(-0.2,1.2,by=0.2) )
  axis( 2 )
  box()

  linmod <- lm( y_x0 ~ fapar3, data=dplyr::filter( overview, finalcluster %in% c(1,2) ) )
  abline( linmod, col="black" )
  abline( c(0,0), c(1,1), col="black", lty=2 )

  ## cluster 1
  ## herbaceous
  with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$herb) ), points( fapar3, y_x0, pch=16, col='tomato', cex=1.0 ) )
  
  ## savannah
  with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$sav) ), points( fapar3, y_x0, pch=18, col='tomato', cex=1.2 ) )
  
  ## evergreen (woody)
  with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_evg) ), points( fapar3, y_x0, pch=17, col='tomato', cex=1.0 ) )

  ## wetland
  with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$wet) ), points( fapar3, y_x0, pch=25, col='tomato', bg='tomato' ) )


  ## cluster 2
  ## herbaceous
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$herb) ), points( fapar3, y_x0, pch=16, col='royalblue3', cex=1.0 ) )
  
  ## savannah
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$sav) ), points( fapar3, y_x0, pch=18, col='royalblue3', cex=1.2 ) )
  
  ## evergreen (woody)
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_evg) ), points( fapar3, y_x0, pch=17, col='royalblue3', cex=1.0 ) )

  ## wetland
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$wet) ), points( fapar3, y_x0, pch=25, col='royalblue3', bg='royalblue3' ) )

  ## shrublands
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$shrub) ), points( fapar3, y_x0, pch=8, col='royalblue3' ) )

  ## deciduous
  with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_dec) ), points( fapar3, y_x0, pch=15, col='royalblue3' ) )

  with( dplyr::filter( overview, finalcluster==1 ), text( fapar3+0.02, y_x0, mysitename, adj=c(0,0.5), col='tomato', cex=0.6 ) )
  with( dplyr::filter( overview, finalcluster==2 ), text( fapar3+0.02, y_x0, mysitename, adj=c(0,0.5), col='royalblue3', cex=0.6 ) )

  text( -0.15, 1.1, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=0.8 )

  text( c(-0.15, -0.15), c(0.9,0.82), c("cluster 1", "cluster 2"), col=c("tomato", "royalblue3"), cex=0.8,  adj=0.0 )

  legend( "left", c("herbaceous", "savannah", "woody evergreen", "wetlands", "shrublands", "woody deciduous"), pch=c(16,18,17,25,8,15), bty="n", cex=0.8, inset=c(0.05,0) )

dev.off()


##------------------------------------------------
## PLOT CLUSTERS ON A MAP
##------------------------------------------------
  growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

  if (nclust==5){
    cols <- c( "springgreen4", "royalblue3", "tomato", "orchid4", "goldenrod" )
  } else if (nclust==2) {
    cols <- c( "springgreen3", "red" )
  } else if (nclust==3) {
    cols <- c( "springgreen4", "tomato", "goldenrod" )
  }

  # library(ncdf4)
  library(fields)
  library(sp)
  library(maptools)
  library(rgeos)
  library(raster)

  greve <- raster( paste( myhome, "data/greve/arid_trans_hum_world.nc", sep="" ) )

  pdf( "fig_nn_fluxnet2015/cluster_map.pdf", width=8, height=6 )
    
    par(las=1)
    plot( greve , col=c( add_alpha("royalblue3", 0.7),  add_alpha("wheat3", 0.7),  add_alpha("tomato", 0.7)), ext=extent(c(-180,180,-65,85)), ylim=c(-65,85) )
    
    map( interior=FALSE, resolution=0, lwd=0.5, add=TRUE  ) # 

    ## cluster 4: no low soil moisture data
    with( dplyr::filter( overview, alignedcluster==4 ), points( lon, lat, col='black', bg='black', pch=4, cex=0.4 ) )

    ## cluster 3: no reduction in LUE
    with( dplyr::filter( overview, quadfitcluster==2 ), points( lon, lat, col='black', bg='black', pch=3, cex=0.4 ) )
    
    ## cluster 1
    with( dplyr::filter( overview, alignedcluster==1 & quadfitcluster!=2 & finalcluster!=4 ), points( lon, lat, col='black', pch=21, bg='tomato', cex=0.9 ) )

    ## cluster 2
    with( dplyr::filter( overview, alignedcluster==2 & quadfitcluster!=2 & finalcluster!=4 ), points( lon, lat, col='black', pch=21, bg='royalblue3', cex=0.9 ) )

  dev.off()


##------------------------------------------------
## Save to file
##------------------------------------------------
save( overview, file="data/overview_data_fluxnet2015_L6.Rdata" ) # loads 'overview'


