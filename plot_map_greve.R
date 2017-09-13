library(ncdf4)
library(fields)
library(sp)
library(maptools)
library(dplyr)

load( "data/overview_data_fluxnet2015_L5.Rdata" ) # loads 'overview'

# load("metainfo_FLUXNET2015.Rdata")
# overview <- overview %>% left_join( dplyr::select( siteinfo, mysitename, REFERENCE_PAPER, REFERENCE_DOI ), by="mysitename" )

load( "data/greve_fluxnet2015.Rdata" )  # df_greve
load( "data/ai_fluxnet2015.Rdata" )  # df_ai
load( "data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

overview <- overview %>% left_join( df_ai,          by="mysitename" )
overview <- overview %>% left_join( df_alpha,       by="mysitename" )

  magn <- 4
  ncols <- 3
  nrows <- 2
  widths <- rep(1.6*magn,ncols)
  widths[1] <- 0.35*widths[3]
  widths[2] <- 0.12*widths[3]
  heights <- rep(magn,nrows)*0.5
  order <- matrix( c(1,2,3,3,4,4), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  ncfiln <- "./data/greve/ep_over_p_cru_ncep.nc"
  if (!file.exists(ncfiln)) {
    epop <- array( 1, dim=c(720,360) )
  } else {
    nc <- nc_open( ncfiln )
    epop <- ncvar_get( nc, varid="EP_OVER_P_CRU_NCEP" )
  }

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  savepar_fig <- par()$fig

  plotfiln <- "./fig_nn_fluxnet2015/map_greve.pdf"
  print(paste("plotting", plotfiln))
  pdf( plotfiln, width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show(panel)

    ## Boxplots for clusters vs. aridity index and alpha  
    par( mar=c(3,4,1,0.5))
    boxplot( ai ~ finalcluster,    data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="aridity index", names=c("cDD", "cGR", "cLS", "cNA"), col=c("tomato", "royalblue3", "grey70", "grey70"), las=1 )
    text( 0.5, 2.55, "a)", font=2, cex=1.2, adj=c(0,0) )

    boxplot( alpha ~ finalcluster, data=dplyr::filter( overview, !is.na(finalcluster) ), ylab="annual mean AET/PET", names=c("cDD", "cGR", "cLS", "cNA"), col=c("tomato", "royalblue3", "grey70", "grey70"), las=1 )
    text( 0.5, 0.92, "b)", font=2, cex=1.2, adj=c(0,0) )

    ## Color key
    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    color <- rev( c( "royalblue3", "wheat", "tomato" ))
    lev <- c(0,0.2,0.3,0.4,0.6,0.7,1,1.3,1.6,2,2.5,3)
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=150 )

    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    image(
          seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
          1/epop,
          ylim=c(-60,85), 
          # zlim=range(lev), 
          yaxt="n", xaxt="n",
          col=out.mycolorbar$colors, breaks=out.mycolorbar$margins
          )
    map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )

    axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
    axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )

    axis( 4, at=lat.labels, lab=F, lwd=1.5 )
    axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )

    axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
    axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )

    axis( 3, at=lon.labels, lab=F, lwd=1.5 )
    axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )

    text( -170, 80, "c)", font=2, cex=1.2 )

    growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

    ## cluster 4: no low soil moisture data
    with( dplyr::filter( overview, finalcluster==4 ), points( lon, lat, col='red', pch=4, cex=0.8 ) )

    ## cluster 3: no reduction in LUE
    with( dplyr::filter( overview, finalcluster==3 ), points( lon, lat, col='springgreen', pch=3, cex=0.8 ) )
    
    ## cluster 1
    with( dplyr::filter( overview, finalcluster==1 ), points( lon, lat, col='black', pch=21, bg='tomato', cex=1.2 ) )

    ## cluster 2
    with( dplyr::filter( overview, finalcluster==2 ), points( lon, lat, col='black', pch=21, bg='royalblue3', cex=1.2 ) )
    

    # legend( "bottomleft", c("cDD", "cGR", "cLS", "cNA" ), pch=c(21,21,3,4), bty="n", cex=c(1.2,1.2,0.8,0.8), bg=c( "tomato", "royalblue3", "black", "black" ), col=c( "black", "black", "springgreen", "red" ) )
    legend( "left", c("cDD", "cGR", "cLS", "cNA" ), pch=c( 21, 21, 3, 4 ), bty="n", col=c( "black", "black", "springgreen", "red" ), pt.bg=c( "tomato", "royalblue3", "black", "black" ), cex=1.2, box.col="white"  )

    ## Inset 1
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v_orig <- v
    # v <- c( (v[1]+v[2])/2, v[2], (v[3]+v[4])/2, v[4] )
    v <- c( v[1], v[1]+0.14*v[2], v[3]+0*v[4], v[3]+0.3*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0) )
    image(
      seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
      1/epop,
      ylim=c(35,60), 
      xlim=c(-10,20),
      yaxt="n", xaxt="n",
      col=out.mycolorbar$colors, breaks=out.mycolorbar$margins
    )
    map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
    with( dplyr::filter( overview, finalcluster==4 ), points( lon, lat, col='red', pch=4, cex=0.6 ) )
    with( dplyr::filter( overview, finalcluster==3 ), points( lon, lat, col='springgreen', pch=3, cex=0.6 ) )
    with( dplyr::filter( overview, finalcluster==1 ), points( lon, lat, col='black', pch=21, bg='tomato', cex=1.2 ) )
    with( dplyr::filter( overview, finalcluster==2 ), points( lon, lat, col='black', pch=21, bg='royalblue3', cex=1.2 ) )
    box()
   
  dev.off()


 
