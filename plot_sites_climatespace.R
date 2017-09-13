library(gplots)
library(graphics)
library(dplyr)

source( paste( "analyse_modobs.R", sep="" ) )
source( paste( "get_spi_spei.R", sep="" ) )

load( "data/overview_data_fluxnet2015_L5.Rdata" ) # loads 'overview'

# load("metainfo_FLUXNET2015.Rdata")
# overview <- overview %>% left_join( dplyr::select( siteinfo, mysitename, REFERENCE_PAPER, REFERENCE_DOI ), by="mysitename" )

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

growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

##------------------------------------------------
## MULTIPANEL
##------------------------------------------------
  magn <- 4
  ncols <- 2
  nrows <- 2
  heights <- c(1,1)*magn
  widths <- c(1.2,1.2)*magn

  plotfiln <- "fig_nn_fluxnet2015/sites_climatespace.pdf"
  print(paste("plotting", plotfiln))
  pdf( plotfiln, width=sum(widths), height=sum(heights) )
    
    par( las=1, mar=c(4,4,1,1) )      
    panel <- layout(
                    matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                    # matrix( order, nrows, ncols, byrow=TRUE ),
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show( panel )

    ##------------------------------------------------
    ## correlation of GPP loss and alpha
    ##------------------------------------------------
      with( dplyr::filter( overview, finalcluster %in% c(1,2) ), 
            plot( alpha, fgpp_total, pch=16, col='tomato', type='n',
                  axes=TRUE, xlab="AET/PET", ylab="% GPP loss", xlim=c(0.2,1.1), ylim=c(-10,50)
                  )
        )

      linmod <- lm( fgpp_total ~ alpha, data=dplyr::filter( overview, finalcluster %in% c(1,2,3,4) ) )
      abline( linmod, col="black" )

      ## cluster 1
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$herb) ), points( alpha, fgpp_total, pch=16, col='tomato', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$sav) ), points( alpha, fgpp_total, pch=18, col='tomato', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_evg) ), points( alpha, fgpp_total, pch=17, col='tomato', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$wet) ), points( alpha, fgpp_total, pch=25, col='tomato', bg='tomato' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_dec) ), points( alpha, fgpp_total, pch=15, col='tomato' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$shrub) ), points( alpha, fgpp_total, pch=8, col='tomato' ) )


      ## cluster 2
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$herb) ), points( alpha, fgpp_total, pch=16, col='royalblue3', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$sav) ), points( alpha, fgpp_total, pch=18, col='royalblue3', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_evg) ), points( alpha, fgpp_total, pch=17, col='royalblue3', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$wet) ), points( alpha, fgpp_total, pch=25, col='royalblue3', bg='royalblue3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$shrub) ), points( alpha, fgpp_total, pch=8, col='royalblue3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_dec) ), points( alpha, fgpp_total, pch=15, col='royalblue3' ) )


      ## cluster 3
        ## herbaceous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$herb) ), points( alpha, fgpp_total, pch=16, cex=1.0, col='springgreen3' ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$sav) ), points( alpha, fgpp_total, pch=18, cex=1.2, col='springgreen3' ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_evg) ), points( alpha, fgpp_total, pch=17, cex=1.0, col='springgreen3' ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$wet) ), points( alpha, fgpp_total, pch=25, bg='springgreen3', col='springgreen3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$shrub) ), points( alpha, fgpp_total, pch=8, col='springgreen3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_dec) ), points( alpha, fgpp_total, pch=15, col='springgreen3' ) )


      ## cluster 4
        ## herbaceous
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$herb) ), points( alpha, fgpp_total, pch=16, cex=1.0, col='grey70' ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$sav) ), points( alpha, fgpp_total, pch=18, cex=1.2, col='grey70' ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$woody_evg) ), points( alpha, fgpp_total, pch=17, cex=1.0, col='grey70' ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$wet) ), points( alpha, fgpp_total, pch=25, bg='grey70', col='grey70' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$shrub) ), points( alpha, fgpp_total, pch=8, col='grey70' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$woody_dec) ), points( alpha, fgpp_total, pch=15, col='grey70' ) )

      with( dplyr::filter( overview, finalcluster==1 ), text( alpha+0.022, fgpp_total, mysitename, adj=c(0,0.5), col='tomato', cex=0.6 ) )
      with( dplyr::filter( overview, finalcluster==2 ), text( alpha+0.022, fgpp_total, mysitename, adj=c(0,0.5), col='royalblue3', cex=0.6 ) )
      # with( dplyr::filter( overview, finalcluster==3 ), text( alpha+0.022, fgpp_total, mysitename, adj=c(0,0.5), col='springgreen3', cex=0.6 ) )
      # with( dplyr::filter( overview, finalcluster==4 ), text( alpha+0.022, fgpp_total, mysitename, adj=c(0,0.5), col='grey70', cex=0.6 ) )

      text( 0.25, 45, "a)", font=2, adj=0.0, cex=1.2 )
      text( 0.90, 45, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
      cf <- coef(linmod) %>% round( 2 )
      eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
      text( 0.93, 40, eq )

    ##------------------------------------------------
    ## Sites in Budyko space
    ##------------------------------------------------
      with( dplyr::filter( overview, finalcluster %in% c(1,2,3,4) ), 
            plot( pet_o_p, aet_o_p, pch=16, col='tomato', type='n',
                  axes=TRUE, xlab="PET/P", ylab="AET/P", xlim=c(0.3,4), ylim=c(0.3,1.1)
                  )
        )

      lines( c(0,1), c(0,1), col="black", lty=2 )
      lines( c(1,10), c(1,1), col="black", lty=2 )


      ## cluster 1
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$herb) ), points( pet_o_p, aet_o_p, pch=16, col='tomato', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$sav) ), points( pet_o_p, aet_o_p, pch=18, col='tomato', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_evg) ), points( pet_o_p, aet_o_p, pch=17, col='tomato', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$wet) ), points( pet_o_p, aet_o_p, pch=25, col='tomato', bg='tomato' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_dec) ), points( pet_o_p, aet_o_p, pch=15, col='tomato' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$shrub) ), points( pet_o_p, aet_o_p, pch=8, col='tomato' ) )


      ## cluster 2
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$herb) ), points( pet_o_p, aet_o_p, pch=16, col='royalblue3', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$sav) ), points( pet_o_p, aet_o_p, pch=18, col='royalblue3', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_evg) ), points( pet_o_p, aet_o_p, pch=17, col='royalblue3', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$wet) ), points( pet_o_p, aet_o_p, pch=25, col='royalblue3', bg='royalblue3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$shrub) ), points( pet_o_p, aet_o_p, pch=8, col='royalblue3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_dec) ), points( pet_o_p, aet_o_p, pch=15, col='royalblue3' ) )


      ## cluster 3
        ## herbaceous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$herb) ), points( pet_o_p, aet_o_p, pch=16, cex=1.0, col='springgreen3' ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$sav) ), points( pet_o_p, aet_o_p, pch=18, cex=1.2, col='springgreen3' ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_evg) ), points( pet_o_p, aet_o_p, pch=17, cex=1.0, col='springgreen3' ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$wet) ), points( pet_o_p, aet_o_p, pch=25, bg='springgreen3', col='springgreen3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$shrub) ), points( pet_o_p, aet_o_p, pch=8, col='springgreen3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_dec) ), points( pet_o_p, aet_o_p, pch=15, col='springgreen3' ) )


      ## cluster 4
        ## herbaceous
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$herb) ), points( pet_o_p, aet_o_p, pch=16, cex=1.0, col='grey70' ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$sav) ), points( pet_o_p, aet_o_p, pch=18, cex=1.2, col='grey70' ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$woody_evg) ), points( pet_o_p, aet_o_p, pch=17, cex=1.0, col='grey70' ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$wet) ), points( pet_o_p, aet_o_p, pch=25, bg='grey70', col='grey70' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$shrub) ), points( pet_o_p, aet_o_p, pch=8, col='grey70' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster %in% c(4) & classid %in% (growtype$woody_dec) ), points( pet_o_p, aet_o_p, pch=15, col='grey70' ) )

        # with( dplyr::filter( overview, finalcluster==1 & pet_o_p > 1.0 ), text( pet_o_p+0.022, aet_o_p-0.02, mysitename, adj=c(0,0.5), col='tomato', cex=0.6 ) )
        # with( dplyr::filter( overview, finalcluster==2 & pet_o_p > 1.0 ), text( pet_o_p+0.022, aet_o_p-0.02, mysitename, adj=c(0,0.5), col='royalblue3', cex=0.6 ) )
        # with( dplyr::filter( overview, finalcluster==3 & pet_o_p > 1.0 ), text( pet_o_p+0.022, aet_o_p-0.02, mysitename, adj=c(0,0.5), col='springgreen3', cex=0.6 ) )
        # # with( dplyr::filter( overview, finalcluster==4 & pet_o_p > 1.2 ), text( pet_o_p+0.022, aet_o_p-0.02, mysitename, adj=c(0,0.5), col='grey70', cex=0.6 ) )

        text( 0.6, 1.04, "b)", font=2, adj=0.0, cex=1.2 )

        text( c(3, 3, 3, 3), seq( 0.7, 0.55, by=-0.04)-0.04, c("cDD", "cGR", "cLS", "cNA"), col=c("tomato", "royalblue3", "springgreen3", "grey70"), cex=0.8,  adj=0.0 )
        legend( "bottomright", c("herbaceous", "savannah", "woody evergreen", "wetlands", "shrublands", "woody deciduous"), pch=c(16,18,17,25,8,15), bty="n", cex=0.8, inset=c(0.05,0) )


    ##------------------------------------------------
    ## correlation of fapar change and LUE reduction
    ##------------------------------------------------
      growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

      with( dplyr::filter( overview, finalcluster %in% c(1,2) ), 
            plot( evi_norm3, y_x0, pch=16, col='tomato', type='n', xlim=c(0,1.3), ylim=c(-0.2,1.2), 
                  axes=FALSE, xlab=expression(paste( "1 - EVI / EVI"[0] )), ylab=expression(paste( "fLUE"[0] ))
                  )
        )
      xaxslab <- seq(1.2,-0.2,by=-0.2)
      xaxslab[7] <- 0.0 
      axis( 1, labels=xaxslab, at=seq(-0.2,1.2,by=0.2) )
      axis( 2 )
      box()

      linmod <- lm( y_x0 ~ evi_norm3, data=dplyr::filter( overview, finalcluster %in% c(1,2) ) )
      abline( linmod, col="black" )
      # abline( c(0,0), c(1,1), col="black", lty=2 )

      ## cluster 1
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$herb) ), points( evi_norm3, y_x0, pch=16, col='tomato', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$sav) ), points( evi_norm3, y_x0, pch=18, col='tomato', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_evg) ), points( evi_norm3, y_x0, pch=17, col='tomato', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$wet) ), points( evi_norm3, y_x0, pch=25, col='tomato', bg='tomato' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_dec) ), points( evi_norm3, y_x0, pch=15, col='tomato' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$shrub) ), points( evi_norm3, y_x0, pch=8, col='tomato' ) )


      ## cluster 2
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$herb) ), points( evi_norm3, y_x0, pch=16, col='royalblue3', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$sav) ), points( evi_norm3, y_x0, pch=18, col='royalblue3', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_evg) ), points( evi_norm3, y_x0, pch=17, col='royalblue3', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$wet) ), points( evi_norm3, y_x0, pch=25, col='royalblue3', bg='royalblue3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$shrub) ), points( evi_norm3, y_x0, pch=8, col='royalblue3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_dec) ), points( evi_norm3, y_x0, pch=15, col='royalblue3' ) )

      with( dplyr::filter( overview, finalcluster==1 ), text( evi_norm3+0.02, y_x0, mysitename, adj=c(0,0.5), col='tomato', cex=0.6 ) )
      with( dplyr::filter( overview, finalcluster==2 ), text( evi_norm3+0.02, y_x0, mysitename, adj=c(0,0.5), col='royalblue3', cex=0.6 ) )

      text( 0.05, 1.1, "c)", font=2, adj=0.0, cex=1.2 )
      text( 0.05, 1.0, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )


    ##------------------------------------------------
    ## correlation of alpha and LUE reduction
    ##------------------------------------------------
      with( dplyr::filter( overview, finalcluster %in% c(1,2,3) ), 
            plot( alpha, y_x0, pch=16, col='tomato', type='n', xlim=c(-0.2,1.2), ylim=c(-0.2,1.2), 
                  axes=FALSE, xlab=expression(paste( "AET/PET" )), ylab=expression(paste( "fLUE"[0] ))
                  )
        )

      axis( 1 )
      axis( 2 )
      box()

      linmod <- lm( y_x0 ~ alpha, data=dplyr::filter( overview, finalcluster %in% c(1,2,3) ) )
      abline( linmod, col="black" )
      # abline( c(0,0), c(1,1), col="black", lty=2 )

      ## cluster 1
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$herb) ), points( alpha, y_x0, pch=16, col='tomato', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$sav) ), points( alpha, y_x0, pch=18, col='tomato', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_evg) ), points( alpha, y_x0, pch=17, col='tomato', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$wet) ), points( alpha, y_x0, pch=25, col='tomato', bg='tomato' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$woody_dec) ), points( alpha, y_x0, pch=15, col='tomato' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==1 & classid %in% (growtype$shrub) ), points( alpha, y_x0, pch=8, col='tomato' ) )


      ## cluster 2
        ## herbaceous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$herb) ), points( alpha, y_x0, pch=16, col='royalblue3', cex=1.0 ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$sav) ), points( alpha, y_x0, pch=18, col='royalblue3', cex=1.2 ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_evg) ), points( alpha, y_x0, pch=17, col='royalblue3', cex=1.0 ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$wet) ), points( alpha, y_x0, pch=25, col='royalblue3', bg='royalblue3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$shrub) ), points( alpha, y_x0, pch=8, col='royalblue3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster==2 & classid %in% (growtype$woody_dec) ), points( alpha, y_x0, pch=15, col='royalblue3' ) )


      ## cluster 3
        ## herbaceous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$herb) ), points( alpha, y_x0, pch=16, cex=1.0, col='springgreen3' ) )
        
        ## savannah
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$sav) ), points( alpha, y_x0, pch=18, cex=1.2, col='springgreen3' ) )
        
        ## evergreen (woody)
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_evg) ), points( alpha, y_x0, pch=17, cex=1.0, col='springgreen3' ) )

        ## wetland
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$wet) ), points( alpha, y_x0, pch=25, bg='springgreen3', col='springgreen3' ) )

        ## shrublands
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$shrub) ), points( alpha, y_x0, pch=8, col='springgreen3' ) )

        ## deciduous
        with( dplyr::filter( overview, finalcluster %in% c(3) & classid %in% (growtype$woody_dec) ), points( alpha, y_x0, pch=15, col='springgreen3' ) )


      with( dplyr::filter( overview, finalcluster==1 ), text( alpha+0.02, y_x0, mysitename, adj=c(0,0.5), col='tomato', cex=0.6 ) )
      with( dplyr::filter( overview, finalcluster==2 ), text( alpha+0.02, y_x0, mysitename, adj=c(0,0.5), col='royalblue3', cex=0.6 ) )

      text( -0.15, 1.1, "d)", font=2, adj=0.0, cex=1.2 )
      text( -0.15, 1.0, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )


  dev.off()
