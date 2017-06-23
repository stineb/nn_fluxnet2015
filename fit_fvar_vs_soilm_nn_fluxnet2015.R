rm(list=ls(all=TRUE))

library( dplyr )
library( tidyr )
library( minpack.lm )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
# do.sites   = "FR-Pue"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
nrep       = 5
dotrain    = FALSE
overwrite_modis = FALSE
overwrite_mte = FALSE
##---------------------------------

siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

##------------------------------------------------
## Load overview L1 and initialise additional columns created here
##------------------------------------------------
load( "data/overview_data_fluxnet2015_L1.Rdata" )

## fLUE in the lower quartile [bysite]
overview$fvar_min <- rep( NA, nrow(overview) )

## fLUE in the upper 10% quantile [bysite]
overview$fvar_max <- rep( NA, nrow(overview) )


## initialise matrix of WUE medians by soil moisture-quantile
nintervals <- 5
fvar_vs_soilm <- matrix( NA, length(do.sites), nintervals )

fitparams <- data.frame()

##------------------------------------------------
## Define functions for the fitting and visualisation
##------------------------------------------------
stress_exp <- function( x, ymin, xmin, par_nonlin ){
  outstress <- (1.0 - ymin) * (1.0 - exp( -par_nonlin * (x - xmin) ) ) + ymin
  return( outstress )
}

# stress_quot <- function( x, ymin, xmin, par_nonlin ){
#   outstress <- (1.0 - ymin) * (x - xmin) / (par_nonlin + (x - xmin)) + ymin
#   return( outstress )
# }

stress_quot <- function( x, xmin, xmax, par_nonlin ){
  if ( x > xmax ){
    outstress <- 1.0
  } else {
    outstress <- ( ( x - xmin ) / ( xmax - xmin ) ) ^ par_nonlin
  }
  return( outstress )
}

stress_quad <- function( x, x0, off, apar ){
  outstress <- 1.0 + off - apar * ( x - x0 ) ^ 2
  return( outstress )
}

mycurve <- function( func, from, to, col='black', add=FALSE, lwd=1, lty=1 ){
  range_x <- seq( from, to, by=(to-from)/100 )
  range_y <- sapply( range_x, func )
  if (add){
    lines( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
  } else {
    plot( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
  }
}

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

print( "Fitting functional relationship for all sites ..." )

jdx <- 0
for (sitename in do.sites){
  jdx <- jdx + 1

  infil     <- paste( "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

    ## For SD-Dem, overwrite measured soil moisture - it's probably wrong as it does not fall below 0.25 or so
    if (sitename == "SD-Dem"){
      varnams_swc_mod <- varnams_swc[ !is.element( varnams_swc, varnams_swc_obs ) ]
      nice$soilm_mean <- apply( dplyr::select( nice, one_of(varnams_swc_mod)), 1, FUN=mean, na.rm=TRUE )
      nice$soilm_mean[ is.nan( nice$soilm_mean ) ] <- NA
    }

    ## add row to aggregated data
    mysitename <- data.frame( mysitename=rep( sitename, nrow(nice) ) )

  ##------------------------------------------------
  ## Get functional form of fVAR vs. soil moisture, treating (soilm_xxx, fvar) pairs as individual observations
  ##------------------------------------------------
    ## Get median by interval and get fvar_vs_soilm for this site (used for clustering)
    intervals <- seq( 0.0, 1.0, 1.0/nintervals )
    mid <- intervals[1:(length(intervals)-1)] + (intervals[2]-intervals[1])/2
    nice$ininterval <- NULL
    nice <- nice %>% mutate( ininterval = cut( soilm_mean , breaks = intervals ) ) %>% group_by( ininterval )
    tmp <- nice %>% dplyr::summarise( median=median( fvar, na.rm=TRUE ) ) %>% complete( ininterval, fill = list( median = NA ) ) %>% dplyr::select( median )
    fvar_vs_soilm[ jdx, ] <- unlist(tmp)[1:nintervals]

    ## Get median by interval and get fvar_vs_soilm for this site (used for clustering)
    nsintervals <- 5
    sintervals <- seq( 0.0, 1.0, 1.0/nsintervals )
    xvals <- sintervals[1:(length(sintervals)-1)] + (sintervals[2]-sintervals[1])/2
    nice$insinterval <- NULL
    nice <- nice %>% mutate( insinterval = cut( soilm_mean , breaks = sintervals ) ) %>% group_by( insinterval )
    tmp <- nice %>% dplyr::summarise( median=median( fvar, na.rm=TRUE ) ) %>% complete( insinterval, fill = list( median = NA ) ) %>% dplyr::select( median )
    yvals <- unlist(tmp)[1:nsintervals]

    ## Fit by medians in bis
    df_tmp <- data.frame( xvals=xvals, yvals=yvals )
    gpp_stressfit <- try( 
                          nlsLM( 
                                yvals ~ stress_quad( xvals, x0, off, apar ),
                                data=df_tmp,
                                start=list( x0=1.0, off=0.0, apar=1.0 ),
                                lower=c( 0.01, -1.0, 0.01 ),
                                algorithm="port"
                                ) 
                          )

    ## add fit parameters to aggregated data frame
    addrow <- data.frame( x0=coef(gpp_stressfit)[[ "x0" ]], off=coef(gpp_stressfit)[[ "off" ]], apar=coef(gpp_stressfit)[[ "apar" ]] )
    fitparams <- rbind( fitparams, addrow )

    if (class(gpp_stressfit)!="try-error"){ rmse_gpp_stressfit <- mean(summary(gpp_stressfit)$residuals^2) }

    ##------------------------------------------------
    ## Plot fvar vs. soil moisture
    ##------------------------------------------------
    pdf( paste( "fig_nn_fluxnet2015/fvar_vs_soilm/fvar_vs_soilm_mean_", sitename, ".pdf", sep="" ), width=6, height=5 )

      par(las=1)
      plot( nice$soilm_mean, nice$fvar, xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("royalblue3", 0.2) )
      bp <- boxplot( fvar ~ insinterval, data=nice, main=sitename, col=NA, las=1, outline = FALSE, na.rm=TRUE, add=TRUE, at=(sintervals[1:nsintervals]+(1.0/(2*nsintervals))), boxwex=0.05, axes=FALSE )
      abline( h=1.0, lwd=0.5 )
      points( xvals, yvals, pch=16, col='red' )

      if (class(gpp_stressfit)!="try-error"){ 
        mycurve( function(x) stress_quad( x, coef(gpp_stressfit)[[ "x0" ]], coef(gpp_stressfit)[[ "off" ]], coef(gpp_stressfit)[[ "apar" ]] ), from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )

        ## alternative tested:
        # mycurve( function(x) stress_quot( x, coef(gpp_stressfit)[[ "xmin" ]], coef(gpp_stressfit)[[ "xmax" ]], coef(gpp_stressfit)[[ "par_nonlin" ]] ), from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
        # abline( v=coef(gpp_stressfit)[[ "xmin" ]], lty=3 )
        # abline( v=coef(gpp_stressfit)[[ "xmax" ]], lty=3 )
        # mtext( bquote( "RMSE" == .(format( rmse_gpp_stressfit, digits=3 ) ) ), side=3, line=0, adj=1 )
      }

    dev.off()


    ##------------------------------------------------
    ## Plot the same again but now not into PDFs so that it will be included in knitted file
    ##------------------------------------------------
    if ( sitename %in% c( "AU-DaP", "FR-Pue", "IT-PT1") ){

      par(las=1)
      plot( nice$soilm_mean, nice$fvar, xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("royalblue3", 0.2) )
      bp <- boxplot( fvar ~ insinterval, data=nice, main=sitename, col=NA, las=1, outline = FALSE, na.rm=TRUE, add=TRUE, at=(sintervals[1:nsintervals]+(1.0/(2*nsintervals))), boxwex=0.05, axes=FALSE )
      abline( h=1.0, lwd=0.5 )
      points( xvals, yvals, pch=16, col='red' )

      if (class(gpp_stressfit)!="try-error"){ 
        mycurve( function(x) stress_quad( x, coef(gpp_stressfit)[[ "x0" ]], coef(gpp_stressfit)[[ "off" ]], coef(gpp_stressfit)[[ "apar" ]] ), from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
      }

    }

    
    ##------------------------------------------------
    ## Add to overview data table: fvar in lowest soil moisture quartile
    ##------------------------------------------------
    overview$fvar_min[ which(overview$mysitename==sitename) ] <- fvar_vs_soilm[ jdx, 1 ]

    ##------------------------------------------------
    ## Add to overview data table: fvar in upper 10% soil moisture quantile
    ##------------------------------------------------
    ## Get median by interval and get fvar_vs_soilm for this site (used for clustering)
    nsintervals <- 10
    sintervals <- seq( 0.0, 1.0, 1.0/nsintervals )
    xvals <- sintervals[1:(length(sintervals)-1)] + (sintervals[2]-sintervals[1])/2
    nice <- nice %>% mutate( inssinterval = cut( soilm_mean , breaks = sintervals ) ) %>% group_by( inssinterval )
    tmp <- nice %>% dplyr::summarise( median=median( fvar, na.rm=TRUE ) ) %>% complete( inssinterval, fill = list( median = NA ) ) %>% dplyr::select( median )
    yvals <- unlist(tmp)[1:nsintervals]

    overview$fvar_max[ which(overview$mysitename==sitename) ] <- yvals[nsintervals]

}

print( "... done." )

if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  ## aggregated response of fLUE to soil moisture
  rownames(fvar_vs_soilm) <- do.sites
  save( fvar_vs_soilm, file="data/fvar_vs_soilm.Rdata")

  rownames(fitparams) <- do.sites
  save( fitparams, file="data/fitparams.Rdata")

  save( overview, file="data/overview_data_fluxnet2015_L2.Rdata" )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

