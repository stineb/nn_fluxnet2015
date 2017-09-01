library(abind)
library(dplyr)
library(tidyr)
library(broom)
library(cluster)
library(ggplot2)

source( "add_alpha.R")

load( file="data/overview_data_fluxnet2015_L3.Rdata" )
load( file="data/fvar_vs_soilm.Rdata" )
load( file="data/fitparams.Rdata")

siteinfo <- dplyr::select( overview, mysitename, lon, lat, year_start, year_end, classid )

##------------------------------------------------
## 1. Get sites not affected by very low soil moisture (group 4)
##------------------------------------------------
rowid <- apply( fvar_vs_soilm[,1:2], 1, function(x) sum(!is.na(x)) ) < 2
print( paste( "Total number of sites in group 4 (no low soil moisture):", sum(rowid) ) )
tmp <- data.frame( mysitename=rownames(fvar_vs_soilm), finalcluster=rep(NA, nrow(fvar_vs_soilm)) )
tmp$finalcluster[ rowid ] <- 4

## add to overview
overview <- overview %>% left_join( tmp, by="mysitename" )


##------------------------------------------------
## 2. Get data frame containing fit parameters
##------------------------------------------------
userows <- apply( fvar_vs_soilm, 1, function(x) sum(!is.na(x)) ) == 5
# print( paste( "Total number of sites with data in lowest soil moisture bin:", sum(userows) ) )

## add y-value of quadratic function at x=0 and x=1
stress_quad <- function( x, x0, off, apar ){
  outstress <- 1.0 + off - apar * ( x - x0 ) ^ 2
  return( outstress )
}

fitparams$y_x0 <- rep( NA, nrow(fitparams) )
fitparams$y_x1 <- rep( NA, nrow(fitparams) )
for (idx in 1:nrow(fitparams)){
  fitparams$y_x0[idx] <- stress_quad( 0, fitparams$x0[idx], fitparams$off[idx], fitparams$apar[idx] )
  fitparams$y_x1[idx] <- stress_quad( 1, fitparams$x0[idx], fitparams$off[idx], fitparams$apar[idx] )
}

## replace y_x0 and y_x1 values >1 with 1
fitparams$y_x0[ which(fitparams$y_x0>1) ] <- 1.0
fitparams$y_x1[ which(fitparams$y_x1>1) ] <- 1.0

##------------------------------------------------
## 3. Do clustering - BY HAND
##------------------------------------------------
## do clustering by hand: quadfitcluster = 1 for sites with low sensitivity
tmp <- fitparams[ which(userows), ]
tmp$quadfitcluster <- rep( NA, nrow(tmp) )
tmp$quadfitcluster <- with( tmp, ifelse( y_x0 > 0.8, quadfitcluster <- 1, quadfitcluster <- 2 ) )
tmp$mysitename <- rownames(tmp)

## plot histogram of fLUE_0 (y_x0) values 
out <- hist( tmp$y_x0, breaks=50, plot=FALSE )
par( las=1 )
hist( dplyr::filter( tmp, quadfitcluster==1 )$y_x0, breaks=out$breaks, col=add_alpha("springgreen3", 0.5), xlim=c(0,1.1), ylim=c(0,10), xlab=expression(paste("fLUE"[0])), main="" )
hist( dplyr::filter( tmp, quadfitcluster==2 )$y_x0, breaks=out$breaks, col=add_alpha("grey30", 0.5), add=TRUE )
abline( v=0.8, col="red" )

## add to overview table
overview <- overview %>% left_join( dplyr::select( tmp, mysitename, y_x0, quadfitcluster ), by="mysitename" )
overview$finalcluster[ which( overview$quadfitcluster == 1 ) ] <- 3

# ###------------------------------------------------
# ## 3. Do clustering - K-MEANS
# ##------------------------------------------------
# nclust_best <- 4
# set.seed(2)
# outkmeans <- kmeans( fitparams$y_x0[ which(userows) ], nclust_best )
# tmp <- fitparams[ which(userows), ]
# tmp$quadfitcluster <- outkmeans$cluster
# tmp$mysitename <- rownames(tmp)
# # tmp <- augment( outkmeans, fitparams$y_x0[ which(userows) ] ) %>% rename( mysitename=.rownames, quadfitcluster=.cluster  )

# ## plot histogram of fLUE_0 (y_x0) values 
# out <- hist( tmp$y_x0, breaks=50, plot=FALSE )
# par( las=1 )
# hist( dplyr::filter( tmp, quadfitcluster==which.max(outkmeans$centers) )$y_x0,    breaks=out$breaks, col=add_alpha("springgreen3", 0.5), xlim=c(0,1.1), ylim=c(0,10), xlab=expression(paste("fLUE"[0])), main="" )
# hist( dplyr::filter( tmp, quadfitcluster!=which.max(outkmeans$centers) )$y_x0, breaks=out$breaks, col=add_alpha("grey30", 0.5), add=TRUE )

# ## add to overview table
# overview <- overview %>% left_join( dplyr::select( tmp, mysitename, y_x0, quadfitcluster ), by="mysitename" )
# overview$finalcluster[ which( overview$quadfitcluster == which.max(outkmeans$centers) ) ] <- 3

##------------------------------------------------
## Save all the data with cluster information
##------------------------------------------------
print( paste( "Total number of sites in group 3 (low fLUE sensitivity):", sum(overview$finalcluster==3, na.rm = TRUE) ) )
save( overview, file="data/overview_data_fluxnet2015_L4.Rdata" )
