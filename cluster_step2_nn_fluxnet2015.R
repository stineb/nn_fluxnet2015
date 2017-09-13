library(abind)
library(dplyr)
library(tidyr)
library(broom)
library(cluster)

load( file="data/overview_data_fluxnet2015_L4.Rdata" )
load( file="data/fvar_vs_soilm.Rdata" )

source( "get_aggresponse_binned.R" )

## Use only sites not grouped into cluster 3 or 4 by 'cluster_step1_nn_fluxnet2015' and sites with a significant number of drought days (successcode == 1)
do.sites <- filter( overview, is.na( finalcluster ) & successcode == 1 ) %>% select( mysitename )

## Manual settings ----------------
# do.sites   = "FR-Pue"
nam_target = "lue_obs_evi"
use_fapar  = FALSE
fapar_data = "evi"
package    = "nnet"
nrep       = 5
dotrain    = FALSE
##--------------------------

## check and override if necessary
if ( nam_target=="lue_obs" || nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
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

zerobin <- 1

if (length(do.sites$mysitename)>=10){
  fvar_agg     <- c()
  fapar_agg    <- c()
  # iwue_agg     <- c()
  df_dday_agg  <- c()
  sitename_agg <- c()
} else {
  print("I'm assuming that no clustering is done with so few sites!!!")
}

##------------------------------------------------
## Collect binned response in different variables (only fapar used finally)
##------------------------------------------------
print( paste( "total number of sites left for remaining clusters:", nrow(do.sites) ) )
for (sitename in do.sites$mysitename){

  out <- get_aggresponse_binned( sitename, verbose=FALSE )

  if (!is.null(out$fapar       ))  fapar_agg        <- rbind( fapar_agg,        out$fapar    )
  if (!is.null(out$fvar        ))  fvar_agg         <- rbind( fvar_agg,         out$fvar     )
  if (!is.null(out$sitename    ))  sitename_agg     <- rbind( sitename_agg,     out$sitename )
  if (!is.null(out$fapar0      ))  fapar0_agg       <- rbind( fapar0_agg,       out$fapar0   )
  # if (!is.na(out$iwue    ))  iwue_agg     <- rbind( iwue_agg,     out$iwue    )

}

## add row names (site name) to naked array
rownames(fvar_agg)  <- sitename_agg
rownames(fapar_agg) <- sitename_agg
# rownames(iwue_agg)  <- sitename_agg


##------------------------------------------------
## DO CLUSTER ANALYSIS ON SHAPE OF FAPAR (AND PTENTIALLY OTHERS) DURING DROUGHT
##------------------------------------------------
## Combine variables based on which clustering is done into a single array 'mega'
mega <- cbind( fapar_agg, fvar_agg )    # this leads to ZM-Mon and AU-ASM being classified in drought-deciduous cluster

# ## Use fapar level only in third bin
# mega <- fapar_agg[,3]    # this leads to ZM-Mon and AU-ASM being classified in evergreen cluster

# ## get optimal number of clusters using the Gap Statistics, from http://www.sthda.com/english/wiki/print.php?id=239
# set.seed(1982)
# gap_stat <- clusGap( mega[,2:7], FUN = kmeans, nstart = 2, K.max = 10, B = 500 )
# nclust_best <- with( gap_stat, maxSE( Tab[,"gap"], Tab[,"SE.sim"] ) )  # from http://widequestion.com/question/retrieving-the-optimal-number-of-clusters-in-r/
# print( paste( "optimal number of clusters:", nclust_best ) )
# plot( gap_stat, xlab = "Number of clusters k" )
# abline( v = nclust_best, lty = 2 )

## override number of clusters
nclust_best <- 2

## do clustering again with best number of clusters
set.seed(1982)
outkmeans <- kmeans( mega, nclust_best )

if (is.null(dim(mega))){

  ## Data used for clustering is actually just a vector (one-dimensional) [USED FOR TESTING ALTERNATIVES]
  df_cluster <- data.frame( mysitename=names( outkmeans$cluster ), alignedcluster=outkmeans$cluster, evi_norm3=mega )

  #  ## Complement for visualisation
  # mega_vis <- cbind( fvar_agg, fapar_agg )
  # mega_vis <- as.data.frame( mega_vis )
  # mega_vis$mysitename <- rownames( mega_vis )
  # mega_vis <- mega_vis %>% left_join( df_cluster )

} else {
 
  df_cluster <- augment( outkmeans, mega )
  df_cluster <- df_cluster %>% rename( mysitename=.rownames, alignedcluster=.cluster ) # %>% rename( dfapar=unrowname.x. )
  # mega_vis <- augment( outkmeans, mega_vis )
  # mega_vis <- mega_vis %>% rename( mysitename=.rownames, alignedcluster=.cluster ) # %>% rename( dfapar=unrowname.x. )
 
}

##------------------------------------------------
## Visualisation of within-SS
##------------------------------------------------
print( "The following figure shows the dependence of the total within-cluster sum of squared differences from the cluster mean:")
## For visualisation of within-cluster SS vs. number of clusters, do K-means clustering for a set of numbers of clusters (1:9)
kclusts <- data.frame( k=1:9 ) %>% group_by(k) %>% do( kclust=kmeans( as.array(mega), .$k ) )

## use library "broom" to get statistics of kmeans clustering, from https://cran.r-project.org/web/packages/broom/vignettes/kmeans.html
# clusters    <- kclusts %>% group_by(k) %>% do(tidy(.$kclust[[1]]))
# assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], mega ) )
clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))

## plot total within-SS vs. number of clusters
library( ggplot2 )
plt <- ggplot( clusterings, aes(k, tot.withinss) ) + geom_line()
print( plt )

# ##------------------------------------------------
# ## Visualisation of cluster response
# ##------------------------------------------------
# print("The following figure shows a visualisation of the mean cluster response (centers):")
# ## see clusters center fVAR response across soil moisture intervals
# plot( outkmeans$centers[1,1:3], type='l', ylim=c(0,1.2), col='springgreen4' )
# points( outkmeans$centers[1,1:3], col='springgreen4', pch=16 )

# lines( outkmeans$centers[1,4:6], col='tomato' )
# points( outkmeans$centers[1,4:6], col='tomato', pch=16 )

# lines( outkmeans$centers[2,1:3], col='springgreen4', lty=2 )
# points( outkmeans$centers[2,1:3], col='springgreen4', pch=16 )

# lines( outkmeans$centers[2,4:6], col='tomato', lty=2 )
# points( outkmeans$centers[2,4:6], col='tomato', pch=16 )

# legend( "bottomleft", c("cluster 1", "cluster 2"), lty=c(1,2), bty="n" )
# legend( "bottom", c("fLUE", expression( paste("EVI / EVI"[0]))), col=c('tomato', 'springgreen4'), lty=1, bty="n" )

##------------------------------------------------
## Histogram of greenness change (median of third bin, representing days 21-40 into fLUE droughts).
##------------------------------------------------
print("The following figure shows a histogram of the greenness changes (median of third bin, representing days 21-40 into fLUE droughts) by cluster 1 and 2:")
## plot histogram of fLUE_0 (y_x0) values 
out <- hist( df_cluster$evi_norm3, plot=FALSE )
par( las=1 )
hist( dplyr::filter( df_cluster, alignedcluster==1 )$evi_norm3, breaks=out$breaks, col=add_alpha("tomato", 0.5), xlab=expression(paste("1 - ", Delta, "fAPAR")), main="", ylim=c(0,max(out$counts)) )
hist( dplyr::filter( df_cluster, alignedcluster==2 )$evi_norm3, breaks=out$breaks, col=add_alpha("royalblue2", 0.5), add=TRUE )


##------------------------------------------------
## Complement overview table
##------------------------------------------------
overview <- overview %>% left_join( df_cluster, by="mysitename" )
overview$finalcluster[ which(is.na(overview$finalcluster)) ] <- overview$alignedcluster[ which(is.na(overview$finalcluster)) ]

##------------------------------------------------
## Save all the data with cluster information
##------------------------------------------------
# save( siteinfo_sub, file="siteinfo_alignedcluster.Rdata" )
save( overview, file="data/overview_data_fluxnet2015_L5.Rdata" )


