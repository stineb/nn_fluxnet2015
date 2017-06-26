source( paste( "plot_bysite_nn_fluxnet2015.R", sep="" ) )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
nrep       = 5
dotrain    = FALSE
overwrite_modis = FALSE
overwrite_mte = FALSE
##---------------------------------

print("creating time series plots for all sites ...")
for (sitename in do.sites){

  plot_bysite_nn_fluxnet2015( sitename, nam_target=nam_target, use_fapar=use_fapar, use_weights=use_weights, makepdf=TRUE )

}
print("... done.")
