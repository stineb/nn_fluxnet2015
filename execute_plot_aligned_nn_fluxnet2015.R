source( "/plot_aligned_nn_fluxnet2015.R" )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
siteinfo <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( siteinfo, successcode==1 )$mysitename

## Manual settings ----------------
# do.sites   = c("FR-Pue", "IT-Cpz", "IT-Ro1")
# do.sites   = "FR-Pue"
##--------------------------

print("create aligned plots for all sites ...")

for (sitename in do.sites){

  # print(paste("site:", sitename))
  plot_aligned_nn_fluxnet2015( sitename, nam_target="lue_obs_evi", use_fapar=FALSE, use_weights=FALSE, makepdf=TRUE, verbose=FALSE )

}

print("... done.")
