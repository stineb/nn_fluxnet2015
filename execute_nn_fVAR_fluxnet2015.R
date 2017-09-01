source( "nn_fVAR_fluxnet2015.R" )

## Set this to true if this is to be executed on local machine (not through R CMD BATCH call by job submitted to cluster)
local_execution <- TRUE

if (local_execution) {

  ## For local execution:
  siteinfo <- read.csv( "soilm_data_usability_fluxnet2015.csv", as.is=TRUE )
  do.sites <- dplyr::filter( siteinfo, code!=0 )$mysitename
  nam_target  = "lue_obs_evi"
  use_weights = FALSE
  use_fapar   = FALSE

} else {

  ##----------------------------------------------------------
  ## Evaluate arguments provided by R CMD BATCH call
  ##----------------------------------------------------------
  ## First read in the arguments listed at the command line.
  ## Call this by 
  ## 'R CMD BATCH --no-save --no-restore '--args sitename="FR-Pue" nam_target="lue_obs"' nn_fVAR_fluxnet2015.R nn_fVAR_fluxnet2015.out &'
  args=(commandArgs(TRUE))

  ## args is now a list of character vectors
  ## First check to see if arguments are passed.
  ## Then cycle through each element of the list and evaluate the expressions.
  if (length(args)==0){
    print("No arguments supplied. Provide at least sitename.")
  } else {
    for (i in 1:length(args)){
       eval( parse( text=args[[i]] ) )
    }
  }
  ##----------------------------------------------------------

}

for ( sitename in do.sites ){
  nn_fVAR_fluxnet( sitename, nam_target=nam_target, use_weights=use_weights, use_fapar=use_fapar, overwrite=TRUE )
}
