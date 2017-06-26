library( dplyr )

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
## Initialise overview table
##------------------------------------------------
overview <- siteinfo

## Total fractional GPP loss due to soil moisture
overview$fgpp_total <- rep( NA, nrow(overview) )
overview$fgpp_drought <- rep( NA, nrow(overview) )  # during drought days only

## Fraction of all days classified as drought (fLUE below threshold)
overview$fdroughtdays <- rep( NA, nrow(overview) )

## perr_XXX: percentage error of P-model, MODIS and MTE during drouht days
overview$perr_pmodel <- rep( NA, nrow(overview) )


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


print( "Getting data for overview table for all sites ..." )

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  infil <- paste( myhome, "/data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
  load( infil ) ## gets list 'nn_fluxnet'
  nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
  
  idxs_drought <- which( nice$is_drought_byvar )
  overview[ which(overview$mysitename==sitename), ]$perr_pmodel <- 100 * sum( nice$gpp_pmodel[idxs_drought] - nice$gpp_obs[idxs_drought], na.rm=TRUE ) / sum( nice$gpp_obs[idxs_drought], na.rm=TRUE )

  if (nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar"){
    nice <- nice %>% mutate( gpp_nn_act = var_nn_act * evi * ppfd, gpp_nn_pot = var_nn_pot * evi * ppfd )
  } else {
    nice <- nice %>% mutate( gpp_nn_act = var_nn_act, gpp_nn_pot = var_nn_pot )
  }

  ##----------------------------------------------------------------------------------------
  ## Data for overview data table
  ##----------------------------------------------------------------------------------------
  idxs <- which( !is.na(nice$var_nn_act) & !is.na(nice$var_nn_pot) )
  overview$fgpp_total[ which(overview$mysitename==sitename) ] <- 100 * ( 1.0 - ( sum( nice$gpp_nn_act[idxs] , na.rm=TRUE ) / sum( nice$gpp_nn_pot[idxs], na.rm=TRUE )))

  idxs <- which( !is.na(nice$var_nn_act) & !is.na(nice$var_nn_pot) & nice$is_drought_byvar==1 )
  overview$fgpp_drought[ which(overview$mysitename==sitename) ] <- 100 * (1.0 - sum( nice$gpp_nn_act[idxs], na.rm=TRUE ) / sum( nice$gpp_nn_pot[idxs], na.rm=TRUE ))

  idxs_drought <- which( !is.na(nice$var_nn_act) & !is.na(nice$var_nn_pot) & nice$is_drought_byvar==1 )
  idxs_nondrought <- which( !is.na(nice$var_nn_act) & !is.na(nice$var_nn_pot) )
  overview$fdroughtdays[ which(overview$mysitename==sitename) ] <- 100 * ( length(idxs_drought) / length(idxs_nondrought))

}

print("... done.")

if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## Get data from 'successcodes' (see 'nn_getfail_fluxnet2015.R') and merge into overview dataframe
  ##------------------------------------------------
  overview <- overview %>% left_join( dplyr::select( successcodes, mysitename, NNall_rsq, NNgood_rsq, successcode ), by="mysitename" )

  ##------------------------------------------------
  ## Save collected data
  ##------------------------------------------------
  filn <- "data/overview_data_fluxnet2015_L1.Rdata"
  print( paste( "Saving to file", filn ) )
  save( overview, file=filn )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

