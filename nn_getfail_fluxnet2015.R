.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

source( paste( "function_nn_getfail_fluxnet2015.R", sep="" ) )

siteinfo <- read.csv( "soilm_data_usability_fluxnet2015.csv", as.is=TRUE )
nam_target  = "lue_obs_evi"
use_weights = FALSE
use_fapar   = FALSE

print( "getting failure info for all sites ..." )
successcodes <- data.frame()
for (idx in 1:nrow(siteinfo)){
  sitename <- siteinfo$mysitename[idx]
  add <- nn_getfail_fluxnet( sitename, code=siteinfo$code[idx], nam_target=nam_target, use_weights=use_weights, use_fapar=use_fapar )
  successcodes <- rbind( successcodes, add )
}
print( "... done" )

# Force exclude by hand - no good criterium found
force_exclude2 <- c( "CN-Din", "CZ-BK1", "IT-CA1", "IT-CA2", "IT-CA3", "IT-Ro2", "US-Me6", "US-Tw1" )
print( "Force exclusion of the following sites: " )
print( force_exclude2 )
successcodes$successcode[ which( is.element(successcodes$mysitename, force_exclude2)) ] <- 3

print( paste( "number of sites with code 1:", sum(successcodes$successcode==1) ) )
print( paste( "number of sites with code 2:", sum(successcodes$successcode==2) ) )
print( paste( "number of sites with code 3:", sum(successcodes$successcode==3) ) )
print( paste( "number of sites with code 0:", sum(successcodes$successcode==0) ) )

write.csv( successcodes, file=paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/successcodes.csv", sep="" ), row.names=FALSE )




