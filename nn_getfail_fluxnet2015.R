source( paste( "/function_nn_getfail_fluxnet2015.R", sep="" ) )

siteinfo <- read.csv( paste( "/soilm_data_usability_fluxnet2015.csv", sep="" ), as.is=TRUE )
nam_target  = "lue_obs_evi"
use_weights = FALSE
use_fapar   = FALSE

print( "getting failure info for all sites ..." )
successcodes <- data.frame()
for (idx in 1:nrow(siteinfo)){
  sitename <- siteinfo$mysitename[idx]
  add <- nn_getfail_fluxnet( sitename, nam_target=nam_target, use_weights=use_weights, use_fapar=use_fapar )
  successcodes <- rbind( successcodes, add )
}
print( "... done" )

# Force exclude by hand - no good criterium found
force_exclude2 <- c( "AU-Cpr", "AU-GWW", "AU-Tum", "CN-Dan", "CN-Din", "DE-SfN", "DE-SfN", "ES-LgS", "IT-CA1", "IT-CA2", "IT-CA3", "IT-Ro2", "US-ARb", "US-ARc", "US-Me6" )
print( "Force exclusion of the following sites: " )
print( force_exclude2 )
successcodes$successcode[ which( is.element(successcodes$mysitename, force_exclude2)) ] <- 3

print( paste( "number of sites with code 1:", sum(successcodes$successcode==1) ) )
print( paste( "number of sites with code 2:", sum(successcodes$successcode==2) ) )
print( paste( "number of sites with code 3:", sum(successcodes$successcode==3) ) )
print( paste( "number of sites with code 0:", sum(successcodes$successcode==0) ) )

write.csv( successcodes, file=paste( "successcodes.csv", sep="" ), row.names=FALSE )




