## save fLUE data for public distribution
library( dplyr )
library( lubridate )
library( readr )

load( "data/nice_agg_lue_obs_evi.Rdata" )
load( "data/overview_data_fluxnet2015_L5.Rdata" )

out <- nice_agg %>% select( mysitename, year, doy, fvar, is_drought_byvar ) %>%
                    mutate( date=ymd( paste0( as.character(year), "-01-01" ) ) + days( doy - 1 ) ) %>%
                    left_join( select( overview, mysitename, finalcluster ) ) %>%
                    mutate( cluster=ifelse( finalcluster==1, "cDD", 
                                      ifelse( finalcluster==2, "cGR", 
                                        ifelse( finalcluster==3, "cLS", 
                                          ifelse( finalcluster==4, "cNA", NA ) 
                                          ) 
                                        ) 
                                      ) 
                                    ) %>%
                    select( site=mysitename, date, year, doy, flue=fvar, is_flue_drought=is_drought_byvar, cluster )

write_csv( out, path="data/flue_stocker18nphyt.csv" )


# Description:

# This dataset contains fLUE data as described in Stocker et al., 2018 (New Phytologist). fLUE is derived from the FLUXNET 2015 dataset, Tier 1, daily. Only sites are included where the method for quantifying fLUE satisfied performance criteria (see Stocker et al. 2018). 
# site: Site name (ID) from the FLUXNET network
# date: DD/MM/YY
# year: year
# doy: day of year
# fLUE: unitless, fraction of actual over potential light use efficiency, derived from artificial neural networks. This quantifies the fractional reduction in light use efficiency due to soil moisture (1 = no reduction).
# is_flue_drought: TRUE if the data is identified as a 'drought' based on deviation of fLUE from 1 (see Stocker et al., 2018)
# cluster: sites are assigned to clusters based on their typical parallel evolution of greenness and fLUE throughout drought events. cDD: 'drought deciduous', cGR: 'evergreen', cLS: 'low sensitivity', cNA: 'not affected'. 