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
