##------------------------------------------------
## Get Standardised Precipitation Index and Standardised Precipitation Evapotranspiration Index
## Returns a monthly data frame with monthly total precipitation, PET, (prec-PET), and SPI, SPEI values
##------------------------------------------------
get_spi_bysite <- function( sitename, lat ){

  require( dplyr )
  require( SPEI )

  ## debug ---------------
  # load( "./data/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata" )
  # sitename = "FR-Pue"
  # lat = 43.7414
  ##----------------------

  ## get precipitation
  path <- paste( "./input_fluxnet2015_sofun/sitedata/climate/", sitename, "/clim_daily_byst_", sitename, ".csv", sep="" )
  ddf <- try( read.csv( path ) )

  ## do spi only with watch data
  ddf$prec <- ddf$prec_watch
  ddf$pet  <- thornthwaite( ddf$temp_watch, lat, na.rm=TRUE )

  # ## get PET in mm
  # iset <- "s13"
  # dirnam_mod <- paste( myhome, "sofun/output_fluxnet2015_sofun/", iset, "/", sep="" )
  # tmp <- get_daily_modelout( sitename, dirnam_mod, c("pet","econ") )
  # tmp <- tmp %>% mutate( pet_mm = pet * econ * 1e-9 )

  # ddf <- ddf %>% left_join( dplyr::select( tmp, year, moy, doy, pet, econ, pet_mm ), by=c("year", "doy", "moy") )
  # adf <- ddf %>% group_by( year ) %>% summarise( pet=sum( pet_mm, na.rm=TRUE ), prec=sum( prec_meteo, na.rm=TRUE ) )

  ## aggregate to monthly values and calculate balance
  mdf <- ddf %>%  group_by( year, moy ) %>% 
                  summarise( pet=sum( pet ), prec=sum( prec ) ) %>%
                  mutate( bal=prec-pet, year_dec=year+(moy-1)/12+15/365 )

  idx_na <- which( is.na(mdf$prec) | is.na(mdf$pet) )


  if (length(unique(mdf$year))>5){

    ## add SPI and SPEI based on 1 and 3 months
    mdf$spi1 <- as.vector( spi( mdf$prec, 1, na.rm=TRUE )$fitted )
    mdf$spi3 <- as.vector( spi( mdf$prec, 3, na.rm=TRUE )$fitted )

    mdf$spei1 <- as.vector( spei( mdf$bal, 1, na.rm=TRUE )$fitted )
    mdf$spei3 <- as.vector( spei( mdf$bal, 3, na.rm=TRUE )$fitted )

    mdf$spi1[ idx_na ] <- NA
    mdf$spi3[ idx_na ] <- NA

    mdf$spei1[ idx_na ] <- NA
    mdf$spei3[ idx_na ] <- NA

  } else {

    mdf$spi1 <- rep( NA, nrow(mdf) )
    mdf$spi3 <- rep( NA, nrow(mdf) )

    mdf$spei1 <- rep( NA, nrow(mdf) )
    mdf$spei3 <- rep( NA, nrow(mdf) )

  }

  return( mdf )

}
