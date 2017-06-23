get_iwue <- function( et_obs, gpp_obs, vpd, prec ){
  ##--------------------------------------------------------------------
  ## returns inherent water use efficiency after Beer et al. (2009)
  ##
  ## arguments: 
  ## et_obs ( J m-2 d-1 )
  ## gpp_obs ( gC m-2 d-1 )
  ## vpd ( Pa )
  ## prec (mm d-1, required to determine whether it's a dry day for assumption E=0)
  ## 
  ## return:
  ## IWUE* (1e-3 * Pa g J-1)
  ##--------------------------------------------------------------------
  source( paste( myhome, "/get_consecutive.R", sep="" ) )

  ## inherent water use efficiency after Beer et al., 2009
  iwue <- gpp_obs * vpd / et_obs * 1e3

  ## Prune days with where E=0 does not hold
  ## indentify instances of dry days row
  ## dry is when precip = 0 and number of days since last day with prec > 0 is > 1. 
  prec_threshold <- 0
  tmp <- ( prec <= prec_threshold )
  instances <- get_consecutive( tmp, do_merge=FALSE, leng_threshold=3 )

  dry <- rep( FALSE, length(et_obs) )
  for (idx in 1:nrow(instances)){
    dry[(instances$idx_start[idx]+2):(instances$idx_start[idx]+instances$len[idx]-1)] <- TRUE
  }

  iwue[which(!dry)] <- NA

  return( iwue )

}
