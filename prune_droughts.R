prune_droughts <- function( instances, is_drought_byvar, fvar_smooth, fvar_smooth_filled, mod_pot=NULL, mod_act=NULL, obs=NULL, soilm=NULL, soilm_mod=NULL, apar=NULL, rmse_cutoff=NA ){

  # ## xxx debug ---------------------------
  # instances      = droughts
  # is_drought_byvar = nice$is_drought_byvar
  # fvar_smooth    = nice$fvar_smooth
  # fvar_smooth_filled = nice$fvar_smooth_filled
  # mod_pot        = nice$var_nn_pot
  # mod_act        = nice$var_nn_act 
  # obs            = nice[[ nam_target ]]
  # soilm     = soilm
  # soilm_mod = soilm_mod
  # apar      = nice$ppfd * nice[[ fapar_data ]]
  # ##---------------------------------------

  if (nrow(instances)>0){

    print(paste("number of events before pruning:", nrow(instances)))

    ##-------------------------------------------------------------------------------
    ## Trim drought event based on missing values in fvar_smooth
    ##-------------------------------------------------------------------------------
    idx_drop <- c()
    for ( idx in 1:nrow(instances) ){

      ## Case 1: NAs at the end of the drought period. Find first non-NA starting from end of drought. Define this as the new end of drought
      if ( is.na( fvar_smooth[ instances$idx_start[idx]+instances$len[idx]-1 ] ) ){
        jdx_look <- instances$idx_start[idx]+instances$len[idx]-2
        while ( is.na( fvar_smooth[ jdx_look ] ) ) jdx_look <- jdx_look - 1
        instances$len[idx] <- jdx_look - instances$idx_start[idx] + 1
      }

      ## Case 2: NAs a the beginning of the drought period (identified based on the interpolated, gap-filled fvar)
      ## ==> drop the event
      if ( is.na( fvar_smooth[ instances$idx_start[idx] ] ) ) idx_drop <- c( idx_drop, idx )

    }

    if (length(idx_drop)>0) instances <- instances[ -idx_drop, ]

    ## keep only droughts longer than 5 days after this trimming
    instances <- instances[ which( instances$len >= 5 ), ]

    print(paste("number of events after trimming:", nrow(instances)))


    ##-------------------------------------------------------------------------------
    ## Prune drought events by whether soil moisture is (relatively) low
    ##-------------------------------------------------------------------------------
    idx_drop <- c()
    if (!is.null(soilm) && nrow(instances)>0 ){

      soilm_threshold <- quantile( soilm, probs=0.75, na.rm=TRUE )
      # soilm_threshold <- 0.5

      idx_drop <- c()

      for ( idx in 1:nrow(instances) ){

        soilm_sub <- soilm[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        if ( !any( !is.na(soilm_sub) ) ) soilm_sub <- soilm_mod[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        if ( mean( soilm_sub, na.rm=TRUE ) > soilm_threshold ) idx_drop <- c( idx_drop, idx )

      }

      if (length(idx_drop)>0) instances <- instances[ -idx_drop, ]
      print(paste("number of events after pruning by soil moisture:", nrow(instances)))

    }


    ##-------------------------------------------------------------------------------
    ## Prune drought events by quality of modelled vs. observed during drought period
    ##-------------------------------------------------------------------------------
    idx_drop <- c()
    if (!is.null(mod_act) && !is.null(mod_pot) && !is.null(obs)  && nrow(instances)>0 ){

      require( hydroGOF )
      idx_drop <- c()

      for ( idx in 1:nrow(instances) ){

        ## extract data for this event
        mod_act_sub <- mod_act[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        mod_pot_sub <- mod_pot[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        obs_sub     <- obs[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        
        ## remove NAs
        tmp <- mod_act_sub + mod_pot_sub + obs_sub
        mod_act_sub <- mod_act_sub[ which(!is.na(tmp)) ]
        mod_pot_sub <- mod_pot_sub[ which(!is.na(tmp)) ]
        obs_sub <- obs_sub[ which(!is.na(tmp)) ]

        ## get RMSE
        # rmse_act_sub    <- rmse( obs_sub, mod_act_sub )
        # rmse_pot    <- rmse( obs_sub, mod_pot_sub )
        rmse_act_sub    <- rmse( mod_act_sub, obs_sub )
        rmse_pot_sub    <- rmse( mod_pot_sub, obs_sub )

        rmse_act    <- rmse( mod_act, obs )
        rmse_pot    <- rmse( mod_pot, obs )

        # ## drop event (instance) if RMSE of NN-pot is smaller than RMSE of NN-act
        # print(paste("rmse_pot_sub", rmse_pot_sub))
        # print(paste("rmse_act_sub", rmse_act_sub))

        # boxplot( mod_act_sub - obs_sub, at=1, xlim=c(0,3) )
        # boxplot( mod_pot_sub - obs_sub, at=2, xlim=c(0,3), add=TRUE )

        if ( rmse_pot_sub < 0.9 * rmse_act_sub ) idx_drop <- c( idx_drop, idx )

        # # print(paste("RMSE during event", idx, "=", rmse))
        # # print("--------------------------------")
        # if (!is.na(rmse_cutoff)){
        #   if (is.nan(rmse) || rmse>rmse_cutoff) { idx_drop <- c( idx_drop, idx ) }
        # }

      }

      if (length(idx_drop)>0) instances <- instances[ -idx_drop, ]
      print(paste("number of events after pruning by mod vs obs:", nrow(instances)))

    }

    ##-------------------------------------------------------------------------------
    ## Prune drought events if observed GPP is very low (due to low light)
    ##-------------------------------------------------------------------------------
    idx_drop <- c()
    if ( !is.null(apar) && nrow(instances)>0  && nrow(instances)>0 ){

      ## multiplying NN-modelled potential LUE with APAR
      gpp <- obs * apar

      for ( idx in 1:nrow(instances) ){

        ## extract data for this event
        gpp_sub <- gpp[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ]
        
        if ( mean( gpp_sub, na.rm=TRUE ) < 0.05 * max( gpp, na.rm=TRUE ) ){
          ## insignificant "drought" - probably just noise under low light levels
          idx_drop <- c( idx_drop, idx )
        }

      }

      if (length(idx_drop)>0) instances <- instances[ -idx_drop, ]
      print(paste("number of events after pruning by GPP level:", nrow(instances)))

    }

    ##-------------------------------------------------------------------------------
    ## Update 'is_drought_byvar'
    ##-------------------------------------------------------------------------------
    print(paste("number of events after pruning:", nrow(instances)))

    is_drought_byvar[] <- FALSE
    if (nrow(instances)>0){
      for ( idx in 1:nrow(instances) ){
        is_drought_byvar[ instances$idx_start[idx]:(instances$idx_start[idx]+instances$len[idx]-1) ] <- TRUE
      }
    }

    # test_instances <- get_consecutive( is_drought_byvar )
    # print("original:")
    # print( instances )
    # print("recalculated from is_drought_byvar:")
    # print(test_instances)

  }

  out <- list( instances=instances, is_drought_byvar=is_drought_byvar )
  return( out )

}
