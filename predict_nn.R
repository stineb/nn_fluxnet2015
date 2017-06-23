predict_nn <- function( data, predictors, nam_target, weights=NULL, nn=NULL, threshold=0.03, do_predict=TRUE, do_modobs=FALSE, trainfrac=1.0, lifesign="none", package="neuralnet", nrep=1, seed=1, hidden=NULL ){

  # ## xxx debug
  # data <- train
  # predictors <- c("temp", "ppfd", "vpd")
  # nam_target <- "le_f_mds"
  # lifesign <- "full"
  # package <- "neuralnet"
  # nrep <- 1

  if (package=="neuralnet"){

    require( neuralnet )
    source( "repeat_neuralnet.R" )

    ## debug
    # data=data
    # predictors=c( "ppfd", "temp", "vpd", "soilm" )
    # nn=NA
    # nam_target=c("gpp_obs", "et_obs")
    # # nam_target="gpp_obs"
    # do_predict=TRUE
    # do_modobs=TRUE
    
    ## sample training data
    index <- sample( 1:nrow( data ), round( trainfrac * nrow( data ) ) )  
    train <- data[index,]
    if (trainfrac==1.0){
      test <- data[index,]
    } else {
      test  <- data[-index,]
    }

    ## scale variables to within [0,1]
    maxs       <- apply(data, 2, max) 
    mins       <- apply(data, 2, min)
    scaled     <- as.data.frame( scale( data, center = mins, scale = maxs - mins ) )
    train_     <- scaled[index,]
    if (trainfrac==1.0){
      test_ <- scaled[index,]
    } else {
      test_  <- scaled[-index,]
    }

    ## train the NN
    names  <- names(train_)
    forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )
    if (is.null(nn)){
      # print( "getting NN for formula:", paste( paste( nam_target, collapse=" + "), "~", paste( predictors, collapse=" + " ) ) )
      # print( nam_target )
      nn <- neuralnet( forml, data=train_, hidden=6, linear.output=TRUE, lifesign="none", threshold=threshold, rep=3 )
      # nn <- repeat_neuralnet( forml, train_, test_, predictors, hidden=6, linear.output=TRUE, lifesign=lifesign, threshold=threshold, rep=nrep )
    }

    if (do_predict){
      ## predicting
      pred_nn_ <- try( compute( nn, subset( scaled, select=predictors ) ) )

      if (class(pred_nn_)!="try-error"){
        
        ## scaling back
        range   <- max( data[[ nam_target ]] ) - min( data[[ nam_target ]] )
        offset  <- min( data[[ nam_target ]] )
        pred_nn <- pred_nn_$net.result * range + offset

        # if (do_modobs){
        #   ##---------------------------------------------------------
        #   ## Modelled vs. observed, and statistics
        #   ##---------------------------------------------------------
        #   ## plot modelled (NN) vs. observed
        #   if (trainfrac==1.0){
        #     mod <- pred_nn[index]
        #     obs <- data[[ nam_target ]][index]
        #   } else {
        #     mod <- pred_nn[-index]
        #     obs <- data[[ nam_target ]][-index]            
        #   }
        #   out <- analyse_modobs( 
        #     mod, obs,
        #     plot.fil=paste( "fig_nn_fluxnet2015/modobs_", nam_target, "_", paste( predictors, collapse="_" ), "_", sitename, ".pdf", sep="" ), 
        #     plot.xlab=expression(paste("observed GPP (gC m"^{-2}, " d"^{-1}, ")")), 
        #     plot.ylab=expression(paste("modelled GPP (gC m"^{-2}, " d"^{-1}, ")")), 
        #     plot.xylim=c(0,max(obs)), 
        #     plot.title="", plot.col=rgb(0,0,0,0.3)
        #   )
        # }

      } else {
        pred_nn=NA
      }
    } else {
      pred_nn=NA
    }

    return( list( nn=nn, vals=as.vector( pred_nn ), hidden_best=6 ) )
    
  } else if (package=="nnet"){

    ## using 'caret' package

    require( nnet )
    require( caret )

    if (is.null(nn)){

      forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )

      # ## Data partitioning 
      # index <- createDataPartition( y = data$gpp_obs,  ## the outcome data are needed
      #                               p = 0.75,          ## The percentage of data in the training set
      #                               list = FALSE       ## The format of the results
      #                             )
      # training <- data[ index,]
      # testing  <- data[-index,]

      # print("formula:")
      # print(forml)
      # print("--------")

      ## this has caused a problem before due to values being too hight -> weird
      if (nam_target=="le_f_mds"|| nam_target=="et_obs"){ data[[ nam_target ]] <- data[[ nam_target ]] * 1e-6 }

      preprocessParams <- preProcess( data, method=c("range") )
      traincotrlParams <- trainControl( method="repeatedcv", number=5, repeats=5, verboseIter=FALSE, p=0.75 ) # take best of 10 repetitions of training with 75% used for training (25% for testing)
      # traincotrlParams <- trainControl( method="cv", number=10, verboseIter=FALSE ) # 5-fold cross-validation

      if (is.null(hidden)){
        tune_grid <- expand.grid( .decay = c(0.1), .size = seq(4,20,2) )
      } else {
        tune_grid <- expand.grid( .decay = c(0.1), .size = c(hidden) )
      }

      set.seed(seed)
      nn <- train(
                  forml,
                  data      = data, #training,
                  weights   = weights,
                  method    = "nnet",
                  linout    = TRUE,
                  tuneGrid  = tune_grid,
                  preProc   = "range", # preProc  = preprocessParams
                  trControl = traincotrlParams,
                  trace     = FALSE
                  )
      # pdf("fig_nn_fluxnet2015/caret_profile.pdf")
      # plot( nn )
      # dev.off()

    }

    if (do_predict){
      vals <- as.vector( predict( nn, data ) )  # try( predict( nn, newdata=testing ) )
    } else {
      vals <- rep( NA, nrow(data) )
    }
    
    if (nam_target=="le_f_mds"|| nam_target=="et_obs"){ vals <- vals * 1e6 }

  }

  # } else if (package=="nnet"){

  #   require( nnet )

  #   # ## sample training data
  #   # index <- sample( 1:nrow( data ), round( trainfrac * nrow( data ) ) )  
  #   # train <- data[index,]
  #   # if (trainfrac==1.0){
  #   #   test <- data[index,]
  #   # } else {
  #   #   test  <- data[-index,]
  #   # }

  #   # ## scale variables to within [0,1]
  #   # maxs       <- apply( data, 2, max ) 
  #   # mins       <- apply( data, 2, min )
  #   # scaled     <- as.data.frame( scale( data, center = mins, scale = maxs - mins ) )
  #   # train_     <- scaled[index,]
  #   # if (trainfrac==1.0){
  #   #   test_ <- scaled[index,]
  #   # } else {
  #   #   test_  <- scaled[-index,]
  #   # }

  #   # ## train the NN
  #   # names  <- names(train_)
  #   # forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )
  #   # if (is.na(nn)){
  #   #   nn_nnet <- nnet( forml, data = train_, linout = TRUE, size = 20 )
  #   # }

  #   # ## scaling back
  #   # range   <- max( data[[ nam_target ]] ) - min( data[[ nam_target ]] )
  #   # offset  <- min( data[[ nam_target ]] )
  #   # vals    <- as.vector( nn_nnet$fitted.values ) * range + offset    

  #   # analyse_modobs( vals, data[[ nam_target ]] )

  #   ##-------------------------------
  #   ## train the NN
  #   forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )

  #   ## using 'nnet'
  #   nn_nnet <- nnet( forml, data = data, linout = TRUE, size = 20 )

  #   analyse_modobs( as.vector( nn_nnet$fitted.values ), data$gpp_obs )
  #   vals <- as.vector( nn_nnet$fitted.values )

  #   analyse_modobs( vals, data[[ nam_target ]] )

  # }

  return( list( nn=nn, vals=vals, hidden_best=nn$bestTune$size ) )

}
