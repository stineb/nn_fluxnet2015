predict_nn <- function( data, predictors, nam_target, weights=NULL, nn=NULL, threshold=0.03, do_predict=TRUE, do_modobs=FALSE, trainfrac=1.0, lifesign="none", package="neuralnet", nrep=1, seed=1, hidden=NULL ){

  # ## xxx debug
  # data <- train
  # predictors <- c("temp", "ppfd", "vpd")
  # nam_target <- "le_f_mds"
  # lifesign <- "full"
  # package <- "neuralnet"
  # nrep <- 1

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


  return( list( nn=nn, vals=vals, hidden_best=nn$bestTune$size ) )

}
