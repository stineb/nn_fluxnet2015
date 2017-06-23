cutna_headtail <- function( vec ){

  ## Remove (cut) NAs from the head and tail of a vector.
  ## Returns the indexes to be dropped from a vector

  ## remove NAs from head
  if (is.na(vec[1])){
    idx <- 0
    while ( idx < length(vec) ){
      idx <- idx + 1
      test <- head( vec, idx )
      if (any(!is.na(test))){
        ## first non-NA found at position idx
        cuthead <- idx - 1
        break
      }
    }
    # vec <- vec[ -(1:cuthead) ]
    idxs_head <- 1:cuthead
  } else {
    idxs_head <- c()
  }


  ## remove NAs from tail
  if (is.na(vec[length(vec)])){
    idx <- 0
    while ( idx < length(vec) ){
      idx <- idx + 1
      test <- tail( vec, idx )
      if (any(!is.na(test))){
        ## first non-NA found at position idx, counting from tail
        cuttail <- idx - 1
        break
      }
    }
    # vec <- vec[ -((length(vec)-cuttail+1):length(vec)) ]
    idxs_tail <- (length(vec)-cuttail+1):length(vec)
  } else {
    idxs_tail <- c()
  }

  idxs <- c( idxs_head, idxs_tail )

  return(idxs)

}

