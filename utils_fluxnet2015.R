ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350) # day of year of middle-month-day

interpol_field <- function( time0, time1, field0, field1, outtime ){
  outfield <- field0 + (field1-field0)/(time1-time0) * (outtime - time0)
  return(outfield)
}

add_alpha <- function( col, alpha ){
  ## add alpha to color given as a name
  col    <- col2rgb( col, alpha=TRUE )/255
  col[4] <- alpha
  col    <- rgb(col[1,],col[2,],col[3,],col[4,])
  return( col )
}

myrmse <- function( actual, predicted ){
  library(Metrics)
  actual    <- actual[ !is.na(actual) & !is.na(predicted) ]
  predicted <- predicted[ !is.na(actual) & !is.na(predicted) ]
  return( rmse( actual, predicted ) )
}

get_doy <- function( year_dec ){
  doy <- rep( NA, length(year_dec))
  for (idx in 1:length(year_dec)){
    doy[idx] <- round(( year_dec[idx] - floor(year_dec[idx]) )*365 + 1)    
  }
  return( doy )
}

get_year <- function( year_dec ){
  year <- floor( year_dec)
  return( year )
}

get_dom <- function( doy ){
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  cumndaymonth <- cumsum( ndaymonth )
  dom <- rep( NA, length(doy))
  for (idx in 1:length(doy)){
    dom[idx] <- doy[idx] - sum( ndaymonth[ which( doy[idx] > cumndaymonth ) ] )
  }
  return( dom )
}

get_moy <- function( doy ){
  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  cumndaymonth <- cumsum( ndaymonth )
  moy <- rep( NA, length(doy))
  for (idx in 1:length(doy)[1]){
    moy[idx] <- min( which( doy[idx] <= cumndaymonth ) )
  }
  return( moy )
}

mycurve <- function( func, from, to ){

  range_x <- seq( from, to, by=(to-from)/100 )
  range_y <- sapply( range_x, func )

  plot( range_x, range_y, type="l" )

}

init_daily_dataframe <- function( yrstart, yrend ){

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)

  ndayyear <- sum(ndaymonth) 
  nmonth   <- length(ndaymonth)
  nyrs    <- length( yrstart:yrend )
  dm   <- rep( NA, sum(ndaymonth)*length(yrstart:yrend) )
  jdx <- 0
  for (yr in yrstart:yrend ){
    for (imoy in 1:nmonth){
      for (idm in 1:ndaymonth[imoy]){
        jdx <- jdx + 1 
        dm[jdx]   <- idm
      }
    }
  }
  ddf <- data.frame( 
    doy=rep( seq(ndayyear), nyrs ), 
    moy=rep( rep( seq(nmonth), times=ndaymonth ), times=nyrs ),
    dom=dm,
    year=rep( yrstart:yrend, each=ndayyear ) 
  )
  ddf$date <- as.POSIXlt( as.Date( paste( as.character(ddf$year), "-01-01", sep="" ) ) + ddf$doy - 1 )
  ddf$year_dec <- ddf$year + ( ddf$doy - 1 ) / sum( ndaymonth )

  return( ddf )
}

init_monthly_dataframe <- function( yrstart, yrend ){

  ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  middaymonth <- c(16,44,75,105,136,166,197,228,258,289,319,350)
  nyrs    <- length( yrstart:yrend )

  mdf <- data.frame( 
    year=rep(yrstart:yrend,each=12) , 
    moy=rep(1:12,nyrs), 
    doy=rep(middaymonth,nyrs)
    )
  mdf$date <- as.POSIXlt( as.Date( paste( as.character(mdf$year), "-01-01", sep="" ) ) + mdf$doy - 1 )

  return( mdf)

}

## /////////////////////////////////////////////////////////////////////////
## Function 'mycolorbar' draws a colorbar based on the image function
## 'layout' must be called before to arrange the colorbar in a subplot as
## defined by 'layout'.
## Beni Stocker, 3.6.2013
## -------------------------------------------------------------------------
mycolorbar <- function( col,           # a vector of colors from which to interpolate
                       lev,            # levels of colorbar either in the form of c(min,max, N levels) or by a vector of length > 3 containing margins of levels
                       plot=TRUE,      # if false, then no colorbar is plotted but colors, levels and margins are returned
                       alpha=NA,       # transparency value, 0=100% transparent, 1=0% transparent
                       orient="h",     # orentation of colorbar
                       maxval=NA,      # maximum value, overrides upper margin
                       minval=NA,      # minimum value, overrides lower margin 
                       mincol=NA,      # overrides color for lowest level
                       dolabels=TRUE,  # add labels for margins to colorbar
                       doticks=TRUE,   # add tick marks at margins to colorbar
                       cex.axis=1.0,   # magnification of axis tickmarks
                       cex.lab=1.0     # magnification of axis labels
                       ) {
  library(gplots)

  if (length(lev)>3){
    explicit <- TRUE
  } else {
    explicit <- FALSE
  }

  if (explicit){

    ## Assume 'lev' declares explicit margins
    print("MYCOLORBAR: assuming explicit margins provided")

    len <- length(lev)
    # print(paste("len=",len))

    # margins.eff is used for labels at color key and is returned as $margins
    margins.eff <- lev
    # print(paste("length of margins.eff",length(margins.eff)))
    # print(margins.eff)

    # margins defines where margins.eff are to be labelled
    margins <- seq( from=0, to=(len-1), by=1 )
    margins.lab <- margins.eff
    # print(paste("length of margins",length(margins)))
    # print(margins)

    if (!is.na(maxval)) {
      margins.eff[length(margins)] <- maxval
    } 
    if (!is.na(minval)){
      margins.eff[1] <- minval
    }    

    ## Define color key centers (mid-points between margins)
    centers <- seq( from=0.5, to=(len-1.5), by=1 )
    # print(paste("length of centers",length(centers)))
    # print(centers)

    ## Define color range
    colors  <- colorRampPalette( col )( length(centers) )

  } else {
    
    ## Assume 'lev' declares (min,max,number of levels)
    ## Define color key margins
    len <- lev[3]
    margins <- seq( from=lev[1], to=lev[2], by=(lev[2]-lev[1])/len )
    margins.eff <- margins
    margins.lab <- margins
    if (!is.na(maxval)) {
      margins.eff[length(margins)] <- maxval
      margins.lab[length(margins)] <- margins[length(margins)]
    } 
    if (!is.na(minval)){
      margins.eff[1] <- minval
      margins.lab[1] <- margins[1]
    }
    
    ## Define color key centers (mid-points between margins)
    centers <- seq( from=lev[1]+(lev[2]-lev[1])/(2*len), to=lev[2]-(lev[2]-lev[1])/(2*len), by=(lev[2]-lev[1])/len )

    ## Define color range
    colors  <- colorRampPalette( col )( lev[3] )
  
  }

  if (!is.na(mincol)){
    colors[1] <- col2hex(mincol)
  }

  ## Alpha is transparency value
  if (!is.na(alpha)){
    colors <- col2rgb( colors, alpha=TRUE )/255
    colors[4,] <- alpha
    colors <- rgb(colors[1,],colors[2,],colors[3,],colors[4,])
  }

  if (plot) {

    if (dolabels==FALSE) {
      labels=FALSE
    } else {
      if (orient=="h"){
        if (explicit) {
          labels <- as.character(margins.lab)
        } else {
          labels <- as.character(margins.lab)
        }
      } else if (orient=="v") {
        if (explicit) {
          labels <- as.character(margins.lab) 
        } else {
          labels <- as.character(margins.lab)
        }
      } else {
        print("argument 'orient' must be either 'v' for vertical or 'h' for horizontal.")
      }
    }

    if (orient=="h"){
      if (explicit) {
        # xlim <- c(lev[1],lev[len])
        xlim <- c(margins[1],margins[length(margins)])
        image( centers, 0.5, as.matrix(centers), xlim=xlim, col=colors, axes=FALSE, ylab="", xlab="", cex.axis=cex.axis, cex.lab=cex.lab )
        box()
      } else {
        xlim <- c(lev[1],lev[2])
        image( centers, 0.5, as.matrix(centers), col=colors, axes=FALSE, xlim=xlim, ylab="", xlab="", cex.axis=cex.axis, cex.lab=cex.lab )
        box()
      }
      box()
    } else if (orient=="v") {
      if (explicit) {
        ylim <- c(margins[1],margins[length(margins)])
        image( 0.5, centers, as.matrix(t(centers)), col=colors, axes=FALSE, ylim=ylim, xlab="",ylab="", cex.axis=cex.axis, cex.lab=cex.lab )
        box()
        # Hack by substracting 1 in the following 2 lines (commented out original lines above)
        # ylim <- c(lev[1],lev[len])
        # image( 0.5, centers-1, as.matrix(t(centers)), col=colors, axes=FALSE, ylim=ylim, xlab="",ylab="", cex.axis=cex.axis, cex.lab=cex.lab )
        # axis( 2, at=margins-1, labels=as.character(lev) )
        } else {
        ylim <- c(lev[1],lev[2])
        image( 0.5, centers, as.matrix(t(centers)), col=colors, axes=FALSE, ylim=ylim, xlab="",ylab="", cex.axis=cex.axis, cex.lab=cex.lab )
        box()
      }
      box()      
    } else {
      print("argument 'orient' must be either 'v' for vertical or 'h' for horizontal.")
    }

    if (doticks) {
      if (orient=="h"){
        axis( 1, at=margins, labels=labels, cex.axis=cex.axis, cex.lab=cex.lab )
      } else if (orient=="v") {
        axis( 2, at=margins, labels=labels, cex.axis=cex.axis, cex.lab=cex.lab )
      }
    }

  }
  
  out.mycolorbar <- list()
  out.mycolorbar$colors <- colors
  out.mycolorbar$margins <- margins.eff
  out.mycolorbar$centers <- centers
  
  return(out.mycolorbar)

}



