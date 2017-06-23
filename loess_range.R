## loess and error curves almost just like ggplot2
## taken from http://stackoverflow.com/questions/3175033/recreate-ggplots-geom-smooth-ci-background-in-r-basic
loess_range <- function( x, y, col=rgb(0.1,0.1,0.1,0.25), dens=1e4 ){

  idxs <- which(!is.na(x)&!is.na(y))
  x <- x[idxs]
  y <- y[idxs]

  m   <- loess(y~x)
  xx  <- seq(min(x, na.rm=T), max(x, na.rm=T), (max(x, na.rm=T)-min(x, na.rm=T))/dens) #increase density of values to predict over to increase quality of curve
  f   <- predict(m, xx, se = TRUE)
  ci  <- f$se * qt( 0.975, f$df )
  cih <- f$fit + ci
  cil <- f$fit - ci
  # plot(x,y, ylim = c(min(cil,y), max(cih,y)), cex.axis = 0.85, xlab = '', ylab = '', type = 'n')
  # title(xlab = 'x', ylab = 'y',line = 2)
  # grid(col = 'gray')
  # points(x,y, pch = 19, cex = 0.65)
  lines( xx, f$fit, col = col, lwd = 1.2)
  xx <- c(xx, rev(xx))

  yy <- c(cil, rev(cih))

  polygon(xx, yy, col=add_alpha(col,0.4), border = NA)

}