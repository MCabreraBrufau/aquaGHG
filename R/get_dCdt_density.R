#' Title
#'
#' @param dataframe blabla
#' @param gastype blabla
#'
#' @return blabla
#'
#' @examples
#' blabla
#'
#' @export
get_dCdt_density <- function(dataframe, gastype){
  mydf <- data.frame(POSIX.time = dataframe$POSIX.time,
                     time = as.numeric(dataframe$POSIX.time-first(dataframe$POSIX.time)),
                     conc = dataframe[[gastype]])

  mydf <- mydf[!duplicated(mydf$time),]

  # smooth signal
  concsmooth <- smooth.spline(x = mydf$time, y = mydf$conc, nknots = round(length(mydf$time)/3), spar = 0.8)
  mydf$concsmooth <- approx(concsmooth$x, concsmooth$y, xout = mydf$time, rule = 2)$y

  # computing first derivative
  mydf$dydt <- get_dxdy(mydf$time, mydf$concsmooth)

  # computing density probability of first derivative
  d <- density(mydf$dydt)
  return(list(d, mydf))
}
