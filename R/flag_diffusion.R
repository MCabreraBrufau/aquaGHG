

flag_diffusion <- function(dataframe, kstar){
  dataframe <- dataframe[which(dataframe$UniqueID==first(dataframe$UniqueID)),]
  mych4 <- data.frame(time = as.numeric(dataframe$POSIX.time-first(dataframe$POSIX.time)),
                      ch4 = dataframe$CH4dry_ppb)

  ch4smooth <- smooth.spline(x = mych4$time, y = mych4$ch4, nknots = round(length(mych4$time)/3), spar = 0.6)

  mych4$ch4smooth <- approx(ch4smooth$x, ch4smooth$y, xout = mych4$time, rule = 2)$y

  # computing first derivative
  mych4$dydt <- get_dxdy(mych4$time, mych4$ch4smooth)

  # computing density probability of first derivative
  d <- density(mych4$dydt)
  half_dens_max <- kstar*max(d$y)
  ind_over_half_dens_max <- which(d$y>half_dens_max)

  # we define lower and upper boundaries of the "main slope" as the slopes with
  # a density of probability higher than half of the maximum density
  lower_bound <- d$x[first(ind_over_half_dens_max)]
  upper_bound <- d$x[last(ind_over_half_dens_max)]

  avg_slope <- mean(d$x[ind_over_half_dens_max])
  sd_slope <- sd(d$x[ind_over_half_dens_max])


  mych4_sel <- mych4[mych4$dydt>lower_bound & mych4$dydt<upper_bound,]
  my_c0 <- min(mych4$ch4smooth[seq(1,100)])

  dataframe$is_diffusion <- FALSE
  dataframe$is_diffusion[mych4$dydt>lower_bound & mych4$dydt<upper_bound] <- TRUE

  return(dataframe$is_diffusion)
}








