

#' Separation between diffusion and ebullition patterns in CH4 measurements
#'
#' @param my_incub
#' @param UniqueID
#' @param kstar
#' @param doPlot
#'
#' @return
#' @export
#'
#' @examples
#'
flux_separator <- function(my_incub, UniqueID, kstar, doPlot){
  my_incub <- my_incub[which(my_incub$UniqueID==UniqueID),]
  mych4 <- data.frame(time = as.numeric(my_incub$POSIX.time-first(my_incub$POSIX.time)),
                      ch4 = my_incub$CH4dry_ppb)

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

  my_incub$is_diffusion <- FALSE
  my_incub$is_diffusion[mych4$dydt>lower_bound & mych4$dydt<upper_bound] <- TRUE

  my_incub$Etime <- as.numeric(my_incub$POSIX.time-first(my_incub$POSIX.time))


  if(doPlot){

    p_density <- ggplot(mych4, aes(dydt))+
      geom_rect(aes(xmin = lower_bound, ymin = -Inf, xmax = upper_bound, ymax = Inf), fill = "#FF7F50", alpha=0.1)+
      geom_hline(yintercept = half_dens_max, alpha=0.5)+
      geom_density()+
      xlab("Rate of change of CH4 [ppb/secs]")+
      theme_article()+
      ggtitle(UniqueID)

    p_fit <- ggplot(my_incub, aes(Etime, CH4dry_ppb, colour = is_diffusion))+geom_point()+theme_article()

    ggarrange(p_density, p_fit, nrow = 2, ncol = 1)

  }
  return(my_incub)
}








