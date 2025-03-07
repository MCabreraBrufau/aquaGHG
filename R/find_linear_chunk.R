


find_linear_chunk <- function(dataframe, gastype, kstar = 0.4, which.chunk = "first"){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[3],]

  mydf <- data.frame(POSIX.time = dataframe$POSIX.time,
                     time = as.numeric(dataframe$POSIX.time-first(dataframe$POSIX.time)),
                     conc = dataframe[[gastype]])

  # smooth signal
  concsmooth <- smooth.spline(x = mydf$time, y = mydf$conc, nknots = round(length(mydf$time)/3), spar = 0.8)
  mydf$concsmooth <- approx(concsmooth$x, concsmooth$y, xout = mydf$time, rule = 2)$y

  # computing first derivative
  mydf$dydt <- get_dxdy(mydf$time, mydf$concsmooth)

  # computing density probability of first derivative
  d <- density(mydf$dydt)
  half_dens_max <- kstar*max(d$y)
  ind_over_half_dens_max <- which(d$y>half_dens_max)

  # we define lower and upper boundaries of the "main slope" as the slopes with
  # a density of probability higher than half of the maximum density
  lower_bound <- d$x[first(ind_over_half_dens_max)]
  upper_bound <- d$x[last(ind_over_half_dens_max)]

  avg_slope <- mean(d$x[ind_over_half_dens_max])
  sd_slope <- sd(d$x[ind_over_half_dens_max])

  delta_C <- last(mydf$concsmooth) - min(mydf$concsmooth)
  duration <- last(mydf$time)


  df_out <- data.frame(duration = last(mydf$time),
                       delta_conc = last(mydf$concsmooth) - min(mydf$concsmooth),
                       avg_diff_slope = avg_slope,
                       sd_diff_slope = sd_slope)

  mydf$is_linear <- F
  mydf$is_linear[mydf$dydt>lower_bound & mydf$dydt<upper_bound] <- T

  ggplot(mydf, aes(dydt))+
    geom_rect(aes(xmin = lower_bound, ymin = -Inf, xmax = upper_bound, ymax = Inf), fill = "#FF7F50", alpha=0.1)+
    geom_hline(yintercept = half_dens_max, alpha=0.5)+
    geom_density()+
    xlab("Rate of change")+
    theme_article()

  ggplot(mydf)+
    geom_path(aes(time, dydt), linewidth=1)+
    theme_article()

  ggplot(mydf, aes(time, dydt))+geom_point(aes( colour = is_linear))+
    theme_article()

  ggplot(mydf, aes(time, conc))+geom_point(aes( colour = is_linear))+
    geom_path(aes(time, concsmooth), linewidth=2)+
    theme_article()

  vect <- mydf$is_linear
  jumps <- NULL
  if(vect[1]==T){jumps <- c(1)}
  for(i in seq(2,length(vect))){
    if(vect[i] != vect[i-1]){
      jumps <- c(jumps, i)
    }
  }

  chunks <- data.frame(start = jumps[seq(1,length(jumps),2)],
                       end = jumps[seq(2,length(jumps),2)])
  chunks$obs.length <- chunks$end-chunks$start

  chunks$start.time <- mydf$POSIX.time[chunks$start]


  if (which.chunk == "first"){chunk.out <- chunks[1,]}
  if (which.chunk == "last"){chunk.out <- chunks[length(chunks$start),]}
  if (which.chunk == "longest"){chunk.out <- chunks[which.max(chunks$obs.length),]}

  return(chunk.out)
}
