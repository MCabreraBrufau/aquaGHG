






flux.separator <- function(dataframe, gastype, kstar){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[3],]


  # computing density probability of first derivative
  d_df <- get_dCdt_density(dataframe)
  d <- d_df[[1]]
  mydf <- d_df[[2]]

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

  df_out <- data.frame(UniqueID = unique(dataframe$UniqueID),
                       duration = last(mydf$time),
                       C0 = min(mydf$concsmooth),
                       Ct = last(mydf$concsmooth),
                       delta_conc = last(mydf$concsmooth) - min(mydf$concsmooth),
                       avg_diff_slope = avg_slope,
                       sd_diff_slope = sd_slope)

  return(df_out)
}

flux.separator.loop <-  function(x, list_of_dataframes, gastype, kstar) {

  # Function to apply in the loop. Adapt parameters to your needs.
  chunk.out <- flux.separator(dataframe = list_of_dataframes[[x]], gastype, kstar)

  return(chunk.out)
}
