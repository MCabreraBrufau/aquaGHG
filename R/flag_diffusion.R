

flag_diffusion <- function(dataframe, kstar){
  dataframe <- dataframe[which(dataframe$UniqueID==first(dataframe$UniqueID)),]

  # initialization
  dataframe$is_diffusion <- TRUE


  # computing density probability of first derivative
  d_df <- get_dCdt_density(dataframe)
  d <- d_df[[1]]
  mych4 <- d_df[[2]]

  # if dydt doesn't exceed 10 ppb/sec at any point, we consider ebullition is unlikely
  if(sum(mych4$dydt>10)>1){

    # computing density probability of first derivative
    d <- d_df[[1]]
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

    dataframe$is_diffusion <- FALSE # if not diffusion, it is necessarily ebullition i.e. is_diffusion = FALSE
    dataframe$is_diffusion[mych4$dydt>lower_bound & mych4$dydt<upper_bound] <- TRUE
  }

  return(dataframe$is_diffusion)
}








