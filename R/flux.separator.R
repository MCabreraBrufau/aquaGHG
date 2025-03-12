






flux.separator <- function(dataframe, gastype, auxfile, criteria){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[2],]
  id = unique(dataframe$UniqueID)

  first_lin_chunk <- find_first_linear_chunk(dataframe = dataframe, gastype = gastype, length.min=30)


  # computing density probability of first derivative
  d_df <- get_dCdt_density(dataframe, gastype)
  d <- d_df[[1]]
  mydf <- d_df[[2]]

  kstar = 0.5

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

  t.win <- 30
  C0 = min(mydf$concsmooth[mydf$time<t.win])
  Cf = max(mydf$concsmooth[mydf$time>max(mydf$time)-t.win])
  incubation_time = last(mydf$time)

  # for each incubation, extraCf data selected at previous step
  auxfile_corr <- auxfile[auxfile$UniqueID==id,]
  auxfile_corr$start.time <- first_lin_chunk$start.time
  auxfile_corr$obs.length <- first_lin_chunk$obs.length
  auxfile_corr$end.time <- auxfile_corr$start.time + auxfile_corr$obs.length

  mydata_ow_corr <- obs.win(inputfile = dataframe, auxfile = auxfile_corr, shoulder = 0)

  # Join mydata_ow with info on start end incubation
  mydiffusion_auto <- lapply(seq_along(mydata_ow_corr), join_auxfile_with_data.loop, flux.unique = mydata_ow_corr) %>%
    map_df(., ~as.data.frame(.x))

  # Additional auxiliary data required for flux calculation.
  mydiffusion_auto <- mydiffusion_auto %>%
    left_join(myauxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))


  # Calculate fluxes
  flux_diffusion <- goFlux(mydiffusion_auto, gastype)

  best.flux_diffusion <- best.flux(flux_diffusion, criteria)

  best.flux_auto <- best.flux_diffusion
  # names(best.flux_auto)[-1] <- paste0(names(best.flux_auto)[-1],"_diffusion")

  best.flux_auto$total.flux <- (Cf-C0)/incubation_time*best.flux_diffusion$flux.term # nmol/m2/s
  best.flux_auto$ebullition.flux <- best.flux_auto$total.flux - best.flux_diffusion$best.flux # nmol/m2/s
  best.flux_auto$diffusion.flux <- best.flux_diffusion$best.flux # nmol/m2/s

  # adding information of data used for diffusion model
  best.flux_auto$obs.length_diffusion <- first_lin_chunk$obs.length

  # Error propagation (expressed as SD)
  deltaC0 <- sd(mydf$conc[mydf$time<t.win])
  deltaCf <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
  deltaconcs = Cf-C0
  SD_deltaconcs <- sqrt(deltaC0^2+deltaCf^2)
  SD_total.flux <- abs(best.flux_auto$total.flux) * SD_deltaconcs/deltaconcs
  if(best.flux_auto$model == "LM"){SE_diffusion.flux = best.flux_auto$LM.SE} else {SE_diffusion.flux = best.flux_auto$HM.SE}
  SD_diffusion.flux <- best.flux_auto$LM.SE*sqrt(best.flux_auto$nb.obs)
  SD_ebullition.flux <- sqrt(SD_diffusion.flux^2+SD_total.flux^2)

  best.flux_auto$total.flux.SD <- SD_total.flux
  best.flux_auto$diffusion.flux.SD <- SD_diffusion.flux
  best.flux_auto$ebullition.flux.SD <- SD_ebullition.flux

  # warnings
  if( best.flux_auto$ebullition.flux < abs(best.flux_auto$diffusion.flux+best.flux_auto$diffusion.flux.SD)){
    warning(paste0("for ",id, ", ebullition term is within range of uncertainty of diffusion."))
    best.flux_auto$quality.check <- "ebullition too low to be trusted"
  }

  if( best.flux_auto$ebullition.flux < 0){
    warning(paste0("for ",id, ", negative ebullition term. It was forced to 0."))
    best.flux_auto$ebullition.flux <- 0
  }

  if( best.flux_auto$diffusion.flux > best.flux_auto$total.flux){
    warning(paste0("for ",id, ", diffusion term is larger than total flux estimated."))
    best.flux_auto$quality.check <- "diffusion > total flux"
  }

  return(best.flux_auto)
}

flux.separator.loop <-  function(x, list_of_dataframes, gastype, auxfile, criteria) {

  # function to apply in the loop. Adapt parameters to your needs.
  best.flux_auto <- flux.separator(dataframe = list_of_dataframes[[x]], gastype, auxfile, criteria)

  return(best.flux_auto)
}
