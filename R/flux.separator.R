



flux.separator <- function(dataframe, gastype, auxfile, criteria, force.separation, mybest.flux){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[2],]
  id = unique(dataframe$UniqueID)
  auxfile_corr <- auxfile[auxfile$UniqueID==id,]
  mybest.flux <- mybest.flux[mybest.flux$UniqueID == id,]

  # Initializing variables
  mybest.flux$total.flux <- mybest.flux$ebullition.flux <- mybest.flux$diffusion.flux <-
    mybest.flux$obs.length_diffusion <- mybest.flux$total.flux.SD <-
    mybest.flux$diffusion.flux.SD <- mybest.flux$ebullition.flux.SD <- NA



  # if CH4 is available, first we check if it makes sense trying to separate fluxes based on find_bubbles()
  # if find_bubbles returns a NULL, flux.separator() automatically skips the separation
  # and calculate total fluxes as diffusion only.
  # if user specifies `force.separation = TRUE`, the separation is done regardless
  # of the result from find_bubbles()

  # if CH4 in dataframe, check bubbles
  if(any(grepl(paste0("\\<CH4dry_ppb\\>"), names(dataframe)))){
    bubbles <- find_bubbles(time = dataframe$Etime,
                            conc = dataframe$CH4dry_ppb, window.size = 10)
  } else {
    warning(paste0("For ",id, ", CH4 measurements are unavailable, `find_bubbles()` cannot look for potential bubbles.
                   `bubbles` was set to NULL"))
    bubbles = NULL
  }

  if (is.null(bubbles) & force.separation==FALSE){
    warning(paste0("For ",id, ", no bubbles were found. Total flux is attributed to diffusion only"))

    mybest.flux$obs.length_diffusion <- auxfile_corr$obs.length

    mybest.flux$total.flux <- mybest.flux$best.flux
    mybest.flux$ebullition.flux <- 0
    mybest.flux$diffusion.flux <- mybest.flux$best.flux

    SD_total.flux <- SD_diffusion.flux <- mybest.flux$LM.SE*sqrt(mybest.flux$nb.obs)
    SD_ebullition.flux <- NA

  } else {
    if (is.null(bubbles)){warning(paste0("For ",id, ", no bubbles were found but flux separation was forced with `force.separation = TRUE`.
                                         Consider setting force.separation to FALSE."))}

    # Automatically identify the first linear chunk
    linear_chunk <- find_first_linear_chunk(dataframe = dataframe, gastype = gastype, length.min = 10)

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
    auxfile_corr$start.time <- linear_chunk$start.time
    auxfile_corr$obs.length <- linear_chunk$obs.length
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

    mybest.flux_diffusion <- best.flux(flux_diffusion, criteria)

    mybest.flux <- mybest.flux_diffusion
    # names(mybest.flux)[-1] <- paste0(names(mybest.flux)[-1],"_diffusion")

    # adding information of data used for diffusion model
    mybest.flux$obs.length_diffusion <- linear_chunk$obs.length

    mybest.flux$total.flux <- (Cf-C0)/incubation_time*mybest.flux_diffusion$flux.term # nmol/m2/s
    mybest.flux$ebullition.flux <- mybest.flux$total.flux - mybest.flux_diffusion$best.flux # nmol/m2/s
    mybest.flux$diffusion.flux <- mybest.flux_diffusion$best.flux # nmol/m2/s

    # Error propagation (expressed as SD)
    SD_C0 <- sd(mydf$conc[mydf$time<t.win])
    SD_Cf <- sd(mydf$conc[mydf$time>max(mydf$time)-t.win])
    deltaconcs = Cf-C0
    SD_deltaconcs <- sqrt(SD_C0^2+SD_Cf^2)
    SD_total.flux <- abs(mybest.flux$total.flux) * SD_deltaconcs/deltaconcs
    if(mybest.flux_diffusion$model == "LM"){SE_diffusion.flux = mybest.flux_diffusion$LM.SE} else {SE_diffusion.flux = mybest.flux_diffusion$HM.SE}
    SD_diffusion.flux <- mybest.flux_diffusion$LM.SE*sqrt(mybest.flux_diffusion$nb.obs)
    SD_ebullition.flux <- sqrt(SD_diffusion.flux^2+SD_total.flux^2)

    # warnings
    if( mybest.flux$ebullition.flux < abs(mybest.flux$diffusion.flux+SD_diffusion.flux)){
      warning(paste0("for ",id, ", ebullition term is within range of uncertainty of diffusion."))
      mybest.flux$quality.check <- "ebullition too low to be trusted"
    }

    if( mybest.flux$ebullition.flux < 0){
      warning(paste0("for ",id, ", negative ebullition term. It was forced to 0."))
      mybest.flux$ebullition.flux <- 0
    }

    if( mybest.flux$diffusion.flux > mybest.flux$total.flux){
      warning(paste0("for ",id, ", diffusion term is larger than total flux estimated."))
      mybest.flux$quality.check <- "diffusion > total flux"
    }
  }

  mybest.flux$total.flux.SD <- SD_total.flux
  mybest.flux$diffusion.flux.SD <- SD_diffusion.flux
  mybest.flux$ebullition.flux.SD <- SD_ebullition.flux

  return(mybest.flux)
}

flux.separator.loop <-  function(x, list_of_dataframes, gastype, auxfile, criteria, force.separation, mybest.flux) {

  # function to apply in the loop. Adapt parameters to your needs.
  mybest.flux <- flux.separator(dataframe = list_of_dataframes[[x]], gastype, auxfile, criteria, force.separation, mybest.flux)

  return(mybest.flux)
}
