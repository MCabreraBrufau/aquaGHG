


automaticflux <-  function(dataframe, myauxfile,
                           shoulder = 0,
                           gastype,
                           fluxSeparation,
                           displayPlots,
                           method){
  ## Check dataframe ####
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.null(dataframe) & !is.data.frame(dataframe)){
    stop("'dataframe' must be of class data.frame")}

  ## Check myauxfile ####
  if(missing(myauxfile)) stop("'myauxfile' is required")
  if(!is.null(myauxfile) & !is.data.frame(myauxfile)){
    stop("'myauxfile' must be of class data.frame")}

  ### Is there any match between dataframe and myauxfile?
  if(!any(unique(dataframe$UniqueID) %in% unique(myauxfile$UniqueID))){
    stop("'UniqueID' in 'myauxfile' has no match for 'UniqueID' in 'dataframe'")}

  ## Check shoulder ####
  if(is.null(shoulder)) stop("'shoulder' is required") else{
    if(!is.numeric(shoulder)) stop("'shoulder' must be of class numeric") else{
      if(shoulder < 0) stop("'shoulder' cannot be a negative value")}}

  ### gastype and match in dataframe ####
  if(missing(gastype)) stop("'gastype' is required")
  if(!is.null(gastype) & !is.character(gastype)) stop("'gastype' must be a character string")
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""),
                c("CO2dry_ppm", "COdry_ppb", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb", "H2O_ppm")))){
    stop("'gastype' must be of class character and one of the following: 'CO2dry_ppm', 'COdry_ppm', 'CH4dry_ppb', 'N2Odry_ppb', 'NH3dry_ppb' or 'H2O_ppm'")}
  if(!any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe)))){
    stop("'dataframe' must contain a column that matches 'gastype'")}
  if(any(grepl(paste("\\<", gastype, "\\>", sep = ""), names(dataframe))) &
     !is.numeric(dataframe[,gastype][[1]])){
    stop("The column that matches 'gastype' in 'dataframe' must be of class character")}

  ## Check fluxSeparation ####
  if(missing(fluxSeparation)) {
    fluxSeparation = FALSE
    warning("in automaticflux(), 'fluxSeparation' was not provided and automatically set to FALSE")
  }
  if(fluxSeparation != TRUE & quality.check != FALSE){
    stop("'fluxSeparation' must be TRUE or FALSE")}

  ## Check displayPlots ####
  if(missing(displayPlots)) {
    displayPlots = FALSE
    warning("in automaticflux(), 'displayPlots' was not provided and automatically set to FALSE")
  }
  if(displayPlots != TRUE & quality.check != FALSE){
    stop("'displayPlots' must be TRUE or FALSE")}

  ## Check method ####
  if(missing(method)) stop("in automaticflux(), 'method' is required")
  if(!any(grepl(paste("\\<", method, "\\>", sep = ""),
                c("trust.it.all", "focus.on.linear")))){
    stop("'method' must be of class character and one of the following: 'trust.it.all', 'focus.on.linear'")}


  # list of criteria for model selection
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")


  mydata_ow <- obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)

  # Join mydata_ow with info on start end incubation
  mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
    map_df(., ~as.data.frame(.x))

  # Additional auxiliary data required for flux calculation.
  mydata_auto <- mydata_auto %>%
    left_join(myauxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))

  if (method == "trust.it.all"){
    if(!fluxSeparation){
      # Calculate fluxes
      flux_auto <- goFlux(dataframe = mydata_auto, gastype)

      # Use best.flux to select the best flux estimates (LM or HM)
      # based on a list of criteria
      best.flux_auto <- best.flux(flux_auto, criteria)

      if(displayPlots){
        p <- flux.plot(
          flux.results = best.flux_auto, dataframe = mydata_auto,
          gastype = gastype, quality.check = TRUE,
          plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
          plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
        print(p)
      }
    } else { # in  that case we proceed with the separation between diffusion and ebullition

      best.flux_auto <- lapply(seq_along(mydata_ow), flux.separator.loop,
                                 list_of_dataframes = mydata_ow, gastype = gastype, auxfile = myauxfile, criteria)%>%
        map_df(., ~as.data.frame(.x))


      if(displayPlots){
        p <- separated.flux.plot(
          flux.results = best.flux_auto, dataframe = mydata_auto,
          gastype = gastype, quality.check = TRUE,
          plot.legend = c("SD"),
          plot.display = c("nb.obs"))
        print(p)
      }
    }
  } else if (method == "focus.on.linear"){

    if(fluxSeparation){
      warning("fluxSeparation ignored because method was set to 'focus.on.linear' and discards non-linear patterns in the observations.
... Change 'method' to 'trust.it.all' to activate fluxSeparation.")
      }

    lin.chunks <- lapply(seq_along(mydata_ow), find_first_linear_chunk.loop,
                         list_of_dataframe = mydata_ow, gastype = gastype, length.min=30) %>%
      map_df(., ~as.data.frame(.x))

    # for each incubation, extract data selected at previous step
    myauxfile_corr <- NULL
    for(id in unique(myauxfile$UniqueID)){
      myauxfile_corr.tmp <- myauxfile[myauxfile$UniqueID==id,]
      ind <- which(lin.chunks$UniqueID==id)
      myauxfile_corr.tmp$start.time <- lin.chunks$start.time[ind]
      myauxfile_corr.tmp$obs.length <- lin.chunks$obs.length[ind]
      myauxfile_corr.tmp$end.time <- lin.chunks$start.time[ind] + lin.chunks$obs.length[ind]
      myauxfile_corr <- rbind(myauxfile_corr, myauxfile_corr.tmp)
    }
    mydata_ow_corr <- obs.win(inputfile = dataframe, auxfile = myauxfile_corr, shoulder = 0)

    # Join mydata_ow with info on start end incubation
    mydiffusion_auto <- lapply(seq_along(mydata_ow_corr), join_auxfile_with_data.loop, flux.unique = mydata_ow_corr) %>%
      map_df(., ~as.data.frame(.x))

    # Calculate fluxes
    flux_diffusion <- goFlux(mydiffusion_auto, gastype)

    best.flux_auto <- best.flux(flux_diffusion, criteria)

    if(displayPlots){
      p <- flux.plot(
        flux.results = best.flux_auto, dataframe = mydata_auto,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
        plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
      print(p)
    }

  }

  return(best.flux_auto)
}

