
# if CO2, selection of start stop
# if CH4, selection of start stop + selection of diffusive pattern


clickflux <- function(dataframe, myauxfile,
                      shoulder = 0,
                      gastype, plot.lim,
                      fluxSeparation,
                      displayPlots){
  ## Check arguments ####
  if(is.null(shoulder)) stop("'shoulder' is required") else{
    if(!is.numeric(shoulder)) stop("'shoulder' must be of class numeric") else{
      if(shoulder < 0) stop("'shoulder' cannot be a negative value")}}

  ## Check dataframe ####
  if(missing(dataframe)) stop("'dataframe' is required")
  if(!is.null(dataframe) & !is.data.frame(dataframe)){
    stop("'dataframe' must be of class data.frame")}

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
    stop("The column that matches 'gastype' in 'dataframe' must be of class numeric")}






  n_incubations <- length(unique(dataframe$UniqueID))

  mydata_ow <- obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)
  mydata_manID <- click.peak2(mydata_ow, seq = seq(1,n_incubations),
                              gastype = gastype,
                              plot.lim = plot.lim)

  # Calculate fluxes
  flux_manID <- goFlux(mydata_manID, gastype)

  # Use best.flux to select the best flux estimates (LM or HM)
  # based on a list of criteria
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

  best.flux_manID <- best.flux(flux_manID, criteria)

  if(displayPlots){
    p <- flux.plot(
      flux.results = best.flux_manID, dataframe = mydata_manID,
      gastype = gastype, quality.check = TRUE,
      plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
      plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
    print(p)
  }

  if(!fluxSeparation){
    return(best.flux_manID)
  } else { # in  that case we proceed with the separation between diffusion and ebullition
    message("Please now select diffusion patterns only")

    # for each incubation, extract data selected at previous step
    myauxfile_corr <- NULL
    for(id in unique(myauxfile$UniqueID)){
      myauxfile_corr.tmp <- myauxfile[myauxfile$UniqueID==id,]
      myauxfile_corr.tmp$start.time <- first(mydata_manID$start.time_corr[mydata_manID$UniqueID==id])
      myauxfile_corr.tmp$obs.length <- first(mydata_manID$obs.length_corr[mydata_manID$UniqueID==id])
      myauxfile_corr.tmp$end.time <- first(mydata_manID$end.time_corr[mydata_manID$UniqueID==id])
      myauxfile_corr <- rbind(myauxfile_corr, myauxfile_corr.tmp)
    }
    mydata_ow_corr <- obs.win(inputfile = dataframe, auxfile = myauxfile_corr, shoulder = 0)

    mydata_diffusionID <- click.peak2(mydata_ow_corr, seq = seq(1,n_incubations),
                                      gastype = gastype,
                                      plot.lim = plot.lim)

    # Calculate fluxes
    flux_diffusion <- goFlux(mydata_diffusionID, gastype)

    best.flux_diffusion <- best.flux(flux_diffusion, criteria)

    best.flux_manID <- best.flux_diffusion
    names(best.flux_manID)[-1] <- paste0(names(best.flux_manID)[-1],"_diffusion")

    # Estimating ebullition component
    for (id in unique(best.flux_manID$UniqueID)){
      i <- which(best.flux_manID$UniqueID == id)

      CH4_initial <-  flux_manID$C0[i]
      CH4_final <- flux_manID$Ct[i]
      incubation_time <- myauxfile_corr$obs.length[which(myauxfile_corr$UniqueID == id)]
      best.flux_manID$totalflux[i] <- (CH4_final-CH4_initial)/incubation_time*flux_manID$flux.term[i] # nmol/m2/s
      best.flux_manID$ebullitionflux[i] <- best.flux_manID$totalflux[i]  - best.flux_manID$LM.flux_diffusion[i] # total flux - diffusive term
      best.flux_manID$diffusionflux[i] <- best.flux_manID$LM.flux_diffusion[i]
    }

    if(displayPlots){
      p <- flux.plot(
        flux.results = best.flux_diffusion, dataframe = mydata_manID,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
        plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
      print(p)
    }

    return(best.flux_manID)
  }
}





