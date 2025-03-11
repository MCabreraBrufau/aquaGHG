


automaticflux <-  function(dataframe, myauxfile,
                           shoulder = 0,
                           gastype,
                           fluxSeparation,
                           displayPlots,
                           method){
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
    stop("The column that matches 'gastype' in 'dataframe' must be of class character")}


  mydata_ow <- obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)

  if (method == "trust.it.all"){

    # Join mydata_ow with info on start end incubation
    mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
      map_df(., ~as.data.frame(.x))

    # Additional auxiliary data required for flux calculation.
    mydata_auto <- mydata_auto %>%
      left_join(myauxfile %>% select(UniqueID, Area, Vtot, Tcham, Pcham))

    # Calculate fluxes
    flux_auto <- goFlux(dataframe = mydata_auto, gastype)

    # Use best.flux to select the best flux estimates (LM or HM)
    # based on a list of criteria
    criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")
    best.flux_auto <- best.flux(flux_auto, criteria)

    if(displayPlots){
      p <- flux.plot(
        flux.results = best.flux_auto, dataframe = mydata_auto,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
        plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
      print(p)
    }

    if(!fluxSeparation){
      return(best.flux_auto)
    } else { # in  that case we proceed with the separation between diffusion and ebullition

      separated.fluxes <- lapply(seq_along(mydata_ow), flux.separator.loop,
                                 list_of_dataframes = mydata_ow, gastype = gastype, kstar = 0.5)%>%
        map_df(., ~as.data.frame(.x))

      # Estimating ebullition component
      for (id in separated.fluxes$UniqueID){
        i <- which(separated.fluxes$UniqueID == id)

        CH4_initial <-  separated.fluxes$C0[i]
        CH4_final <- separated.fluxes$Ct[i]
        incubation_time <- separated.fluxes$duration[i]
        myflux.term <- flux_auto$flux.term[which(flux_auto$UniqueID == id)]
        separated.fluxes$totalflux[i] <- (CH4_final-CH4_initial)/incubation_time*myflux.term # nmol/m2/s
        separated.fluxes$diffusionflux[i] <-  separated.fluxes$avg_diff_slope[i]*myflux.term # nmol/m2/s
        separated.fluxes$ebullitionflux[i] <- separated.fluxes$totalflux[i]  - separated.fluxes$diffusionflux[i] # total flux - diffusive term
      }

      best.flux_auto$totalflux <- separated.fluxes$totalflux
      best.flux_auto$diffusionflux <- separated.fluxes$diffusionflux
      best.flux_auto$ebullitionflux <- separated.fluxes$ebullitionflux


      if(displayPlots){
        p <- flux.plot(
          flux.results = best.flux_auto, dataframe = mydata_auto,
          gastype = gastype, quality.check = TRUE,
          plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
          plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
        print(p)
      }
      return(best.flux_auto)
    }
  } else if (method == "focus on linear"){
        lin.chunks <- lapply(seq_along(mydata_ow), find_linear_chunk.loop,
                             list_of_dataframe = mydata_ow, gastype = gastype, kstar = 0.4, which.chunk = "longest", length.min=30) %>%
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

        best.flux_diffusion <- best.flux(flux_diffusion, criteria)

        if(displayPlots){
          p <- flux.plot(
            flux.results = best.flux_diffusion, dataframe = dataframe,
            gastype = gastype, quality.check = TRUE,
            plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
            plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
          print(p)
        }

  #
  #       # Estimating ebullition component
  #       for (id in unique(best.flux_auto$UniqueID)){
  #         i <- which(best.flux_auto$UniqueID == id)
  #
  #         CH4_initial <-  flux_manID$C0[i]
  #         CH4_final <- flux_manID$Ct[i]
  #         incubation_time <- myauxfile_corr$obs.length[which(myauxfile_corr$UniqueID == id)]
  #         best.flux_auto$totalflux[i] <- (CH4_final-CH4_initial)/incubation_time*flux_manID$flux.term[i] # nmol/m2/s
  #         best.flux_auto$ebullitionflux[i] <- best.flux_auto$totalflux[i]  - best.flux_auto$LM.flux_diffusion[i] # total flux - diffusive term
  #         best.flux_auto$diffusionflux[i] <- best.flux_auto$LM.flux_diffusion[i]
  #       }
  }

}

