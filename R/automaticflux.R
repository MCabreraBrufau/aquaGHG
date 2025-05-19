
#' Automatic flux calculation
#'
#' @name automaticflux
#' @param dataframe a data.frame containing gas measurements (see \code{gastype}
#'                  below), water vapor measurements (see \code{H2O_col} below)
#'                  and the following columns: \code{UniqueID}, \code{Etime}, and
#'                  the precision of the instrument for each gas (see description below).
#'
#' The precision of the instrument is needed to restrict kappa-max
#' (\code{\link[goFlux]{k.max}}) in the non-linear flux calculation
#' (\code{\link[goFlux]{HM.flux}}). Kappa-max is inversely proportional to
#' instrument precision. If the precision of your instrument is unknown, it is
#' better to use a low value (e.g. 1 ppm) to allow for more curvature, especially
#' for water vapor fluxes, or very long measurements, that are normally curved.
#' The default values given for instrument precision are the ones provided by
#' the manufacturer upon request, for the latest model of this instrument
#' available at the time of the creation of this function (11-2023).
#'
#' @param myauxfile a data.frame containing auxiliary information needed with the
#'                  following columns: \code{UniqueID}, \code{start.time}, \code{obs.length},
#'                  \code{Vtot}, \code{Area}, \code{Pcham}, and \code{Tcham}.
#' @param shoulder numerical value; time before and after measurement in observation
#'                 window (seconds). Default is 30 seconds.
#' @param gastype character string; specifies which column should be used for the
#'                flux calculations. Must be one of the following: "CO2dry_ppm",
#'                "CH4dry_ppb", "COdry_ppb", "N2Odry_ppb", "NH3dry_ppb" or "H2O_ppm".
#' @param fluxSeparation logical; if \code{fluxSeparation = TRUE}, the model proceeds with
#'                        automatic separation between diffusion and ebullition fluxes.
#' @param force.separation logical; if TRUE, flux separation is forced to be performed,
#'                          regardless of the output of \code{find_bubbles}
#' @param displayPlots logical; if \code{displayPlots = TRUE}, plots showing how
#'                    the model performs are shown
#' @param method character string; specifies which method should be used to automatically
#'                process the gas measurements. Must be one of the following:
#'                "trust.it.all" (the model keeps the entire incubation time series)
#'                or "keep.it.linear" (the model first automatically selects the
#'                first linear chunk in the data and discard the rest). Note that
#'                \code{fluxSeparation} cannot be set to T if \code{method = 'keep.it.linear'}
#'
#' @return Returns a data frame with: a \code{UniqueID} per
#' measurement, 11 columns for the linear model results (linear flux estimate
#' (\code{\link[goFlux]{LM.flux}}), initial gas concentration
#' (\code{LM.C0}), final gas concentration (\code{LM.Ct}), slope of linear
#' regression (\code{LM.slope}), mean absolute error (\code{LM.MAE}), root mean
#' square error (\code{LM.RMSE}), Akaike's information criterion corrected for
#' small sample size (\code{LM.AICc}), standard error (\code{LM.SE}), relative
#' standard error (\code{LM.se.rel}), coefficient of determination (\code{LM.r2}),
#' and \emph{p-value} (\code{LM.p.val})), 11 columns for the non-linear model
#' results (non-linear flux estimate (\code{\link[goFlux]{HM.flux}}),
#' initial gas concentration (\code{HM.C0}), the assumed concentration of
#' constant gas source below the surface (\code{HM.Ci}), slope at \code{t=0}
#' (\code{HM.slope}), mean absolute error (\code{HM.MAE}), root mean square error
#' (\code{HM.RMSE}), Akaike's information criterion corrected for small sample
#' size (\code{HM.AICc}), standard error (\code{HM.SE}), relative standard error
#' (\code{HM.se.rel}), coefficient of determination (\code{HM.r2}), and curvature
#' (kappa; \code{HM.k}), as well as the minimal detectable flux
#' (\code{\link[goFlux]{MDF}}), the precision of the instrument
#' (\code{prec}), the flux term (\code{\link[goFlux]{flux.term}}),
#' kappa-max (\code{\link[goFlux]{k.max}}), the g factor (g.fact;
#' \code{\link[goFlux]{g.factor}}), the number of observations used
#' (\code{nb.obs}) and the true initial gas concentration (\code{C0}) and final
#' gas concentration (\code{Ct}). In case \code{fluxSeparation} was set to 'TRUE',
#' (\code{dataframe}) also contains information on flux separation, with flux value and
#' standard deviation estimated for total, diffusive, and ebullition fluxes.
#'
#' @examples
#' CH4_flux.auto <- automaticflux(dataframe = mydata, myauxfile = myauxfile, shoulder = 30, gastype = "CH4dry_ppb",
#'                               fluxSeparation = T, displayPlots = T,
#'                               method = "trust.it.all")
#'
#' @export
#'
automaticflux <- function(dataframe, myauxfile, shoulder, gastype,
                           fluxSeparation, force.separation,
                           displayPlots, method){
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

  ## gastype and match in dataframe ####
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
  if(fluxSeparation != TRUE & fluxSeparation != FALSE){
    stop("'fluxSeparation' must be TRUE or FALSE")}

  ## Check force.separation ####
  if(missing(force.separation)) {
    force.separation = FALSE
    # warning("in automaticflux(), 'force.separation' was not provided and automatically set to FALSE")
  }
  if(force.separation != TRUE & force.separation != FALSE){
    stop("'force.separation' must be TRUE or FALSE")}

  ## Check displayPlots ####
  if(missing(displayPlots)) {
    displayPlots = FALSE
    warning("in automaticflux(), 'displayPlots' was not provided and automatically set to FALSE")
  }
  if(displayPlots != TRUE & displayPlots != FALSE){
    stop("'displayPlots' must be TRUE or FALSE")}

  ## Check method ####
  if(missing(method)) stop("in automaticflux(), 'method' is required")
  if(!any(grepl(paste("\\<", method, "\\>", sep = ""),
                c("trust.it.all", "focus.on.linear")))){
    stop("'method' must be of class character and one of the following: 'trust.it.all', 'focus.on.linear'")}


  # list of criteria for model selection
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

  # making sure we work with the auxfile corresponding to measurements
  myauxfile <- myauxfile[myauxfile$UniqueID %in% unique(dataframe$UniqueID),]

  # split measurements by UniqueID
  mydata_ow <- my_obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)

  # Join mydata_ow with info on start end incubation
  mydata_auto <- lapply(seq_along(mydata_ow), join_auxfile_with_data.loop, flux.unique = mydata_ow) %>%
    map_df(., ~as.data.frame(.x))

  # Calculate fluxes
  flux_auto <- goFlux(dataframe = mydata_auto, gastype)

  # Use best.flux to select the best flux estimates (LM or HM)
  # based on a list of criteria
  best.flux_auto <- best.flux(flux_auto, criteria)

  #MIGUEL: When HM.flux results in NA, flux.plot functions fails (error: flux.plot(flux.results = best.flux_auto, dataframe = mydata_auto, :   #'HM.flux' in 'flux.results' must be of class numeric) )
  #MIGUEL: we try to work arround this by forcing as.numeric to the HM result columns that require numeric class
  best.flux_auto<- best.flux_auto %>% 
    mutate(across(c(HM.flux, HM.C0,HM.Ci,HM.k,HM.MAE,HM.RMSE,HM.AICc,HM.se.rel,HM.SE,HM.r2,g.fact), as.numeric))

  
  
  if (method == "trust.it.all"){

    if(!fluxSeparation){

      if(displayPlots){
        p <- flux.plot(
          flux.results = best.flux_auto, dataframe = mydata_auto,
          gastype = gastype, quality.check = TRUE,
          plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
          plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
        print(p)
      }

    } else { # in  that case we proceed with the separation between diffusion and ebullition

      best.flux_auto <- lapply(seq_along(mydata_ow), flux_separator.loop,
                               list_of_dataframes = mydata_ow, gastype = gastype, auxfile = myauxfile,
                               criteria = criteria, force.separation = force.separation,
                               mybest.flux = best.flux_auto)%>%
        map_df(., ~as.data.frame(.x))


      if(displayPlots){
        p <- separated_flux_plot(
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

    linear_chunk <- lapply(seq_along(mydata_ow), find_first_linear_chunk.loop,
                           list_of_dataframe = mydata_ow, gastype = gastype, length.min=30) %>%
      map_df(., ~as.data.frame(.x))

    # for each incubation, extract data selected at previous step
    myauxfile_corr <- NULL
    for(id in unique(myauxfile$UniqueID)){
      myauxfile_corr.tmp <- myauxfile[myauxfile$UniqueID==id,]
      ind <- which(linear_chunk$UniqueID==id)
      myauxfile_corr.tmp$start.time <- linear_chunk$start.time[ind]
      myauxfile_corr.tmp$obs.length <- linear_chunk$obs.length[ind]
      myauxfile_corr.tmp$end.time <- linear_chunk$start.time[ind] + linear_chunk$obs.length[ind]
      myauxfile_corr <- rbind(myauxfile_corr, myauxfile_corr.tmp)
    }
    mydata_ow_corr <- my_obs.win(inputfile = dataframe, auxfile = myauxfile_corr, shoulder = 0)

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

