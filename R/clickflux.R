#' Manual selection of valid data for GHG experts
#'
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
#'
#' @param myauxfile a data.frame containing auxiliary information needed with the
#'                  following columns: \code{UniqueID}, \code{start.time}, \code{obs.length},
#'                  \code{Vtot}, \code{Area}, \code{Pcham}, and \code{Tcham}.
#' @param shoulder numerical value; time before and after measurement in observation
#'                 window (seconds). Default is 30 seconds.
#' @param gastype character string; specifies which column should be used for the
#'                flux calculations. Must be one of the following: "CO2dry_ppm",
#'                "CH4dry_ppb", "COdry_ppb", "N2Odry_ppb", "NH3dry_ppb" or "H2O_ppm".
#' @param plot.lim numerical vector of length 2; sets the Y axis limits in the
#'                 plots. Default values are set for a typical gas measurement
#'                 of "CO2dry_ppm" in aquatic ecosystems: \code{plot.lim = c(380,1000)}.
##' @param fluxSeparation logical; if \code{fluxSeparation = TRUE}, the model proceeds with
#'                        automatic separation between diffusion and ebullition fluxes.
#' @param displayPlots logical; if \code{displayPlots = TRUE}, plots showing how
#'                      the model performs are shown
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
#' @references Rheault et al., (2024). goFlux: A user-friendly way to
#' calculate GHG fluxes yourself, regardless of user experience. \emph{Journal
#' of Open Source Software}, 9(96), 6393, https://doi.org/10.21105/joss.06393
#'
#' @export
#'
#' @examples
#' #' # IMPORTANT! This function makes use of the function graphics::identify()
#' # which is only supported on screen devices such as X11, windows and quartz.
#' # It is therefore essential to verify that your system options are compatible
#' # with this function before running it, to avoid errors. Here is an example
#' # of how to modify your system options for graphics device:
#' \dontrun{
#' default.device <- getOption("device") # save default option
#' options(device = "X11") # change system option to device = "X11"
#' options(device = default.device) # revert back to default option }
#'
#' # Loading data
#' load(data_example_1)
#' mydata$Etime <- as.numeric(mydata$Etime)
#'
#' Loading auxfile table
#' myauxfile = read.csv("data/myauxfile.csv")
#' myauxfile$start.time <- as.POSIXct(myauxfile$start.time, tz = 'UTC', format="%d/%m/%Y %H:%M")
#'
#' CH4_flux.manual <- clickflux(dataframe = mydata_all, myauxfile = myauxfile, shoulder = 0, gastype = "CH4dry_ppb",
#' plot.lim = c(1800,max(mydata_all$CH4dry_ppb)), fluxSeparation = T, displayPlots = T)
#' clickflux <- function(dataframe, myauxfile,
#'                       shoulder = 0,
#'                       gastype, plot.lim,
#'                       fluxSeparation,
#'                       displayPlots)
#'
#'
clickflux <- function(dataframe, myauxfile,
                      shoulder = 0,
                      gastype, plot.lim,
                      fluxSeparation,
                      displayPlots){

  # ----------------------- Check arguments -------------------------###

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

  ## Check plot.lim ####
  if(!is.numeric(plot.lim) | length(plot.lim) != 2){
    stop("'plot.lim' must be numeric and of length 2")}

  ## Check fluxSeparation ####
  if(missing(fluxSeparation)) {
    fluxSeparation = FALSE
    warning("in automaticflux(), 'fluxSeparation' was not provided and automatically set to FALSE")
  }
  if(fluxSeparation != TRUE & fluxSeparation != FALSE){
    stop("'fluxSeparation' must be TRUE or FALSE")}

  ## Check displayPlots ####
  if(missing(displayPlots)) {
    displayPlots = FALSE
    warning("in automaticflux(), 'displayPlots' was not provided and automatically set to FALSE")
  }
  if(displayPlots != TRUE & displayPlots != FALSE){
    stop("'displayPlots' must be TRUE or FALSE")}



  # ----------------------- Function starts here -------------------------###


  # making sure we work with the auxfile corresponding to measurements
  myauxfile <- myauxfile[myauxfile$UniqueID %in% unique(dataframe$UniqueID),]

  n_incubations <- length(unique(dataframe$UniqueID))

  mydata_ow <- obs.win(inputfile = dataframe, auxfile = myauxfile, shoulder = shoulder)

  message("... Please click on start and end points of VALID OBSERVATIONS")
  mydata_manID <- click.peak2(mydata_ow, seq = seq(1,n_incubations),
                              gastype = gastype,
                              plot.lim = plot.lim)

  # Calculate fluxes
  flux_manID <- goFlux(mydata_manID, gastype)

  # Use best.flux to select the best flux estimates (LM or HM)
  # based on a list of criteria
  criteria <- c("g.factor", "kappa", "MDF", "R2", "SE.rel")

  best.flux_manID <- best.flux(flux_manID, criteria)

  if(!fluxSeparation){

    if(displayPlots){
      p <- flux.plot(
        flux.results = best.flux_manID, dataframe = mydata_manID,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
        plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
      print(p)
    }

    return(best.flux_manID)
  } else { # in  that case we proceed with the separation between diffusion and ebullition
    message("... Please now select DIFFUSION patterns only")

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

    best.flux_manID <- lapply(seq_along(mydata_ow_corr), flux.separator.loop,
                              list_of_dataframes = mydata_ow_corr, gastype = gastype, auxfile = myauxfile,
                              criteria = criteria, force.separation = F,
                              mybest.flux = best.flux_manID)%>%
      map_df(., ~as.data.frame(.x))

    if(displayPlots){
      p <- separated.flux.plot(
        flux.results = best.flux_manID, dataframe = mydata_manID,
        gastype = gastype, quality.check = TRUE,
        plot.legend = c("SD"),
        plot.display = c("nb.obs"))
      print(p)
    }

    return(best.flux_manID)
  }
}





