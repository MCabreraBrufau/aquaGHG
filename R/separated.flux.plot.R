#' Plots for quality checking of GHG flux measurements
#'
#' Returns a list of plots, drawn from flux results (output from the functions
#' \code{\link[goFlux]{goFlux}} and \code{\link[goFlux]{best.flux}}).
#' The plots are customizable.
#'
#' @param flux.results a data.frame; output from the function
#'                     \code{\link[goFlux]{best.flux}}
#' @param dataframe a data.frame containing gas measurements (see \code{gastype}
#'                  below) and the following columns: \code{UniqueID}, \code{Etime}
#'                  and \code{flag} (same \code{dataframe} as used with the function
#'                  \code{\link[goFlux]{goFlux}}). \code{chamID} may be used
#'                  instead of \code{UniqueID}.
#' @param gastype character string; specifies which column was used for the
#'                flux calculations. Must be one of the following: "CO2dry_ppm",
#'                "COdry_ppb", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb" or "H2O_ppm".
#' @param shoulder numerical value; time before and after measurement in observation
#'                 window (seconds). Default is 30 seconds.
#' @param plot.legend character vector; specifies which parameters should be
#'                    displayed in a legend above each plot. "flux" is always
#'                    displayed. A maximum of 6 parameters can be displayed in the
#'                    legend (including "flux"). Chose up to five extra parameters
#'                    from the following: "MAE", "RMSE", "AICc", "SErel", "SE",
#'                    "r2", "LM.p.val", "HM.k", "k.max", "k.ratio" and "g.factor".
#'                    Default is \code{plot.legend = c("MAE", "AICc", "k.ratio", "g.factor")}.
#' @param plot.display character vector; specifies which parameters should be
#'                     displayed on the plot. Choose from the following: "C0",
#'                     "Ci", "cham.close", "cham.open", "crop", "MDF", "nb.obs",
#'                     "flux.term" and "prec". Default is
#'                     \code{plot.display = c("MDF", "prec")}.
#' @param quality.check logical; if \code{quality.check = TRUE}, the column
#'                      \code{quality.check} (output from the function
#'                      \code{\link[goFlux]{best.flux}}) is displayed
#'                      below the plot.
#' @param flux.unit character string; flux units to be displayed on the plots.
#'                  By default, the units are
#'                  \ifelse{html}{\out{µmol m<sup>-2</sup>s<sup>-1</sup>}}{\eqn{µmol m^{-2}s^{-1}}{ASCII}}
#'                  (if initial concentration is ppm, e.g. CO2dry_ppm) and
#'                  \ifelse{html}{\out{nmol m<sup>-2</sup>s<sup>-1</sup>}}{\eqn{nmol m^{-2}s^{-1}}{ASCII}}
#'                  (if initial concentration is ppb, e.g. CH4dry_ppb).
#'                  \cr \cr
#'                  For example, one may want to use
#'                  \ifelse{html}{\out{nmol kg<sup>-1</sup>h<sup>-1</sup>}}{\eqn{µmol kg^{-1}h^{-1}}{ASCII}}
#'                  for incubated soil samples. In such a case, write
#'                  \code{flux.unit = "nmol~kg^-1*h^-1"}.
#' @param flux.term.unit character string; units for the flux term to be displayed
#'                       on the plots. By default, the units are
#'                       \ifelse{html}{\out{mol m<sup>-2</sup>}}{\eqn{mol m^{-2}}{ASCII}}:
#'                       \code{flux.term.unit = "mol~m^-2"}
#' @param best.model logical; if \code{best.model = TRUE}, display a star sign
#'                   next to the best model selected by the function
#'                   \code{\link[goFlux]{best.flux}}.
#' @param p.val.disp character string; indicates how the \emph{p-value} should
#'                   be displayed in the legend above the plot. Choose one of
#'                   the following: "star", "round", "value".
#' @param side character string; choose a side to display
#'             \ifelse{html}{\out{C<sub>0</sub>}}{\eqn{mol C[0]}{ASCII}} and
#'             \ifelse{html}{\out{C<sub>i</sub>}}{\eqn{mol C[i]}{ASCII}} values.
#'             By default, they are displayed on the left side of the plot.
#'
#' @details
#' In \code{flux.results}, one may choose to use the output from the
#' \code{\link[goFlux]{goFlux}} function instead. However, in that case,
#' any element produced by the function \code{\link[goFlux]{best.flux}}
#' cannot be displayed on the plots.
#'
#' In \code{plot.legend}, one may choose to display up to five additional parameters
#' in a legend above the plots. Some parameters are displayed for both the linear
#' model (\code{\link[goFlux]{LM.flux}}) and the non-linear model
#' (\code{\link[goFlux]{HM.flux}}): Mean Absolute Error (\code{MAE}),
#' Root Mean Square Error (\code{RMSE}), Aikaike's Information Criterion
#' corrected for small sample size (\code{LM.AICc}), Standard Error (\code{SE}),
#' relative Standard Error (\code{SErel}), and coefficient of determination
#' (\code{r2}). The \emph{p-value} (\code{LM.p.val}) is displayed for the linear
#' model only. The kappa (\code{HM.k}), kappa-max
#' (\code{\link[goFlux]{k.max}}), kappa ratio (\code{k.ratio}) and
#' g-factor (\code{\link[goFlux]{g.factor}}) are displayed for the
#' Hutchinson and Mosier model only. One may choose to display no additional
#' parameter with \code{plot.legend = NULL}.
#'
#' In \code{plot.display}, one may chose to display some parameters on the plot:
#' The initial gas concentration (\code{C0}) for both models, the assumed
#' concentration of constant gas source below the surface (\code{Ci}) calculated
#' from the Hutchinson and Mosier model, the number of observations
#' (\code{nb.obs}) flagged, the Minimal Detectable Flux
#' (\code{\link[goFlux]{MDF}}), the flux term
#' (\code{\link[goFlux]{flux.term}}), the instrument precision
#' (\code{prec}), the chamber closure (\code{cham.close}) and opening
#' (\code{cham.open}) (indicated with a green star), and the data points between
#' chamber closure and opening that have been removed (\code{crop}) (indicated
#' in light red). For manual chamber measurements, because there is no automatic
#' chamber closure and opening, no green stars can be displayed. In addition,
#' \code{crop} is only relevant if data points have been removed with the
#' function \code{crop.meas()} (this function is not available yet).
#' One may choose to display none of these parameters with \code{plot.display = NULL}.
#' The order in which \code{prec}, \code{flux.term}, \code{MDF} and \code{nb.obs}
#' are put in \code{plot.display = c()} decides the order in which they are
#' displayed at the bottom of the plot.
#'
#' In \code{flux.unit}, remember to multiply the flux results with an appropriate
#' factor to convert the results from a unit to another. If kilograms of soil
#' were used to calculate the fluxes (see the details section of the function
#' \code{\link[goFlux]{goFlux}}), the units would be
#' \ifelse{html}{\out{µmol kg<sup>-1</sup>s<sup>-1</sup>}}{\eqn{µmol kg^{-1}s^{-1}}{ASCII}}.
#' To convert the units to
#' \ifelse{html}{\out{µmol kg<sup>-1</sup>h<sup>-1</sup>}}{\eqn{µmol kg^{-1}h^{-1}}{ASCII}}
#' instead, one would need to multiply the flux results by 3600 to convert from
#' seconds to hours. To print non-ASCII characters use Unicode. For example, to
#' print the Greek letter "mu" (\eqn{µ}), use the Unicode \code{\\u00B5}:
#' \code{flux.unit = "\\u00B5mol~kg^-1*h^-1"}.
#'
#' In \code{p.val.disp}, if \code{p.val.disp = "star"}, the \emph{p-values} will
#' be displayed as star symbols (asterisks) as follows: ***, ** or * for
#' \emph{p-values} of p < 0.001, p < 0.01 and p < 0.05, respectively. If
#' \code{p.val.disp = "round"}, the \emph{p-values} are rounded to p < 0.001,
#' p < 0.01 and p < 0.05. If \code{p.val.disp = "value"}, the actual values are
#' displayed, rounded to two significant numbers.
#'
#' In \code{gastype}, the gas species listed are the ones for which this package
#' has been adapted. Please write to the maintainer of this package for
#' adaptation of additional gases.
#'
#' @return A list of plots, one per \code{UniqueID}, drawn from flux results (output
#' from the functions \code{\link[goFlux]{goFlux}} and
#' \code{\link[goFlux]{best.flux}}).
#'
#' @include goFlux-package.R
#'
#' @seealso See also the functions \code{\link[goFlux]{goFlux}},
#'          \code{\link[goFlux]{best.flux}} and
#'          \code{\link[goFlux]{flux2pdf}}
#'          for more information about usage.
#'
#' @seealso Look up the functions \code{\link[goFlux]{g.factor}},
#'          \code{\link[goFlux]{k.max}},
#'          \code{\link[goFlux]{MDF}},
#'          \code{\link[goFlux]{flux.term}},
#'          \code{\link[goFlux]{LM.flux}} and
#'          \code{\link[goFlux]{HM.flux}} for more information about
#'          these parameters.
#'
#' @examples
#' data(manID.UGGA)
#' CO2_flux <- goFlux(manID.UGGA, "CO2dry_ppm")
#' criteria <- c("MAE", "AICc", "g.factor", "MDF")
#' CO2_best <- best.flux(CO2_flux, criteria)
#' CO2_plots <- flux.plot(
#'   flux.results = CO2_best, dataframe = manID.UGGA,
#'   gastype = "CO2dry_ppm", quality.check = TRUE,
#'   plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
#'   plot.display = c("Ci", "C0", "MDF", "prec", "nb.obs", "flux.term"))
#'
#' @export
#'
separated.flux.plot <- function(flux.results, dataframe, gastype, shoulder = 30,
                      plot.legend = c("MAE", "AICc", "k.ratio", "g.factor"),
                      plot.display = c("MDF", "prec"),
                      quality.check = TRUE, flux.unit = NULL,
                      flux.term.unit = NULL, best.model = TRUE,
                      p.val.disp = "round", side = "left") {

  # Check arguments ####
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

  ### Etime and flag ####
  if(!any(grepl("\\<Etime\\>", names(dataframe)))) stop("'dataframe' must contain 'Etime'")
  if(any(grepl("\\<Etime\\>", names(dataframe))) & !is.numeric(dataframe$Etime)){
    stop("'Etime' in 'dataframe' must be of class numeric (or integer)")}

  if(!any(grepl("\\<flag\\>", names(dataframe)))) stop("'dataframe' must contain 'flag'")
  if(any(grepl("\\<flag\\>", names(dataframe))) & !is.numeric(dataframe$flag)){
    stop("'flag' in 'dataframe' must be of class numeric (or integer)")}

  ### POSIX.time ####
  if(!any(grepl("\\<POSIX.time\\>", names(dataframe)))) stop("'dataframe' must contain 'POSIX.time'")
  if(any(grepl("\\<POSIX.time\\>", names(dataframe))) & !is.POSIXct(dataframe$POSIX.time)){
    stop("'POSIX.time' in 'dataframe' must be of class POSIXct")}

  ### UniqueID (or chamID) ####
  if(!any(grepl(paste(c("\\<UniqueID\\>", "\\<chamID\\>"), collapse = "|"), names(dataframe)))){
    stop("'dataframe' must contain 'UniqueID'")}

  ## Check flux.results ####
  if(missing(flux.results)) stop("'flux.results' is required")
  if(!is.null(flux.results) & !is.data.frame(flux.results)){
    stop("'flux.results' must be of class data.frame")}

  ### UniqueID ####
  if(!any(grepl("\\<UniqueID\\>", names(flux.results)))){
    stop("'UniqueID' is required and was not found in 'flux.results'")}

  # Create UniqueID from chamID, if missing
  if(!any(grepl("\\<UniqueID\\>", names(dataframe)))){
    if(any(grepl("\\<chamID\\>", names(dataframe)))){
      dataframe <- dataframe %>% mutate(UniqueID = paste(chamID, DATE, sep = "_"))}
  }

  ### Is there any match between dataframe and flux.results?
  if(!any(unique(dataframe$UniqueID) %in% unique(flux.results$UniqueID))){
    stop("'UniqueID' in 'flux.results' has no match for 'UniqueID' in 'dataframe'")}

  ### Model results ####
  ### total.flux
  if(!any(grepl("\\<total.flux\\>", names(flux.results)))){
    stop("'total.flux' required in 'flux.results'")
  } else if(!is.numeric(flux.results$total.flux)){
    stop("'total.flux' in 'flux.results' must be of class numeric")}
  ### diffusion.flux
  if(!any(grepl("\\<diffusion.flux\\>", names(flux.results)))){
    stop("'diffusion.flux' required in 'flux.results'")
  } else if(!is.numeric(flux.results$diffusion.flux)){
    stop("'diffusion.flux' in 'flux.results' must be of class numeric")}
  ### ebullition.flux
  if(!any(grepl("\\<ebullition.flux\\>", names(flux.results)))){
    stop("'ebullition.flux' required in 'flux.results'")
  } else if(!is.numeric(flux.results$ebullition.flux)){
    stop("'ebullition.flux' in 'flux.results' must be of class numeric")}

  ## Check plot.legend ####
  plot.legend.all <- c("SD")
  if(!is.null(plot.legend)){
    if(!is.character(plot.legend)){
      stop("'plot.legend' must be of class character")
    } else if(length(plot.legend) > 5){
      stop("in 'plot.legend': A maximum of 5 additional parameters can be displayed above the plot.")
    } else if(!any(grepl(paste(paste("\\<", plot.legend.all, "\\>", sep = ""),
                               collapse = "|"), plot.legend))){
      stop("if 'plot.legend' is not NULL, it must contain at least one of the following: 'MAE', 'RMSE', 'AICc', 'SErel', 'SE', 'r2', 'LM.p.val', 'HM.k', 'k.max', 'k.ratio', 'g.factor', 'best.model'")
    }
    ### SD ####
    if(any(grepl("\\<SD\\>", plot.legend))){
      if(!any(grepl("\\<total.flux.SD\\>", names(flux.results)))){
        stop("'SD' selected in 'plot.legend', but 'total.flux.SD' missing in 'flux.results'")
      } else if(!is.numeric(flux.results$total.flux.SD)){
        stop("'total.flux.SD' in 'flux.results' must be of class numeric")}
      if(!any(grepl("\\<diffusion.flux.SD\\>", names(flux.results)))){
        stop("'SD' selected in 'plot.legend', but 'diffusion.SD' missing in 'flux.results'")
      } else if(!is.numeric(flux.results$diffusion.flux.SD)){
        stop("'diffusion.flux.SD' in 'flux.results' must be of class numeric")}
      if(!any(grepl("\\<ebullition.flux.SD\\>", names(flux.results)))){
        stop("'SD' selected in 'plot.legend', but 'ebullition.flux.SD' missing in 'flux.results'")
      } else if(!is.numeric(flux.results$ebullition.flux.SD)){
        stop("'ebullition.flux.SD' in 'flux.results' must be of class numeric")}
    }
  }
  ## Check plot.display ####
  plot.display.all <- c("nb.obs",
                        "prec", "flux.term")
  if(!is.null(plot.display)){
    if(!is.character(plot.display)){
      stop("'plot.display' must be of class character")
    } else if(!any(grepl(paste(paste("\\<", plot.display.all, "\\>", sep = ""),
                               collapse = "|"), plot.display))){
      stop("if 'plot.display' is not NULL, it must contain at least one of the following: 'Ci', 'C0', 'cham.close', 'cham.open', 'nb.obs', 'crop', 'prec', 'flux.term', 'MDF'")
    }
    ### nb.obs ####
    if(any(grepl("\\<nb.obs\\>", plot.display))){
      if(!any(grepl("\\<nb.obs\\>", names(flux.results)))){
        stop("'nb.obs' selected in 'plot.display', but missing in 'flux.results'")
      } else if(!is.numeric(flux.results$nb.obs)){
        stop("'nb.obs' in 'flux.results' must be of class numeric")}
    }
    ### prec ####
    if(any(grepl("\\<prec\\>", plot.display))){
      if(!any(grepl("\\<prec\\>", names(flux.results)))){
        stop("'prec' selected in 'plot.display', but missing in 'flux.results'")
      } else if(!is.numeric(flux.results$prec)){
        stop("'prec' in 'flux.results' must be of class numeric")}
    }
    ### flux.term ####
    if(any(grepl("\\<flux.term\\>", plot.display))){
      if(!any(grepl("\\<flux.term\\>", names(flux.results)))){
        stop("'flux.term' selected in 'plot.display', but missing in 'flux.results'")
      } else if(!is.numeric(flux.results$flux.term)){
        stop("'flux.term' in 'flux.results' must be of class numeric")}
    }
  }
  ## Check quality.check ####
  if(quality.check != TRUE & quality.check != FALSE){
    stop("'quality.check' must be TRUE or FALSE")}
  if(isTRUE(quality.check)){
    if(!any(grepl("\\<quality.check\\>", names(flux.results)))){
      stop("'quality.check' is TRUE, but missing in 'flux.results'")
    } else if(!is.character(flux.results$quality.check)){
      stop("'quality.check' in 'flux.results' must be of class character")}
  }
  ## flux.term and flux units ####
  if(!is.null(flux.unit)){
    if(!is.character(flux.unit)) stop("'flux.unit' must be of class character")}
  if(!is.null(flux.term.unit)){
    if(!is.character(flux.term.unit)) stop("'flux.term.unit' must be of class character")}

  ## side ####
  if(is.null(side)) {stop("'side' cannot be NULL")
  } else if(!is.character(side)){
    stop("'side' must be of class character")} else {
      if(!any(grepl(paste("\\<", side, "\\>", sep = ""), c("left", "right")))){
        stop("'side' must be one of the following: 'left' or 'right'")}
    }

  # Assign NULL to variables without binding ####
  UniqueID <- . <- flag <- start.Etime <-
    end.Etime <- Etime <- x <- y <- content <- color <- POSIX.time <-
    quality.check.display <-
    nb.obs.display <- prec.display <-
    flux.term.display <- DATE <-
    legend.flux <- legend.SD <- GASTYPE <- NULL


  # Function to find decimal places
  nb.decimal = function(x) {
    #length zero input
    if (length(x) == 0) return(numeric())

    #count decimals
    x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% as.numeric()
    x_int = floor(x) %>% abs() %>% nchar()
    x_nchr = x_nchr - 1 - x_int
    x_nchr[x_nchr < 0] = 0

    x_nchr
  }

  # seq.rep: A wrapper function that merges seq and rep, to repeat a sequence
  # or sequence a repeat.
  seq.rep <- function(from, by, n.rep, length.seq, rep.seq = F) {

    if (rep.seq == F) {

      sequence <- seq(from = from, by = by, length.out = length.seq)
      out.ls <- list()
      for (i in 1:length(sequence)) {
        out.ls[[i]] <- rep(sequence[i], n.rep) }
      unlist(out.ls)

    } else {

      repetition <- rep(from, n.rep)
      out.ls <- list()
      for (i in 1:length(repetition)) {
        out.ls[[i]] <- seq(from = repetition[i], by = by, length.out = length.seq) }
      unlist(out.ls)

    }
  }




  # -------------------------- FUNCTION STARTS -------------------------------####

  # Define gas units
  if(gastype == "CO2dry_ppm") gas.unit <- "ppm"
  if(gastype == "CH4dry_ppb") gas.unit <- "ppb"
  if(gastype == "N2Odry_ppb") gas.unit <- "ppb"
  if(gastype == "COdry_ppb") gas.unit <- "ppb"
  if(gastype == "NH3dry_ppb") gas.unit <- "ppb"
  if(gastype == "H2O_ppm") gas.unit <- "ppm"

  # Define y axis legend on plots
  if(gastype == "CO2dry_ppm") ylab <- ylab(expression(bold(CO["2"]*" dry (ppm)")))
  if(gastype == "CH4dry_ppb") ylab <- ylab(expression(CH["4"]*" dry (ppb)"))
  if(gastype == "N2Odry_ppb") ylab <- ylab(expression(N["2"]*"O dry (ppb)"))
  if(gastype == "COdry_ppb") ylab <- ylab(expression(CO*" dry (ppb)"))
  if(gastype == "NH3dry_ppb") ylab <- ylab(expression(NH["3"]*" dry (ppb)"))
  if(gastype == "H2O_ppm") ylab <- ylab(expression(H["2"]*"O (ppm)"))

  # Define flux units
  if(!is.null(flux.unit)) flux.unit <- flux.unit else {
    if(gastype == "CO2dry_ppm") flux.unit <- "\u00B5mol~m^-2*s^-1"
    if(gastype == "CH4dry_ppb") flux.unit <- "nmol~m^-2*s^-1"
    if(gastype == "N2Odry_ppb") flux.unit <- "nmol~m^-2*s^-1"
    if(gastype == "COdry_ppb") flux.unit <- "nmol~m^-2*s^-1"
    if(gastype == "NH3dry_ppb") flux.unit <- "nmol~m^-2*s^-1"
    if(gastype == "H2O_ppm") flux.unit <- "\u00B5mol~m^-2*s^-1"
  }

  # Define flux.term units
  if(!is.null(flux.term.unit)) flux.term.unit <- flux.term.unit else {
    flux.term.unit <- "mol~m^-2"
  }

  # Create a list of dataframe (by UniqueID)
  data_split <- dataframe %>%
    right_join(flux.results, by = c("UniqueID")) %>% group_by(UniqueID) %>%
    # Correct Etime for NAs
    # mutate(start.Etime = POSIX.time[which(Etime == 0)[1]],
    #        Etime = as.numeric(POSIX.time - start.Etime, units = "secs")) %>%
    # Calculate HM_mod
    mutate(HM_mod = HMmod(HM.Ci, HM.C0, HM.k, Etime)) %>%
    #select(!c(start.Etime)) %>%
    group_split()

  # Remove non-measurements (flag == 0)
  data_corr <- lapply(seq_along(data_split), function(f) {
    data_split[[f]] %>% filter(flag == 1) })

  # Loop through list of data frames (by UniqueID)
  pboptions(char = "=")
  plot_list <- pblapply(seq_along(data_split), function(f) {

    # Plot limits
    if(any(grepl("\\<nb.obs\\>", data_split[[f]]$quality.check))){

      xmax <- max(na.omit(data_split[[f]]$Etime))
      xmin <- min(na.omit(data_split[[f]]$Etime))
      xdiff <- xmax - xmin

      ymax <- max(na.omit(data_split[[f]][, gastype]))
      ymin <- min(na.omit(data_split[[f]][, gastype]))
      ydiff <- ymax - ymin

    } else {

      xmax <- max(na.omit(data_corr[[f]]$Etime)) + shoulder
      xmin <- -shoulder
      xdiff <- xmax - xmin

      ymax <- max(na.omit(data_corr[[f]][, gastype]))
      ymin <- min(na.omit(data_corr[[f]][, gastype]))
      ydiff <- ymax - ymin

      if(any(grepl("\\<crop\\>", plot.display))){
        if(nrow(dataframe %>% filter(flag == 2)) == 0){
          cham.close <- unique(na.omit(data_corr[[f]]$cham.close))
          cham.open <- unique(na.omit(data_corr[[f]]$cham.open))
          flag2 <- which(between(data_split[[f]]$POSIX.time, cham.close, cham.open)) %>%
            setdiff(which(data_split[[f]]$flag == 1))
          data_split[[f]] <- data_split[[f]] %>%
            mutate(flag = if_else(row_number() %in% flag2, 2, flag))

          xmax <- data_split[[f]] %>% filter(flag != 0) %>%
            select(Etime) %>% max(na.omit(.)) + shoulder
          xmin <- data_split[[f]] %>% filter(flag != 0) %>%
            select(Etime) %>% min(na.omit(.)) - shoulder
          xdiff <- xmax - xmin

          ymax <- data_split[[f]] %>% filter(flag != 0) %>%
            select(all_of(gastype)) %>% max(na.omit(.))
          ymin <- data_split[[f]] %>% filter(flag != 0) %>%
            select(all_of(gastype)) %>% min(na.omit(.))
          ydiff <- ymax - ymin
        }
      }
    }
    ## plot.legend ####

    # Variables decimals
    flux.dec <- nb.decimal(ifelse(
      nb.decimal(signif(unique(data_corr[[f]]$MDF), 1)) != 0,
      signif(unique(data_corr[[f]]$MDF), 2),
      round(unique(data_corr[[f]]$MDF), 1))) %>%
      ifelse(. > 2, . -1, .)

    gas.dec <- nb.decimal(unique(data_corr[[f]]$prec))

    ### LM and HM flux are always in the legend
    LM.flux <- round(unique(data_corr[[f]]$LM.flux), flux.dec)
    HM.flux <- round(unique(data_corr[[f]]$HM.flux), flux.dec)

    ### LM and HM flux are always in the legend
    total.flux <- round(unique(data_corr[[f]]$total.flux), flux.dec)
    diffusion.flux <- round(unique(data_corr[[f]]$diffusion.flux), flux.dec)
    ebullition.flux <- round(unique(data_corr[[f]]$ebullition.flux), flux.dec)

    legend.flux <- cbind.data.frame(
      content = c("Term", "Total", "Diffusion", "Ebullition", paste("'Flux units:'", "~", flux.unit),
                  "Flux", total.flux, diffusion.flux, ebullition.flux, ""))

    ### Legend length ####
    legend.length <- length(grep(paste(c(
      "\\<SD\\>"),
      collapse = "|"), plot.legend)) +2


    ### SD ####
    if(any(grepl("\\<SD\\>", plot.legend))){
      total.flux.SD <- round(unique(data_corr[[f]]$total.flux.SD), flux.dec)
      diffusion.flux.SD <- round(unique(data_corr[[f]]$diffusion.flux.SD), flux.dec)
      ebullition.flux.SD <- round(unique(data_corr[[f]]$ebullition.flux.SD), flux.dec)

      legend.SD <- cbind.data.frame(content = c("SD", total.flux.SD, diffusion.flux.SD, ebullition.flux.SD, ""))
    }

    # Legends' positions
    seq.x <- seq.rep(0.93, -0.13, 5, legend.length)
    seq.y <- seq.rep(0.28, -0.07, legend.length, 5, rep.seq = T)

    ## Merge legend data frames
    mod.legend <- rbind(legend.flux, legend.SD) %>%
      cbind.data.frame(
        color = rep(c("black", "blue", "red", "black", "black"), legend.length),
        x = xmax - xdiff*seq.x,
        y = ymax + ydiff*seq.y)

    ## Extract quality check ####
    if(quality.check == TRUE){
      # NEW PLOT LIMITS
      ymin <- ymin - ydiff*min(seq.y)*1.8

      # value
      quality <- unique(data_split[[f]]$quality.check)
      # quality.check.display
      quality.check.display <- annotate(
        "text", x = seq(xmin, xmax, length.out=9)[2], colour = "black",
        y = ymin - ydiff*min(seq.y)*0.9, hjust = 0, parse = TRUE, size = 3.2,
        label = paste("'Quality check for diffusion term:'~", paste("'", quality, "'")))
    }

    # Content of plot
    Etime <- data_split[[f]]$Etime
    gas_meas <- Reduce("c", data_split[[f]][, gastype])
    flag <- data_split[[f]]$flag
    plot_data <- cbind.data.frame(gas_meas, Etime, flag)
    plot_data$flag_diff <- F
    plot_data$flag_diff[plot_data$Etime < data_split[[f]]$obs.length_diffusion] <- T

    LM.slope <- unique(data_corr[[f]]$LM.slope)
    LM.C0 <- unique(data_corr[[f]]$LM.C0)
    UniqueID <- unique(data_corr[[f]]$UniqueID)
    HM_mod <- data_split[[f]]$HM_mod

    # Draw plot ####
    plot <- ggplot(plot_data, aes(x = Etime)) +
      geom_point(aes(y = gas_meas, col = as.factor(flag_diff))) +
      scale_color_manual(values = c("darkgrey", "black"), guide = "none") +


      # Linear model
      geom_abline(slope = LM.slope, intercept = LM.C0,
                  linewidth = 1, col = "blue") +

      # Hutchinson and Mosier
      geom_line(aes(y = HM_mod), linewidth = 1, col = "red") +



      # Add a legend with info on the two models
      new_scale_color() +
      geom_text(data = mod.legend, parse = T, size = 3.5,
                aes(x = x, y = y, label = content, hjust = 0, color = color)) +
      scale_color_manual(values = mod.legend$color, guide = "none") +

      # plot.display
      nb.obs.display +
      flux.term.display +
      prec.display +
      quality.check.display +

      # Make the plot pretty
      xlab("Time (sec)") + ylab +
      scale_x_continuous(breaks = seq(-60, max(Etime), 30),
                         minor_breaks = seq(-60, max(Etime)+60, 10)) +
      coord_cartesian(xlim = c(xmin + xdiff*0.05, xmax - xdiff*0.05),
                      ylim = c(ymin - ydiff*0.05, ymax + ydiff*max(seq.y))) +
      theme_bw() +
      ggtitle(UniqueID)+
      theme(axis.title.x = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(size = 10, face = "bold"))

    return(plot)
  })
}
