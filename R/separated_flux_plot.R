
#' Title
#'
#' @param flux.results blabla
#' @param dataframe blabla
#' @param gastype blabla
#' @param shoulder blabla
#' @param plot.legend blabla
#' @param plot.display blabla
#' @param quality.check blabla
#' @param flux.unit blabla
#' @param flux.term.unit blabla
#' @param p.val.disp blabla
#' @param side blabla
#'
#' @return blabla
#'
#' @examples
#' blabla
#'
#' @export
separated_flux_plot <- function(flux.results, dataframe, gastype, shoulder = 30,
                      plot.legend = c("SD"),
                      plot.display = c("prec"),
                      quality.check = TRUE, flux.unit = NULL,
                      flux.term.unit = NULL,
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
  }

  ## Check plot.legend ####
  plot.legend.all <- c("SD")
  if(!is.null(plot.legend)){
    if(!is.character(plot.legend)){
      stop("'plot.legend' must be of class character")
    } else if(!any(grepl(paste(paste("\\<", plot.legend.all, "\\>", sep = ""),
                               collapse = "|"), plot.legend))){
      stop("if 'plot.legend' is not NULL, it must contain at least one of the following: 'SD'")
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
      }
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



  # Hutchinson and Mosier model
  HMmod <- function(Ci, C0, k, x){
    Ci + (C0 - Ci) * exp(-k * x)
  }


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

    ### total.flux, diffusion.flux and ebullition.flux are always in the legend
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


    ## plot.display ####
    if(!is.null(plot.display)){

      ### NEW PLOT LIMITS with nb.obs, flux term and prec ####
      if(any(grepl(paste(c("\\<nb.obs\\>", "\\<flux.term\\>", "\\<prec\\>"),
                         collapse = "|"), plot.display))){
        display.length <- length(grep(paste(c(
          "\\<nb.obs\\>", "\\<flux.term\\>", "\\<prec\\>"),
          collapse = "|"), plot.display))
        if(display.length > 2) ymin <- ymin - ydiff*min(seq.y)*1.8
        if(display.length <= 2) ymin <- ymin - ydiff*min(seq.y)*0.9
      }
      ### nb.obs ####
      if(any(grepl("\\<nb.obs\\>", plot.display))){
        # position
        nb.obs.ord <- which(grep(paste(c("\\<nb.obs\\>",
                                         "\\<flux.term\\>", "\\<prec\\>"), collapse = "|"),
                                 plot.display, value = T) == "nb.obs")
        if(nb.obs.ord == 1 | nb.obs.ord == 3) nb.obs.x <- 6
        if(nb.obs.ord == 2 | nb.obs.ord == 4) nb.obs.x <- 3
        if(nb.obs.ord <= 2) nb.obs.y <- -0.8
        if(nb.obs.ord > 2) nb.obs.y <- 1.6
        # value
        nb.obs <- round(unique(data_corr[[f]]$nb.obs), 0)
        # nb.obs.display
        nb.obs.display <- annotate(
          "text", x = seq(xmin, xmax, length.out=9)[nb.obs.x], colour = "black",
          y = ymin - ydiff*max(seq.y)*nb.obs.y/2, hjust = 0,
          label = paste(nb.obs, "~'data points used for diffusion'"), parse = TRUE, size = 3.2)
      }

      ### flux.term ####
      if(any(grepl("\\<flux.term\\>", plot.display))){
        # position
        flux.term.ord <- which(grep(paste(c("\\<nb.obs\\>",
                                            "\\<flux.term\\>", "\\<prec\\>"), collapse = "|"),
                                    plot.display, value = T) == "flux.term")
        if(flux.term.ord == 1 | flux.term.ord == 3) flux.term.x <- 6
        if(flux.term.ord == 2 | flux.term.ord == 4) flux.term.x <- 3
        if(flux.term.ord <= 2) flux.term.y <- -0.8
        if(flux.term.ord > 2) flux.term.y <- 1.6
        # value
        flux.term <- round(unique(data_corr[[f]]$flux.term), 1)
        # flux.term.display
        flux.term.display <- annotate(
          "text", x = seq(xmin, xmax, length.out=9)[flux.term.x], colour = "black",
          y = ymin - ydiff*max(seq.y)*flux.term.y/2, hjust = 0, parse = TRUE,
          label = paste("'flux.term ='~", flux.term, "~", flux.term.unit), size = 3.2)
      }
      ### prec ####
      if(any(grepl("\\<prec\\>", plot.display))){
        # position
        prec.ord <- which(grep(paste(c("\\<nb.obs\\>",
                                       "\\<flux.term\\>", "\\<prec\\>"), collapse = "|"),
                               plot.display, value = T) == "prec")
        if(prec.ord == 1 | prec.ord == 3) prec.x <- 6
        if(prec.ord == 2 | prec.ord == 4) prec.x <- 3
        if(prec.ord <= 2) prec.y <- -0.8
        if(prec.ord > 2) prec.y <- 1.6
        # value
        prec <- unique(data_corr[[f]]$prec)
        # prec.display
        prec.display <- annotate(
          "text", x = seq(xmin, xmax, length.out=9)[prec.x], colour = "black",
          y = ymin - ydiff*max(seq.y)*prec.y/2, hjust = 0,
          label = paste("'prec ='~", prec, "~", gas.unit), parse = TRUE, size = 3.2)
      }
    }

    ## Extract quality check ####
    if(quality.check == TRUE){
      # NEW PLOT LIMITS
      ymin <- ymin - ydiff*min(seq.y)*1.8

      # value
      quality <- unique(data_split[[f]]$quality.check)
      # quality.check.display
      quality.check.display <- annotate(
        "text", x = seq(xmin, xmax, length.out=9)[2], colour = "black",
        y = ymin - ydiff*min(seq.y)*.9, hjust = 0, parse = TRUE, size = 3.2,
        label = paste("'Quality check for diffusion:'~", paste("'", quality, "'")))
    }

    # Content of plot
    Etime <- data_split[[f]]$Etime
    gas_meas <- Reduce("c", data_split[[f]][, gastype])
    flag <- data_split[[f]]$flag
    plot_data <- cbind.data.frame(gas_meas, Etime, flag)
    plot_data$flag_diff <- F
    plot_data$flag_diff[which(data_split[[f]]$POSIX.time > unique(data_split[[f]]$start_diffusion) &
                                data_split[[f]]$POSIX.time < unique(data_split[[f]]$start_diffusion) + unique(data_split[[f]]$obs.length_diffusion))] <- T

    LM.slope <- unique(data_corr[[f]]$LM.slope)
    # reporting C0 to t = 0
    LM.C0 <- unique(data_corr[[f]]$LM.C0) - LM.slope*(as.numeric(unique(data_split[[f]]$start_diffusion)) - min(as.numeric(data_split[[f]]$POSIX.time)))
    UniqueID <- unique(data_corr[[f]]$UniqueID)
    HM_mod <- data_split[[f]]$HM_mod

    # Draw plot ####
    plot <- ggplot(plot_data, aes(x = Etime)) +
      geom_point(aes(y = gas_meas, col = as.factor(flag_diff))) +
      scale_color_manual(values = c("darkgrey", "black"), guide = "none") +


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


    if(unique(data_corr[[f]]$model) == "LM"){
      plot <- plot +
        # Linear model
        geom_abline(slope = LM.slope, intercept = LM.C0,
                    linewidth = 1, col = "red")
    } else {
      plot <- plot +
        # Hutchinson and Mosier
        geom_line(aes(y = HM_mod), linewidth = 1, col = "red")
    }

    return(plot)
  })
}
