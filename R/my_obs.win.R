

my_obs.win <- function (inputfile, auxfile = NULL, gastype = "CO2dry_ppm", 
                     obs.length = NULL, shoulder = 120) {
  if (missing(inputfile)) 
    stop("'inputfile' is required")
  if (!is.null(inputfile) & !is.data.frame(inputfile)) 
    stop("'inputfile' must be of class data.frame")
  if (!any(grepl(gastype, c("CO2dry_ppm", "COdry_ppb", "CH4dry_ppb", 
                            "N2Odry_ppb", "NH3dry_ppb", "H2O_ppm")))) {
    stop("'gastype' must be of class character and one of the following: 'CO2dry_ppm', 'COdry_ppm', 'CH4dry_ppb', 'N2Odry_ppb', 'NH3dry_ppb' or 'H2O_ppm'")
  }
  if (!is.null(auxfile) & !is.data.frame(auxfile)) 
    stop("'auxfile' must be of class data.frame")
  if (is.null(shoulder)) 
    stop("'shoulder' is required")
  else {
    if (!is.numeric(shoulder)) 
      stop("'shoulder' must be of class numeric")
    else {
      if (shoulder < 0) 
        stop("'shoulder' cannot be a negative value")
    }
  }
  if (is.null(auxfile)) {
    if (!any(grepl("\\<UniqueID\\>", names(inputfile))) & 
        !any(grepl("\\<chamID\\>", names(inputfile)))) {
      stop("'UniqueID' is required and was not found in 'inputfile'. Alternatively, provide chamID in 'inputfile'.")
    }
  }
  else {
    if (!any(grepl("\\<UniqueID\\>", names(auxfile))) & !any(grepl("\\<chamID\\>", 
                                                                   names(auxfile)))) {
      stop("'UniqueID' is required and was not found in 'auxfile'. Alternatively, provide chamID in 'auxfile'.")
    }
  }
  if (is.null(auxfile)) {
    if (!any(grepl("\\<UniqueID\\>", names(inputfile))) & 
        any(grepl("\\<chamID\\>", names(inputfile)))) {
      if (!any(grepl("\\<DATE\\>", names(inputfile)))) {
        stop("The column DATE in required in 'inputfile' to create a UniqueID from chamID.")
      }
    }
  }
  else {
    if (!any(grepl("\\<UniqueID\\>", names(auxfile))) & any(grepl("\\<chamID\\>", 
                                                                  names(auxfile)))) {
      if (!any(grepl("\\<DATE\\>", names(inputfile)))) {
        stop("The column DATE in required in 'inputfile' to create a UniqueID from chamID.")
      }
    }
  }
  if (is.null(auxfile)) {
    if (!any(grepl("\\<start.time\\>", names(inputfile)))) {
      stop("'start.time' is required and was not found in 'inputfile'")
    }
    if (any(grepl("\\<start.time\\>", names(inputfile))) & 
        !is.POSIXct(inputfile$start.time)) {
      stop("'start.time' in 'inputfile' must be of class POSIXct")
    }
    if (any(grepl("\\<start.time\\>", names(inputfile))) & 
        is.POSIXct(inputfile$start.time) & attr(inputfile$start.time, 
                                                "tzone") != attr(inputfile$POSIX.time, "tzone")) {
      stop("'start.time' in 'inputfile' must be in the same time zone as 'POSIX.time'")
    }
  }
  else {
    if (!any(grepl("\\<start.time\\>", names(auxfile)))) {
      stop("'start.time' is required and was not found in 'auxfile'")
    }
    if (any(grepl("\\<start.time\\>", names(auxfile))) & 
        !is.POSIXct(auxfile$start.time)) {
      stop("'start.time' in 'auxfile' must be of class POSIXct")
    }
    if (any(grepl("\\<start.time\\>", names(auxfile))) & 
        is.POSIXct(auxfile$start.time) & attr(auxfile$start.time, 
                                              "tzone") != attr(inputfile$POSIX.time, "tzone")) {
      stop("'start.time' in 'auxfile' must be in the same time zone as 'POSIX.time' in 'inputfile'")
    }
  }
  if (!is.null(obs.length) & !is.numeric(obs.length)) {
    stop("'obs.length' must be of class numeric")
  }
  if (is.null(obs.length)) {
    if (!any(grepl("\\<obs.length\\>", c(names(inputfile), 
                                         names(auxfile))))) {
      if (!is.null(auxfile)) {
        if (length(grep(paste(c("\\<end.time\\>", "\\<cham.open\\>"), 
                              collapse = "|"), c(names(inputfile), names(auxfile)))) < 
            1) {
          stop("'obs.length' missing. 'inputfile' or 'auxfile' must contain alternative arguments to calculate 'obs.length': 'start.time' and 'cham.open' or 'end.time'.")
        }
      }
      else {
        if (length(grep(paste(c("\\<end.time\\>", "\\<cham.open\\>"), 
                              collapse = "|"), names(inputfile))) < 1) {
          stop("'obs.length' missing. 'inputfile' must contain alternative arguments to calculate 'obs.length': start.time' and 'cham.open' or 'end.time'.")
        }
      }
    }
  }
  if (any(grepl("\\<end.time\\>", names(inputfile)))) {
    if (!is.POSIXct(inputfile$end.time)) {
      stop("'end.time' in 'inputfile' must be of class POSIXct")
    }
    else if (attr(inputfile$end.time, "tzone") != attr(inputfile$POSIX.time, 
                                                       "tzone")) {
      stop("'end.time' in 'inputfile' must be in the same time zone as 'POSIX.time'")
    }
  }
  if (any(grepl("\\<end.time\\>", names(auxfile)))) {
    if (!is.POSIXct(auxfile$end.time)) {
      stop("'end.time' in 'auxfile' must be of class POSIXct")
    }
    else if (attr(auxfile$end.time, "tzone") != attr(inputfile$POSIX.time, 
                                                     "tzone")) {
      stop("'end.time' in 'auxfile' must be in the same time zone as 'POSIX.time' in 'inputfile'")
    }
  }
  if (any(grepl("\\<cham.open\\>", names(inputfile)))) {
    if (!is.POSIXct(inputfile$cham.open)) {
      stop("'cham.open' in 'inputfile' must be of class POSIXct")
    }
    else if (attr(inputfile$cham.open, "tzone") != attr(inputfile$POSIX.time, 
                                                        "tzone")) {
      stop("'cham.open' in 'inputfile' must be in the same time zone as 'POSIX.time'")
    }
  }
  if (any(grepl("\\<cham.open\\>", names(auxfile)))) {
    if (!is.POSIXct(auxfile$cham.open)) {
      stop("'cham.open' in 'auxfile' must be of class POSIXct")
    }
    else if (attr(auxfile$cham.open, "tzone") != attr(inputfile$POSIX.time, 
                                                      "tzone")) {
      stop("'cham.open' in 'auxfile' must be in the same time zone as 'POSIX.time' in 'inputfile'")
    }
  }
  
  
  
  POSIX.time <- chamID <- start.time <- UniqueID <- Etime <- flag <- DATE <- CO2dry_ppm <- COdry_ppb <- CH4dry_ppb <- N2Odry_ppb <- NH3dry_ppb <- H2O_ppm <- cham.open <- end.time <- Tcham <- Pcham <- NULL
  inputfile <- inputfile %>% mutate(POSIX.time = as.POSIXct(round(POSIX.time, 
                                                                  "secs")))
  if (is.null(auxfile)) {
    if (any(grepl("\\<chamID\\>", names(inputfile)))) {
      inputfile <- inputfile %>% mutate(UniqueID = paste(chamID, 
                                                         DATE, sep = "_")) %>% mutate(start.time = as.POSIXct(round(start.time, 
                                                                                                                    "secs")))
    }
    aux.data <- inputfile %>% select(UniqueID, start.time) %>% 
      distinct()
  }
  else {
    if (any(grepl("\\<chamID\\>", names(auxfile)))) {
      auxfile <- auxfile %>% mutate(UniqueID = paste(chamID, 
                                                     DATE, sep = "_"))
    }
    
    
    # here
    aux.data <- auxfile %>% select(UniqueID, start.time) %>% 
      distinct()
  
    
    }
  if (any(grepl("\\<cham.open\\>", names(inputfile)))) {
    inputfile <- inputfile %>% mutate(end.time = cham.open)
  }
  if (any(grepl("\\<cham.open\\>", names(auxfile)))) {
    auxfile <- auxfile %>% mutate(end.time = cham.open)
  }
  if (!is.null(obs.length)) {
    aux.data <- aux.data %>% mutate(obs.length = obs.length)
  }
  else {
    if (!is.null(auxfile) & any(grepl("\\<obs.length\\>", names(auxfile)))) {
      
      
      # here
      aux.data <- aux.data %>% left_join(auxfile %>% select(UniqueID, 
                                                            obs.length) %>% distinct(), by = "UniqueID")
      
      
    }
    else if (!is.null(auxfile) & any(grepl("\\<end.time\\>", 
                                           names(auxfile)))) {
      aux.data <- aux.data %>% left_join(auxfile %>% select(UniqueID, 
                                                            end.time) %>% distinct(), by = "UniqueID") %>% 
        group_by(UniqueID) %>% mutate(obs.length = as.numeric(end.time - 
                                                                start.time, units = "secs")) %>% ungroup() %>% 
        select(-end.time)
    }
    else if (any(grepl("\\<obs.length\\>", names(inputfile)))) {
      aux.data <- aux.data %>% left_join(inputfile %>% 
                                           select(UniqueID, obs.length) %>% distinct(), 
                                         by = "UniqueID")
    }
    else if (any(grepl("\\<end.time\\>", names(inputfile)))) {
      aux.data <- aux.data %>% left_join(inputfile %>% 
                                           select(UniqueID, end.time) %>% distinct(), by = "UniqueID") %>% 
        group_by(UniqueID) %>% mutate(obs.length = as.numeric(end.time - 
                                                                start.time, units = "secs")) %>% ungroup() %>% 
        select(-end.time)
    }
  }
  
  
  time_range <- aux.data %>% group_by(UniqueID) %>% reframe(time_min = start.time - 
                                                              shoulder, time_max = start.time + obs.length + shoulder, 
                                                            obs.length = obs.length, start.time = start.time)
  time_filter.ls <- list()
  for (i in 1:nrow(time_range)) {
    time_filter.ls[[i]] <- cbind.data.frame(UniqueID = time_range$UniqueID[[i]], 
                                            start.time = time_range$start.time[[i]], obs.length = time_range$obs.length[[i]], 
                                            POSIX.time = seq(from = time_range$time_min[i], to = time_range$time_max[i], 
                                                             by = "sec"))
  }
  time_filter <- map_df(time_filter.ls, ~as.data.frame(.x)) %>% 
    distinct()
  if (any(grepl("\\<UniqueID\\>", names(inputfile)))) {
    inputfile <- inputfile %>% select(-UniqueID)
  }
  if (any(grepl("\\<start.time\\>", names(inputfile)))) {
    inputfile <- inputfile %>% select(-start.time)
  }
  if (any(grepl("\\<obs.length\\>", names(inputfile)))) {
    inputfile <- inputfile %>% select(-obs.length)
  }
  data.filter <- inputfile %>% right_join(time_filter, relationship = "many-to-many", 
                                          by = "POSIX.time")# %>% drop_na(matches(gastype))
  
  
  if (any(grepl("\\<CO2dry_ppm\\>", names(data.filter)))) {
    data.filter[["CO2dry_ppm"]] <- approx(inputfile$POSIX.time, inputfile[["CO2dry_ppm"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["CO2_prec"]] <- first(inputfile[["CO2_prec"]])
  }
  if (any(grepl("\\<COdry_ppb\\>", names(data.filter)))) {
    data.filter[["COdry_ppb"]] <- approx(inputfile$POSIX.time, inputfile[["COdry_ppb"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["CO_prec"]] <- first(inputfile[["CO_prec"]])
  }
  if (any(grepl("\\<CH4dry_ppb\\>", names(data.filter)))) {
    data.filter[["CH4dry_ppb"]] <- approx(inputfile$POSIX.time, inputfile[["CH4dry_ppb"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["CH4_prec"]] <- first(inputfile[["CH4_prec"]])
  }
  if (any(grepl("\\<N2Odry_ppb\\>", names(data.filter)))) {
    data.filter[["N2Odry_ppb"]] <- approx(inputfile$POSIX.time, inputfile[["N2Odry_ppb"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["N2O_prec"]] <- first(inputfile[["N2O_prec"]])
  }
  if (any(grepl("\\<NH3dry_ppb\\>", names(data.filter)))) {
    data.filter[["NH3dry_ppb"]] <- approx(inputfile$POSIX.time, inputfile[["NH3dry_ppb"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["NH3_prec"]] <- first(inputfile[["NH3_prec"]])
  }
  if (any(grepl("\\<H2O_ppm\\>", names(data.filter)))) {
    data.filter[["H2O_ppm"]] <- approx(inputfile$POSIX.time, inputfile[["H2O_ppm"]], xout = data.filter$POSIX.time, rule = 2)$y
    data.filter[["H2O_prec"]] <- first(inputfile[["H2O_prec"]])
  }
  
  data.filter[["Etime"]] <- as.numeric(data.filter$POSIX.time) - min(as.numeric(data.filter$POSIX.time))
  
  
  if (!is.null(auxfile)) {
    if (any(grepl("\\<obs.length\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(-obs.length)
    }
    if (any(grepl("\\<start.time\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(-start.time)
    }
    if (any(grepl("\\<Etime\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(Etime))
    }
    if (any(grepl("\\<flag\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(flag))
    }
    if (any(grepl("\\<DATE\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(DATE))
    }
    if (any(grepl("\\<CO2dry_ppm\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(CO2dry_ppm))
    }
    if (any(grepl("\\<COdry_ppb\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(CO2dry_ppm))
    }
    if (any(grepl("\\<CH4dry_ppb\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(CH4dry_ppb))
    }
    if (any(grepl("\\<N2Odry_ppb\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(N2Odry_ppb))
    }
    if (any(grepl("\\<NH3dry_ppb\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(N2Odry_ppb))
    }
    if (any(grepl("\\<H2O_ppm\\>", names(auxfile)))) {
      auxfile <- auxfile %>% select(!c(H2O_ppm))
    }
    if (any(grepl("\\<Tcham\\>", names(auxfile)))) {
      Tcham <- first(na.omit(auxfile$Tcham))
      auxfile <- auxfile %>% select(!c(Tcham))
      data.filter$Tcham <- Tcham
    }
    if (any(grepl("\\<Pcham\\>", names(auxfile)))) {
      Pcham <- first(na.omit(auxfile$Pcham))
      auxfile <- auxfile %>% select(!c(Pcham))
      data.filter$Pcham <- Pcham
    }
    auxfile <- auxfile %>% distinct()
    if (any(grepl("\\<POSIX.time\\>", names(auxfile)))) {
      data.filter <- data.filter %>% left_join(auxfile, 
                                               by = c("UniqueID", "POSIX.time"))
    } else {
      data.filter <- data.filter %>% full_join(auxfile, 
                                               by = "UniqueID")
    }
  }
  
  data.filter <- data.filter[order(data.filter$POSIX.time),]
  
  flux.unique <- data.filter %>% group_split(UniqueID) %>% 
    as.list()
  if (length(flux.unique) > 20) {
    message("WARNING! Do not loop through more than 20 measurements at a time to avoid mistakes.", 
            "\nYou have ", length(flux.unique), " measurements in your dataset.", 
            "\nYou should split the next step into at least ", 
            round(length(flux.unique)/20), " loops.")
  }
  return(flux.unique)
}
