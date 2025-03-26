

find_first_linear_chunk <- function(dataframe, gastype, length.min){

  # dataframe <- mydata_all[mydata_all$UniqueID==unique(mydata_all$UniqueID)[3],]

  mydf <- data.frame(POSIX.time = dataframe$POSIX.time,
                     time = as.numeric(dataframe$Etime),
                     conc = dataframe[[gastype]])

  # identify possible local minimum in first half of incubation
  ind_min <- which.min(mydf$conc[mydf$time<max(mydf$time)/2])
  mydf <- mydf[mydf$time >= mydf$time[ind_min],]


  # initialization
  # by default, we select the entire incubation to process diffusion. This is refined below.
  first_linear_chunk <- data.frame(UniqueID = unique(dataframe$UniqueID),
                                   start = mydf$time[1],
                                   end = max(mydf$time))


  # a first Lm fit with a moving window to discard highly fluctuating parts of the measurements
  # moving window
  mywindow <- 100
  start = min(mydf$time)

  df.mw <- NULL
  for(t in unique(mydf$time[mydf$time>=start])){
    ind <- mydf$time >= t & mydf$time <= t+mywindow
    n = sum(ind)
    if(n>length.min){
      lm.tmp <- lm(formula = conc~time, data = mydf[ind,])
      df.mw <- rbind(df.mw, data.frame(t = t,
                                       n = n,
                                       slope = summary(lm.tmp)$coefficients[2,1],
                                       r2 = summary(lm.tmp)$adj.r.squared,
                                       pval = summary(lm.tmp)$coefficients[2,4]))
    }
  }
  # ggplot(df.mw, aes(t, slope))+geom_point()+theme_article()

  # we select only positive slopes
  df.mw <- df.mw[df.mw$slope>0,]
  if(dim(df.mw)[1]==0){
    warning("... moving-window slopes are all negative")
  }

  # ggplot(df.mw, aes(t, r2))+geom_point()+theme_article()

  # take a look at continuity of rising slope
  df.n <- NULL
  for(t in df.mw$t){
    n = sum(mydf$time<=t)
    if(n>length.min){
      df.n <- rbind(df.n, data.frame(t = t,
                                     n = n))
    }
  }


  # now we proceed with a linear fit over a growing window
  df.stats <- NULL
  for(t in df.n$t){
    ind <- mydf$time >=(df.n$t[1]-length.min) & mydf$time<=t
    n = sum(ind)
    if(n>length.min){
      lm.tmp <- lm(formula = conc~time, data = mydf[ind,])
      df.stats <- rbind(df.stats, data.frame(t = t,
                                             n = n,
                                             slope = summary(lm.tmp)$coefficients[2,1],
                                             r2 = summary(lm.tmp)$adj.r.squared,
                                             pval = summary(lm.tmp)$coefficients[2,4]))
    }
  }
  # ggplot(df.stats, aes(t, slope))+geom_point()+theme_article()

  df.stats <- df.stats[df.stats$slope>0,]
  if(is.null(df.stats)){
    warning("... could not fit any linear model")

  } else if(dim(df.stats)[1]==0){
    warning("... growing-window slopes are all negative")
  } else {
    ind_bests <- which(df.stats$r2>=0.7 & df.stats$pval<0.05)

    if(length(ind_bests)==0){
      ind_best <- which.max(df.stats$r2)
      warning(paste0("... linear chunk fitting did not satisfy the criteria, but a maximum for r2 was found: r2 = ", round(max(df.stats$r2)*100)/100))
    } else {
      t_bests <- df.stats$t[ind_bests]
      ind_not_1 <- which(diff(t_bests)>10)
      if(length(ind_not_1)==0){
        ind_best <- ind_bests[which.max(df.stats$r2[ind_bests])]
      } else {
        if(min(ind_not_1)>1){
          ind_best <- ind_bests[min(ind_not_1)-1]
        } else {
          ind_best <- which.max(df.stats$r2[-ind_not_1])
        }

      }
    }
    t_selected <- df.stats$t[ind_best]

    # p_r2 <- ggplot(df.stats, aes(t, r2))+geom_point()+theme_article()+xlab("time [secs]")+
    #   geom_point(data = df.stats[ind_bests,], aes(t, r2), color = "blue")+
    #   geom_vline(xintercept = t_selected, colour = "red")+xlim(c(0,max(df.stats$t)))+ggtitle("LM.r2 over a growing time window")
    # p_slope <- ggplot(df.stats, aes(t, slope))+geom_point()+theme_article()+xlab("time [secs]")+
    #   geom_point(data = df.stats[ind_bests,], aes(t, slope), color = "blue")+
    #   geom_vline(xintercept = t_selected, colour = "red")+xlim(c(0,max(df.stats$t)))+ggtitle("LM.slope over a growing time window")
    # p_pval <- ggplot(df.stats, aes(t, pval))+geom_point()+theme_article()+xlab("time [secs]")+
    #   geom_vline(xintercept = t_selected, colour = "red")+xlim(c(0,max(df.stats$t)))+ggtitle("p-val over a growing time window")+scale_y_log10()
    # p_conc <- ggplot(mydf, aes(time, conc))+geom_point()+theme_article()+
    #   geom_point(data = mydf[mydf$time >= df.stats$t[1] & mydf$time<=t_selected,], aes(time, conc), colour = "red")+
    #   ggtitle(unique(dataframe$UniqueID))+geom_vline(xintercept = t_selected, colour = "red")+xlim(c(0,max(df.stats$t)))
    # ggarrange(p_conc, p_r2, p_slope,  ncol = 1, align = "v")


    mystart = df.n$t[1]-length.min

    first_linear_chunk <- data.frame(UniqueID = unique(dataframe$UniqueID),
                                     start = mydf$time[which.min(abs(mydf$time - mystart))+1],
                                     end = t_selected)
  }

  first_linear_chunk$obs.length <- first_linear_chunk$end - first_linear_chunk$start
  first_linear_chunk$start.time <- mydf$POSIX.time[which.min(abs(mydf$time-first_linear_chunk$start))]

  return(first_linear_chunk)
}




find_first_linear_chunk.loop <- function(x, list_of_dataframes, gastype, length.min) {

  # Function to apply in the loop. Adapt parameters to your needs.
  chunk.out <- find_first_linear_chunk(dataframe = list_of_dataframes[[x]], gastype, length.min)

  return(chunk.out)
}
