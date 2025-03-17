


find_bubbles <- function(time, conc, window.size){
  # standarizing conc
  conc_std <- (conc-min(conc))/max(conc)

  x = seq(min(time),max(time),1)
  conc_std <- approx(time, conc_std, xout = x, method = "linear", rule = 2)$y


  # defining start and end for moving window
  t_start <- x[1]+ceiling(window.size/2)
  t_stop <- max(time)-floor(window.size/2)

  # computing variance of conc_std with a moving window
  df.stats <- NULL
  for(t in seq(t_start, floor(t_stop),1)){
    ind_window <- which(x>=t-floor(window.size/2) & x<=t+floor(window.size/2))
    df.stats <- rbind(df.stats, data.frame(t = t,
                                           var = var(conc_std[ind_window])))
  }

  # defining threshold
  thresh <- max(1e-4, quantile(df.stats$var, 0.7))

  # plot(df.stats$t, df.stats$var)
  # points(df.stats$t[df.stats$var > thresh], df.stats$var[df.stats$var > thresh], col='red')

  # finding chunks with high variance
  vect <- df.stats$var > thresh

  #initializing
  chunks <- NULL
  if(sum(vect)>0){ # in that case, there are some possible bubbling events to check

    jumps <- NULL
    if(vect[1]==T){jumps <- c(1)}
    for(i in seq(2,length(vect))){
      if(vect[i] != vect[i-1]){
        jumps <- c(jumps, i)
      }
    }
    if(last(vect)==T){jumps <- c(jumps, i)}
    chunks <- data.frame(start = jumps[seq(1,length(jumps),2)],
                         end = jumps[seq(2,length(jumps),2)])

    # we discard chunks of less than 5 datapoints
    chunks <- chunks[which(chunks$end-chunks$start > 5),]
    if(dim(chunks)[1] == 0){
      chunks <- NULL
    }

    # plot(x, conc_std)
    # abline(v = c(chunks$start, chunks$end), col='red')
  }

  return(chunks)
}

