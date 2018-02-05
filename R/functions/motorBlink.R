# function to assess motor eyeblink, incidence, onset latency, peak latency
motorBlink <- function(waveform = NULL, button = NULL, us = NULL){
  if(any(is.null(c(waveform, button, us)))){
    stop("Not passed correct data to motorBlink function:\none or more parameters are empty")
  }
  
  # determine average blink size
  us_region <- us + 0.5
  avg_blink <- numeric(0)
  for(i in 1:length(us)){
    ur_waveform <- waveform$waveform[waveform$time > us[i] & waveform$time < us_region[i]]
    avg_blink <- c(avg_blink, max(ur_waveform) - min(ur_waveform))
  }
  avg_blink <- mean(avg_blink)
  
  # dataframe for results
  button$CR <- NA
  button$onset <- NA
  button$peak <- NA
  
  # dataframe for trial waveforms
  trial_wv_df <- as.data.frame(matrix(nrow=0, ncol=8))
  colnames(trial_wv_df) <- c("time", "value", "response", "button", "id", "group", "homepad", "trial.n")
  
  # iterate through each trial
  for(i in 1 : nrow(button)){
    # extract waveform for analysis
    wv <- waveform$waveform[waveform$time >= button$time[i] - 0.5 & waveform$time < button$time[i] + 0.5]
    # determine if it is flat in first 100ms
    not_flat <- abs(max(wv[1 : 100]) - min(wv[1 : 100])) > (avg_blink / 10)
    if(not_flat){
      # not flat
      button$CR[i] <- "not_flat"
    } else{
      # work out average of first 200ms and SD
      avg_wv <- mean(wv[1 : 200])
      sd_wv <- sd(wv[1 : 200])
      # if any wv goes lower than 10th of blink size
      if(any(wv[201 : length(wv)] < avg_wv - (avg_blink * 0.05))){
        # if any is greater than 3 SD of baseline, CR
        if(any(wv < avg_wv - (sd_wv * 10))){
          index <- min(which(wv[201 : length(wv)] < avg_wv - (sd_wv * 10)))
          onset <- (index * 0.001) - 0.3
          button$CR[i] <- ifelse(onset <= 0, "CR", "No_CR")
          button$onset[i] <- onset
          button$peak[i] <- min(which(wv == min(wv))) * 0.001 - 0.5
        } else {
          button$CR[i] <- "No_CR"
        }
      } else {
        button$CR[i] <- "No_CR"
      }
    }
    tmp_wv <- data.frame(time = seq(-0.5, length.out = length(wv), by = 0.001), value = wv - mean(wv[1:200]),
                         response = rep(button$CR[i], length(wv)), button = rep(button$button[i], length(wv)), id = rep(button$id[i], length(wv)), 
                         group = rep(button$group[i], length(wv)), homepad = rep(button$RT[i], length(wv)),
                         trial.n = rep(button$press.N[i], length(wv)))
    trial_wv_df <- rbind(trial_wv_df, tmp_wv)
  }
  
  dfs <- c(list(button), list(trial_wv_df))
  
  return(dfs)
}