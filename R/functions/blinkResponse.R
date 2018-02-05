# function to determine what the response is to the stimulus
blinkResponse <- function(marker, waveform, trial.type, sample.rate){
  
  premarker <- marker - 0.5
  pre_premarker <- premarker - 0.2
  
  if (trial.type == "CS2CS1" | trial.type == "CS2"){
    response_period <- marker + 0.7
  } else {
    response_period <- marker + 0.35
  }
  
  min_wv <- min(waveform$NMR.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  max_wv <- max(waveform$NMR.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  mean_wv <- mean(waveform$NMR.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  
  if((max_wv - min_wv) > 0.1){
    return(c("not flat before", NA, NA, NA))
  }
  
  baseline_wv <- waveform$NMR.position[waveform$Time > premarker & waveform$Time <= marker]
  if(any(baseline_wv >  (mean_wv + 0.15))){
    onset <- min(which(baseline_wv > (mean_wv + 0.3))) * sample.rate
    return(c("not flat baseline", onset, NA, NA))
  } else {
    mean_wv <- mean(baseline_wv)
  }
  
  response_wv <- waveform$NMR.position[waveform$Time > marker & waveform$Time <= response_period]
  if(any(response_wv > (mean_wv + 0.3))){
    onset <- min(which(response_wv > (mean_wv + 0.3))) * sample.rate
    peak_index <- which.max(response_wv)
    peak <- peak_index * sample.rate
    amplitude <- response_wv[peak_index] - mean_wv
    alpha <- onset < 0.035
    if(alpha){
      if(trial.type == "CS2CS1" | trial.type == "CS2"){
        onset <- onset - 0.35
        peak <- peak - 0.35
      }
      return(c("alpha response", onset, peak, amplitude))
    } else{
      if(trial.type == "CS2CS1" | trial.type == "CS2"){
        onset <- onset - 0.35
        peak <- peak - 0.35
      }
      return(c("CR", onset, peak, amplitude))
    }
  } else {
    return(c("NO CR", NA, NA, NA))
  }
}


# function that reads in waveform data and markers and returns a df for response
blinkDataframes <- function(markers=NULL, waveform=NULL, trial.type=NULL, subject=NULL, date=NULL,
                            session.n = NULL, sample.rate = NULL){
  response_df <- as.data.frame(matrix(nrow = length(markers), ncol = 9))
  colnames(response_df) <- c("response", "trial.type", "subject", "trial.n", "session.n", "date", "onset", "peak", 
                    "amplitude")
  
  response_df$trial.type <- trial.type
  response_df$subject <- subject
  response_df$date <- date
  response_df$session.n <- session.n
  
  for (i in 1:length(markers)){
    response_df$trial.n[i] <- i
    response_df[i, ][c("response", "onset", "peak", "amplitude")] <- blinkResponse(marker = markers[i], waveform = waveform,
                                                                          trial.type = trial.type, sample.rate = sample.rate)
  }
  
  response_df[c("onset", "peak", "amplitude")] <- apply(response_df[c("onset", "peak", "amplitude")], 2, as.numeric)
  
  for (i in 1:length(markers)){
    if(trial.type == "CS2CS1" | trial.type == "CS2"){
      start_mark <- markers[i] - 0.15
      end_mark <- markers[i] + 0.8
    } else {
      start_mark <- markers[i] - 0.5
      end_mark <- markers[i] + 0.45
    }
    
    trial_wv <- waveform$NMR.position[waveform$Time > start_mark & waveform$Time <= end_mark]
    baseline_avg <- mean(waveform$NMR.position[waveform$Time > (markers[i]-0.5) & waveform$Time <= markers[i]], na.rm = T)
    trial_wv <- trial_wv - baseline_avg
    
    temp_df <- as.data.frame(matrix(nrow=length(trial_wv), ncol=8))
    colnames(temp_df) <- c("time", "value", "trial.type", "response", "date", "session.n", "trial.n", "subject")
    
    temp_df$value <- trial_wv
    temp_df$time <- seq(-0.5, length.out = length(trial_wv), by = sample.rate)
    temp_df$time <- temp_df$time * 1000
    temp_df$trial.type <- trial.type
    temp_df$response <- response_df$response[i]
    temp_df$trial.n <- i
    temp_df$session.n <- session.n
    temp_df$date <- date
    temp_df$subject <- subject
    
    if(i == 1){
      trial_wv_df <- temp_df
    } else {
      trial_wv_df <- rbind(trial_wv_df, temp_df)
    }
  }
  
  dataframes <- list(response_df, trial_wv_df) 
  return(dataframes)
}