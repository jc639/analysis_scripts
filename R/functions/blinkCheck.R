# function to determine if blink has occurred in unrestrained animals
# takes in markers waveform, trial type and sample rate of sample acquistion

blinkCheck <- function(marker, waveform, trial.type, sample.rate){
  premarker <- marker - 0.4
  pre_premarker <- premarker - 0.1
  
  if (trial.type == "CS2CS1" | trial.type == "CS2"){
    response_period <- marker + 0.8
  } else {
    response_period <- marker + 0.35
  }
  
  min_wv <- min(waveform$Eye.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  max_wv <- max(waveform$Eye.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  mean_wv <- mean(waveform$Eye.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  sd_wv <- sd(waveform$Eye.position[waveform$Time > pre_premarker & waveform$Time <= premarker], na.rm = T)
  
  if((max_wv - min_wv) > 0.1){
    return(c("not flat before", NA, NA, NA))
  }
  
  baseline_wv <- waveform$Eye.position[waveform$Time > premarker & waveform$Time <= marker]
  
  if(any(baseline_wv >  (mean_wv + 0.1))){
    onset <- min(which(baseline_wv > (mean_wv + 5*sd_wv))) * sample.rate
    return(c("not flat baseline", onset, NA, NA))
  } else {
    mean_wv <- mean(baseline_wv)
    sd_wv <- sd(baseline_wv)
  }
  
  response_wv <- waveform$Eye.position[waveform$Time > marker & waveform$Time <= response_period]
  if(any(response_wv > (mean_wv + 0.1))){
    onset <- min(which(response_wv > (mean_wv + 5*sd_wv))) * sample.rate
    peak_index <- which.max(response_wv)
    peak <- peak_index * sample.rate
    amplitude <- response_wv[peak_index] - mean_wv
    alpha <- onset < 0.035
    if(alpha){
      if(trial.type == "CS2CS1" | trial.type == "CS2"){
        onset <- onset - 0.45
        peak <- peak - 0.45
      }
      return(c("alpha response", onset, peak, amplitude))
    } else{
      if(trial.type == "CS2CS1" | trial.type == "CS2"){
        onset <- onset - 0.45
        peak <- peak - 0.45
      }
      return(c("CR", onset, peak, amplitude))
    }
  } else {
    return(c("NO CR", NA, NA, NA))
  }
}
  

