# human second order
# make averages of pre : CS1 ("CS1-temp"), CS2 ("CS2 alone")
# first order: CS1US ("CS1-US"), CS1 ("CS1-temp") 
# pre second: CS2 ("CS2 alone")
# second order: CS2CS1 ("CS2-CS1"), CS1US ("CS1-US"), CS2 ("CS2 alone")
# post: CS1 ("CS1-temp"), CS2 ("CS2 alone")

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(lattice)
library(reshape2)

# stimTime function to split out
stimTime <- function(x){
  as.numeric(strsplit(x, "\t")[[1]][1])
}

# peak finding function
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

# standard errors
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
} 

# setwd
setwd(choose.dir())

# list subject folders
folders <- list.dirs(recursive = F)

# create empty dataframe
waveform_data <- as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(waveform_data) <- c("time", "value", "trial.type", "session", "subject")

# create list of markers for each file type
file_types <- list(c("CS1-temp", "CS2 alone"), c("CS1-US", "CS1-temp"), c("CS2 alone"), c("CS2-CS1", "CS1-US", "CS2 alone"), 
                   c("CS1-temp", "CS2 alone"))
names(file_types) <- c("preconditioning_unpaired", "first_order", "cs2_pre", "second_order", 
                       "postconditioning_unpaired")


# iterate through and compile a dataframe of values to create averages from
for(f in folders){
  # setwd
  setwd(f)
  
  # print the number its working on
  cat("working on:", which(folders==f), "out of", length(folders), "\n")
  cat(round(which(folders==f)/length(folders))*100, "% complete\n")
  
  # list of the text files
  file_list <- list.files()
  
  for(file in file_list){
    text <- readLines(file)
    
    start_point <- grep("NMR", text)[2] + 4
    # positions of "CHANNEL", waveform goes to next one - 2
    chan_pos <- grep("CHANNEL", text)
    # which one is the right index to stop at
    next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
    # read in waveform
    waveform <- as.numeric(text[start_point : next_chan_pos])
    # create time
    wv_df <- data.frame(time=seq(0, length.out = length(waveform), by=0.001), waveform=waveform)
    
    # choose revelant markers 
    file_session <- paste(strsplit(strsplit(file, split = ".", fixed = T)[[1]][1], "_")[[1]][1:2], collapse = "_")
    file_markers <- file_types[[file_session]]
    
    # iterate through the marker types (trial types)
    for(marks in file_markers){
      start_point <- grep(marks, text)[2] + 2
      next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
      next_chan_pos <- ifelse(is.infinite(next_chan_pos), length(text), next_chan_pos)
      markers <- sapply(text[start_point:next_chan_pos], stimTime, USE.NAMES = F)
      markers <- markers[!is.na(markers)]
      pre_time <- ifelse(grepl("CS2", marks), -0.2, -0.75)
      post_time <- ifelse(grepl("CS2", marks), 1.55, 1)
      
      # iterate through each marker and extract waveform with baseline as waveform minus average
      for(m in markers){
        wave <- wv_df$waveform[wv_df$time >= m + pre_time & wv_df$time <= m + post_time]
        wave <- wave - mean(wv_df$waveform[wv_df$time >= m - 0.2 & wv_df$time < m])
        wave_len <- length(wave)
        subject <- as.character(f)
        tmp_df <- data.frame(time = seq(-0.75, length.out = wave_len, by=0.001), waveform = wave, 
                             trial.type = rep(marks, wave_len), session = rep(file_session, wave_len),
                             subject = rep(subject, wave_len)) 
        waveform_data <- rbind(waveform_data, tmp_df)
      }
    }
  }
  setwd("..")
}

# save the waveform data as a csv
write.csv(waveform_data, file = "Human eyeblink waveform.csv", row.names = F)


# now process into averages for each trial type in each session per individual
average_df <- waveform_data %>%
  group_by(subject, session, trial.type, time) %>%
  summarise(value = mean(waveform),
            se = se(waveform))

# want to determine several things:
# based on the average is there a substantial movement before US time
# onset of that movement
# generate average plot of that
cr_response_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(cr_response_df) <- c("subject", "pre.CS1", "pre.CS2", "first.CS1US", "first.CS1temp", "cs2pre", "second.CS1US",
                           "second.CS2CS1", "second.CS2test",  "post.CS1", "post.CS2")

onset_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(onset_df) <- c("subject", "pre.CS1", "pre.CS2", "first.CS1US", "first.CS1temp", "cs2pre", "second.CS1US",
                           "second.CS2CS1", "second.CS2test",  "post.CS1", "post.CS2")


# peak analysis
peak_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(peak_df) <- colnames(onset_df)

# earliest peak latency
min_peak_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(min_peak_df) <- colnames(peak_df)

# peak amplitude
peak_amp_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(peak_amp_df) <- colnames(peak_df)

# earliest peak amplitude
min_peak_amp_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(min_peak_amp_df) <- colnames(peak_df)

# unique subjects
uni_subs <- as.character(unique(average_df$subject))

# iterate through
for(s in uni_subs){
  setwd(s)
  sub_df <- filter(average_df, subject==s, time < 0.8)
  
  # fill with onset times of average
  tmp_onset_df <- data.frame(subject=s, pre.CS1=NA, pre.CS2=NA, first.CS1US=NA, first.CS1temp=NA, cs2pre=NA, second.CS1US=NA,
                                second.CS2CS1=NA, second.CS2test=NA, post.CS1=NA, post.CS2=NA)
  
  tmp_cr_df <- tmp_onset_df
  
  tmp_peak_df <- tmp_onset_df
  
  tmp_min_peak <- tmp_onset_df
  
  tmp_peak_amp <- tmp_onset_df
  
  tmp_min_peak_amp <- tmp_onset_df
  
  # iterate through the relevent sessions
  uni_sessions <- as.character(unique(sub_df$session))
  
  for(sess in uni_sessions){
    dir.create(sess)
    setwd(sess)
    sess_df <- filter(sub_df, session==sess)
    
    # min and max values
    max_value <- max(sess_df$value + sess_df$se, na.rm = T) + 0.1
    min_value <- min(sess_df$value - sess_df$se, na.rm = T) - 0.1
    
    # iterate through the trial types
    uni_trials <- as.character(unique(sess_df$trial.type))
    
    for(trial in uni_trials){
      trial_df <- filter(sess_df, trial.type==trial)
      
      # assign the colour for each trial type: CS1 - purple ("#C77CFF"), CS1US - red ("#F8766D"),
      # CS2CS1 - blue ("#619CFF"), CS2 - green ("#00BA38")
      col <- ifelse(grepl("alone", trial), "#00BA38", ifelse(grepl("CS2", trial), "#619CFF", ifelse(grepl("US", trial), "#F8766D", "#C77CFF")))
      trial_type <- ifelse(grepl("alone", trial), "CS2 alone", ifelse(grepl("CS2", trial), "CS2-CS1", ifelse(grepl("US", trial), "CS1-US", "CS1 alone")))
      
      # determine onset if there is one
      baseline_begin <- ifelse(grepl("CS2", trial_type), -0.75, -0.2)
      baseline_period <- trial_df$value[trial_df$time > baseline_begin & trial_df$time <= baseline_begin + 0.2]
      
      # blink period, if CS1 -> till US time, if CS2 first check up to CS1 time
      blink_period <- trial_df$value[trial_df$time > baseline_begin + 0.2 & trial_df$time <= ifelse(grepl("CS2", trial_type), 0, 0.45)]
      
      # blink amplitude and latency
      blink_amplitude <- trial_df$value[trial_df$time > baseline_begin + 0.2 & trial_df$time <= 0.45]
      
      # absolute peak
      peak_value <- max(blink_amplitude)
      peak_latency <- min(which(blink_amplitude == peak_value)) * 0.001 
      peak_latency <- ifelse(grepl("CS2", trial_type), peak_latency - 0.55, peak_latency)
      
      # min peak
      min_peak <- min(find_peaks(blink_amplitude, m=150))
      min_peak_value <- blink_amplitude[min_peak]
      min_peak_latency <- min_peak * 0.001
      min_peak_latency <- ifelse(grepl("CS2", trial_type), min_peak_latency - 0.55, min_peak_latency)
      
      # threshold is there a blink?
      if(any(blink_period > sd(baseline_period)*15)){
        
        # smooth spline
        blink_smo <- smooth.spline(blink_period, cv=T)$y
        
        # max derivative across the window
        max_deriv <- numeric(0)
        for(i in 0:9){
          window_index <- c(1+i, seq(10+i, length(blink_smo)+1 - (10-(i-1)), 10))
          second_deriv <- diff(blink_smo[window_index], differences = 2)
          max_deriv <- c(max_deriv, which.max(second_deriv))
        }
        # where is the onset
        blink_onset <- (baseline_begin + 0.2) + 0.01*as.numeric(names(which.max(table(max_deriv)))) + 0.005
        blink_value <- mean(trial_df$value[trial_df$time > blink_onset - 0.002 & trial_df$time < blink_onset + 0.002])
        blink_amp <- peak_value - blink_value
        min_blink_amp <- min_peak_value - blink_value
        
        # there is a cr
        cr <- T
        
      } else if(grepl("CS2", trial_type) & !any(blink_period > sd(baseline_period)*15)){
        # check the cs1 period for a response
        blink_period <- trial_df$value[trial_df$time > 0 & trial_df$time <=  0.45]
        
        if(any(blink_period > sd(baseline_period)*15)){
          # smooth spline
          blink_smo <- smooth.spline(blink_period, cv=T)$y
          
          # max derivative across the window
          max_deriv <- numeric(0)
          for(i in 0:9){
            window_index <- c(1+i, seq(10+i, length(blink_smo)+1 - (10-(i-1)), 10))
            second_deriv <- diff(blink_smo[window_index], differences = 2)
            max_deriv <- c(max_deriv, which.max(second_deriv))
          }
          # where is the onset
          blink_onset <- 0.01*as.numeric(names(which.max(table(max_deriv)))) + 0.005
          blink_value <- mean(trial_df$value[trial_df$time > blink_onset - 0.002 & trial_df$time < blink_onset + 0.002])
          blink_amp <- peak_value - blink_value
          min_blink_amp <- min_peak_value - blink_value
          
          # there is a cr
          cr <- T
        } else {
          blink_onset <- NA
          cr <- F
          blink_amp <- NA
          min_blink_amp <- NA
          peak_latency <- NA
          min_peak_latency <- NA
        }
        
      } else {
        blink_onset <- NA
        cr <- F
        blink_amp <- NA
        min_blink_amp <- NA
        peak_latency <- NA
        min_peak_latency <- NA
      }
      
    
      # put the onset time in to a database
      if(sess=="preconditioning_unpaired"){
        if(trial == "CS2 alone"){
          tmp_onset_df$pre.CS2 <- blink_onset
          tmp_cr_df$pre.CS2 <- cr
          tmp_peak_df$pre.CS2 <- peak_latency
          tmp_min_peak$pre.CS2 <- min_peak_latency
          tmp_peak_amp$pre.CS2 <- blink_amp
          tmp_min_peak_amp$pre.CS2 <- min_blink_amp
        } else {
          tmp_onset_df$pre.CS1 <- blink_onset
          tmp_cr_df$pre.CS1 <- cr
          tmp_peak_df$pre.CS1 <- peak_latency
          tmp_min_peak$pre.CS1 <- min_peak_latency
          tmp_peak_amp$pre.CS1 <- blink_amp
          tmp_min_peak_amp$pre.CS1 <- min_blink_amp
        }
      }
      
      if(sess=="first_order"){
        if(trial=="CS1-US"){
          tmp_onset_df$first.CS1US <- blink_onset
          tmp_cr_df$first.CS1US <- cr
          tmp_peak_df$first.CS1US <- peak_latency
          tmp_min_peak$first.CS1US <- min_peak_latency
          tmp_peak_amp$first.CS1US <- blink_amp
          tmp_min_peak_amp$first.CS1US <- min_blink_amp
        } else {
          tmp_onset_df$first.CS1temp <- blink_onset
          tmp_cr_df$first.CS1temp <- cr
          tmp_peak_df$first.CS1temp <- peak_latency
          tmp_min_peak$first.CS1temp <- min_peak_latency
          tmp_peak_amp$first.CS1temp <- blink_amp
          tmp_min_peak_amp$first.CS1temp <- min_blink_amp
        }
      }
      
      if(sess=="cs2_pre"){
        tmp_onset_df$cs2pre <- blink_onset
        tmp_cr_df$cs2pre <- cr
        tmp_peak_df$cs2pre <- peak_latency
        tmp_min_peak$cs2pre <- min_peak_latency
        tmp_peak_amp$cs2pre <- blink_amp
        tmp_min_peak_amp$cs2pre <- min_blink_amp
      }
      
      if(sess=="second_order"){
        if(trial=="CS2-CS1"){
          tmp_onset_df$second.CS2CS1 <- blink_onset
          tmp_cr_df$second.CS2CS1 <- cr
          tmp_peak_df$second.CS2CS1 <- peak_latency
          tmp_min_peak$second.CS2CS1 <- min_peak_latency
          tmp_peak_amp$second.CS2CS1 <- blink_amp
          tmp_min_peak_amp$second.CS2CS1 <- min_blink_amp
        } else if(trial=="CS2 alone"){
          tmp_onset_df$second.CS2test <- blink_onset
          tmp_cr_df$second.CS2test <- cr
          tmp_peak_df$second.CS2test <- peak_latency
          tmp_min_peak$second.CS2test <- min_peak_latency
          tmp_peak_amp$second.CS2test <- blink_amp
          tmp_min_peak_amp$second.CS2test <- min_blink_amp
        } else {
          tmp_onset_df$second.CS1US <- blink_onset
          tmp_cr_df$second.CS1US <- cr
          tmp_peak_df$second.CS1US <- peak_latency
          tmp_min_peak$second.CS1US <- min_peak_latency
          tmp_peak_amp$second.CS1US <- blink_amp
          tmp_min_peak_amp$second.CS1US <- min_blink_amp
        }
      }
      
      if(sess=="postconditioning_unpaired"){
        if(trial=="CS1-temp"){
          tmp_onset_df$post.CS1 <- blink_onset
          tmp_cr_df$post.CS1 <- cr
          tmp_peak_df$post.CS1 <- peak_latency
          tmp_min_peak$post.CS1 <- min_peak_latency
          tmp_peak_amp$post.CS1 <-  blink_amp
          tmp_min_peak_amp$post.CS1 <- min_blink_amp
        } else {
          tmp_onset_df$post.CS2 <- blink_onset
          tmp_cr_df$post.CS2 <- cr
          tmp_peak_df$post.CS2 <- peak_latency
          tmp_min_peak$post.CS2 <- min_peak_latency
          tmp_peak_amp$post.CS2 <- blink_amp
          tmp_min_peak_amp$post.CS2 <- min_blink_amp
        }
      }
      
      # still need a value for blink but off the plot margins if it is NA
      if(is.na(blink_onset)){
        blink_value <- -100
        blink_onset <- -100
      }
      
      
      # plot of the average waveform
      waveform_plot <- ggplot() + geom_line(data=trial_df, aes(x=time, y=value), colour=col, size=1.5) + 
        geom_ribbon(data=trial_df, aes(x=time, ymin = value - se, ymax = value + se), alpha=0.3) +
        theme_classic() + scale_y_continuous(limits = c(min_value, max_value), name = "arbitary value (V)") +
        scale_x_continuous(limits=c(min(trial_df$time), 0.75), name = "time relative to CS1 stimulus onset (ms)",
                           labels = function(x) x*1000, breaks = c(-0.55, 0, 0.45), expand = c(0,0)) + 
        geom_vline(xintercept = c(-0.55, 0, 0.45), alpha=0.4, linetype=2) + ggtitle(trial_type) +
        geom_point(aes(x=blink_onset, y=blink_value), size=6, alpha=0.6, colour="#ff8000")
      
      if(grepl("CS2", trial_type)){
        waveform_plot <- waveform_plot + annotate("text", x = -0.46, y = max_value, label="CS2 onset")
      }
      if(grepl("US", trial_type)){
        waveform_plot <- waveform_plot + annotate("text", x = 0.53, y = max_value, label = "US onset") 
      }
      if(grepl("CS1", trial_type)){
        waveform_plot <- waveform_plot + annotate("text", x = 0.09, y = max_value, label = "CS1 onset")
      }
      
      ggsave(waveform_plot, filename = paste0(sess, "_", trial, ".pdf"), width = 9, height = 6)
      
    }
    setwd("..")
  }
  setwd("..")
  
  # row bind the response dataframe with the temporary dataframe
  onset_df <- rbind(onset_df, tmp_onset_df)
  cr_response_df <- rbind(cr_response_df, tmp_cr_df)
  peak_df <- rbind(peak_df, tmp_peak_df)
  min_peak_df <- rbind(min_peak_df, tmp_min_peak)
  peak_amp_df <- rbind(peak_amp_df, tmp_peak_amp)
  min_peak_amp_df <- rbind(min_peak_amp_df, tmp_min_peak_amp)
}
write.csv(onset_df, file = "onset times for trial types in different sessions.csv", row.names = F)
write.csv(cr_response_df, file = "cr responses.csv", row.names = F)


# want to create onset difference for those that have onsets for both
onsets <- onset_df %>%
  select(subject, second.CS1US, second.CS2CS1) %>%
  na.omit(onsets) 

# melt the data into long  
onsets <- melt(onsets, id.vars="subject")

onset_plot <- ggplot(onsets, aes(x=variable, y=value, colour=subject))+geom_point(size=2) + geom_line(aes(group=subject), size=1) +
  theme_classic() + scale_x_discrete(labels=c("second.CS1US"="CS1-US", "second.CS2CS1"="CS2-CS1"), name="Conditioned Stimuli") +
  scale_y_continuous(name = "Onset Latency Relative to CS1 (ms)", labels=function(x)x*1000) + 
  guides(colour="none")

ggsave(onset_plot, filename = "dot plot of onset latency.pdf", width=6, height=4)

# paired t test
paired_t <- t.test(x=onsets$value[onsets$variable=="second.CS1US"], y=onsets$value[onsets$variable=="second.CS2CS1"],
                   paired=T)

sink("paired t test of onsets.txt")
print(paired_t)
sink()

# write csv for peaks
write.csv(peak_df, file = "peak latencies.csv", row.names = F)
write.csv(min_peak_df, file = "first peak latencies.csv", row.names = F)
write.csv(peak_amp_df, file = "peak amplitudes.csv", row.names = F)
write.csv(min_peak_amp_df, file = "first peak amplitudes.csv", row.names=F)
