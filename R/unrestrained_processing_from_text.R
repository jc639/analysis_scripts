# eyeblink unrestrained
# how to determine blink?
# potentially voltage change on UR
# and SD thresholds from baselines...

# function to check for CR/NO CR, not flat before
source("functions/blinkCheck.R")

# function to remove any markers that appear in gaps
inGaps <- function(markers, start, end){
  in_gaps <- numeric(0)
  for(i in 1:length(start)){
    for(j in 1:length(markers)){
      in_gaps <- append(in_gaps, markers[j] > start[i] & markers[j] < end[i])
    }
    markers <- markers[!in_gaps]
    in_gaps <- numeric(0)
  }
  return(markers)
}

# setwd to folder containing second order files
setwd(choose.dir())
home_dir <- getwd()

# list animals folder
folders <- list.dirs(recursive = F)

# response dataframe
response_df <- as.data.frame(matrix(ncol=9, nrow=0))
colnames(response_df) <- c("trial.type", "trial.n", "response","onset", "peak", "amplitude", "subject", "date", "first.order")

# trial waveform data frame
trial_wave <- as.data.frame(matrix(ncol=6, nrow=0))
colnames(trial_wave) <- c("time", "value", "subject", "trial.type", "response", "date")

# iterate into each folder
for (f in folders){
  # f is subject folder, all text files in here
  setwd(f)
  cat(paste0("working on ", f, "\n"))
  
  # list all text files of sessions
  text_files <- list.files()
  
  # iterate through the text files
  for(t in text_files){
    # date is second element in text title, and strip away .txt
    date <- strsplit(strsplit(t, "_")[[1]][2], "\\.")[[1]][1]
    # subject is first element
    subject <- strsplit(t, "_")[[1]][1]
    cat(paste0("currently on ", date, "\n"))
    
    # IR data start point
    start_point <- which(grepl("START",readLines(t))) + 1
    sample_rate <- as.numeric(strsplit(readLines(t)[start_point-1], split="\t")[[1]][3])
    
    # IR starts on line after "Start" and runs to 3 lines before second "Evt+"
    IR <- readLines(t)[(start_point) : (which(grepl("Evt+", readLines(t)))[2] -3)]
    IR <- IR[IR != ""]
    
    # small number of sessions have gaps where pause button was used
    gaps <- IR[grepl("GAP", IR)]
    
    #if no gaps then length of grepl should be zero
    if(length(gaps) != 0){
      ref_table <- as.data.frame(matrix(ncol = 2,nrow = length(gaps)))
      colnames(ref_table) <- c("Index", "Time")
      
      # reference table for gap starts
      for(i in 1:length(gaps)){
        ref_table$Index[i] <- which(IR == gaps[i])
        ref_table$Time[i] <- as.numeric(strsplit(gaps[i], split="\t")[[1]][2])
        IR <- IR[IR != gaps[i]]
      }
      
      # convert to numeric
      IR <- as.numeric(IR)
      
      # waveform data frame, eye position at time
      waveform_data <- as.data.frame(matrix(nrow=length(IR),ncol=2))
      colnames(waveform_data) <- c("Eye.position", "Time")
      
      # IR is is eyelid position
      waveform_data$Eye.position <- IR
      
      # start time is 0, index 1
      ref_table <- rbind(c(1, 0), ref_table)
      ref_table <- rbind(ref_table, c(nrow(waveform_data), NA))
      
      # gap start and end empty list
      gap_start <- numeric(0)
      gap_end <- numeric(0)
      
	  # iterate through
      for(i in 1:(nrow(ref_table) -1)){
        if (i == (nrow(ref_table) -1)){
          waveform_data$Time[ref_table$Index[i] : ref_table$Index[i + 1]] <- seq(ref_table$Time[i], 
                                                                               length.out = nrow(waveform_data[ref_table$Index[i] : ref_table$Index[i + 1],]),
                                                                               by = sample_rate)
        } else {
          waveform_data$Time[ref_table$Index[i]:(ref_table$Index[i + 1] -1)] <- seq(ref_table$Time[i], 
                                                                                 length.out = ref_table$Index[i + 1]-ref_table$Index[i],
                                                                                 by=sample_rate)
          gap_start <- append(gap_start, waveform_data$Time[ref_table$Index[i + 1] -1])
          gap_end <- append(gap_end, ref_table$Time[i + 1])
        }
      }
    } else {
      # convert to numeric
      IR <- as.numeric(IR)
      
      # waveform data
      waveform_data <- as.data.frame(matrix(nrow = length(IR), ncol = 2))
      colnames(waveform_data) <- c("Eye.position", "Time")
      
      waveform_data$Eye.position <- IR
      waveform_data$Time<-seq(0, length.out = nrow(waveform_data),  by = sample_rate)		
    }
    
    # need to string split lines with markers on and convert to numeric
    # firstElement does this
    firstElement <- function(x){
      element <- numeric(0)
      for (i in 1:length(x)){
        element <- c(element, strsplit(x[i], split = "\t")[[1]][1])
      }
      return(as.numeric(element))
    }
    
    CS1US_markers <- readLines(t)[(which(grepl("CS1-US", readLines(t)))[2] + 2):(which(grepl("CS2-CS1", readLines(t)))[2] - 5)]
    CS1US_markers <- na.omit(firstElement(CS1US_markers))
    if(length(gaps) != 0){
      CS1US_markers <- inGaps(markers = CS1US_markers, start = gap_start, end = gap_end)
    }
    
    CS2CS1_markers <- readLines(t)[(which(grepl("CS2-CS1", readLines(t)))[2] + 2):(which(grepl("CS1-temp", readLines(t)))[2] - 5)]
    CS2CS1_markers <- na.omit(firstElement(CS2CS1_markers))
    if(length(gaps) != 0){
      CS2CS1_markers <- inGaps(markers = CS2CS1_markers, start = gap_start, end = gap_end)
    }
    
    CS1test_markers <- readLines(t)[(which(grepl("CS1test", readLines(t)))[2] + 2):(length(readLines(t)))]
    CS1test_markers <- na.omit(firstElement(CS1test_markers))
    if(length(gaps) != 0){
      CS1test_markers <- inGaps(markers = CS1test_markers, start = gap_start, end = gap_end)
    }
    
    CS2test_markers <- readLines(t)[(which(grepl("CS2 alone", readLines(t)))[2] + 2):(which(grepl("CS1test", readLines(t)))[2] - 5)]
    CS2test_markers <- na.omit(firstElement(CS2test_markers))
    if(length(gaps) != 0){
      CS2test_markers <- inGaps(markers = CS2test_markers, start = gap_start, end = gap_end)
    }
    
    # need to iterate through each marker
    # determine flat baseline
    # and deviation from baseline to blink
    
    for(i in 1:length(CS1US_markers)){
      blinks <- blinkCheck(marker = CS1US_markers[i],
                           waveform = waveform_data,
                           trial.type = "CS1US", 
                           sample.rate = sample_rate)
      tmp_df <- data.frame(trial.type = "CS1US", trial.n = i, response =  blinks[1],
                           onset = as.numeric(blinks[2]), peak = as.numeric(blinks[3]), amplitude = as.numeric(blinks[4]),
                           subject = subject, date = date, first.order = date%in%c("0517", "0518", "0519", "0520"))
      response_df <- rbind(response_df, tmp_df)
      
      if(blinks[1] %in% c("CR", "NO CR")){
        wv <- waveform_data$Eye.position[waveform_data$Time > CS1US_markers[i] - 0.4 & waveform_data$Time <= CS1US_markers[i] + 0.45]
        mean_wv <- mean(waveform_data$Eye.position[waveform_data$Time > CS1US_markers[i] - 0.4 & waveform_data$Time <= CS1US_markers[i]], na.rm = T)
        wv <- wv - mean_wv
        wv_tmp_df <- as.data.frame(matrix(ncol=6, nrow=length(wv)))
        colnames(wv_tmp_df) <- c("time", "value", "subject", "trial.type", "response", "date")
        wv_tmp_df$time <- seq(-0.4, length.out = length(wv), by = sample_rate)
        wv_tmp_df$value <- wv
        wv_tmp_df$subject <- subject
        wv_tmp_df$trial.type <- "CS1US"
        wv_tmp_df$response <- blinks[1]
        wv_tmp_df$date <- date
        
        trial_wave <- rbind(trial_wave, wv_tmp_df)
      }
    }
    
    if(length(CS2CS1_markers) != 0){
      for(i in 1:length(CS2CS1_markers)){
        blinks <- blinkCheck(marker=CS2CS1_markers[i],
                             waveform = waveform_data,
                             trial.type = "CS2CS1",
                             sample.rate = sample_rate)
        tmp_df <- data.frame(trial.type = "CS2CS1", trial.n = i, response =  blinks[1],
                             onset = as.numeric(blinks[2]), peak = as.numeric(blinks[3]), amplitude = as.numeric(blinks[4]),
                             subject = subject, date = date, first.order = date%in%c("0517", "0518", "0519", "0520"))
        response_df <- rbind(response_df, tmp_df)
        
        if(blinks[1] %in% c("CR", "NO CR")){
          wv <- waveform_data$Eye.position[waveform_data$Time > CS2CS1_markers[i] - 0.4 & waveform_data$Time <= CS2CS1_markers[i] + 0.9]
          mean_wv <- mean(waveform_data$Eye.position[waveform_data$Time > CS2CS1_markers[i] - 0.4 & waveform_data$Time <= CS2CS1_markers[i]], na.rm = T)
          wv <- wv - mean_wv
          wv_tmp_df <- as.data.frame(matrix(ncol=6, nrow=length(wv)))
          colnames(wv_tmp_df) <- c("time", "value", "subject", "trial.type", "response", "date")
          wv_tmp_df$time <- seq(-0.85, length.out = length(wv), by = sample_rate)
          wv_tmp_df$value <- wv
          wv_tmp_df$subject <- subject
          wv_tmp_df$trial.type <- "CS2CS1"
          wv_tmp_df$response <- blinks[1]
          wv_tmp_df$date <- date
          
          trial_wave <- rbind(trial_wave, wv_tmp_df)
        }
      }
    }
    
    for(i in 1:length(CS1test_markers)){
      blinks <- blinkCheck(marker=CS1test_markers[i],
                           waveform = waveform_data,
                           trial.type = "CS1",
                           sample.rate = sample_rate)
      tmp_df <- data.frame(trial.type = "CS1", trial.n = i, response =  blinks[1],
                           onset = as.numeric(blinks[2]), peak = as.numeric(blinks[3]), amplitude = as.numeric(blinks[4]),
                           subject = subject, date = date, first.order = date%in%c("0517", "0518", "0519", "0520"))
      response_df <- rbind(response_df, tmp_df)
      
      if(blinks[1] %in% c("CR", "NO CR")){
        wv <- waveform_data$Eye.position[waveform_data$Time > CS1test_markers[i] - 0.4 & waveform_data$Time <= CS1test_markers[i] + 0.45]
        mean_wv <- mean(waveform_data$Eye.position[waveform_data$Time > CS1test_markers[i] - 0.4 & waveform_data$Time <= CS1test_markers[i]], na.rm = T)
        wv <- wv - mean_wv
        wv_tmp_df <- as.data.frame(matrix(ncol=6, nrow=length(wv)))
        colnames(wv_tmp_df) <- c("time", "value", "subject", "trial.type", "response", "date")
        wv_tmp_df$time <- seq(-0.4, length.out = length(wv), by = sample_rate)
        wv_tmp_df$value <- wv
        wv_tmp_df$subject <- subject
        wv_tmp_df$trial.type <- "CS1"
        wv_tmp_df$response <- blinks[1]
        wv_tmp_df$date <- date
        
        trial_wave <- rbind(trial_wave, wv_tmp_df)
      }
    }
    
    if(length(CS2test_markers) != 0){
      for(i in 1:length(CS2test_markers)){
        blinks <- blinkCheck(marker = CS2test_markers[i],
                             waveform = waveform_data,
                             trial.type = "CS2",
                             sample.rate = sample_rate)
        tmp_df <- data.frame(trial.type = "CS2", trial.n = i, response =  blinks[1],
                             onset = as.numeric(blinks[2]), peak = as.numeric(blinks[3]), amplitude = as.numeric(blinks[4]),
                             subject = subject, date = date, first.order = date%in%c("0517", "0518", "0519", "0520"))
        response_df <- rbind(response_df, tmp_df)
        
        if(blinks[1] %in% c("CR", "NO CR")){
          wv <- waveform_data$Eye.position[waveform_data$Time > CS2test_markers[i] - 0.4 & waveform_data$Time <= CS2test_markers[i] + 0.9]
          mean_wv <- mean(waveform_data$Eye.position[waveform_data$Time > CS2test_markers[i] - 0.4 & waveform_data$Time <= CS2test_markers[i]], na.rm = T)
          wv <- wv - mean_wv
          wv_tmp_df <- as.data.frame(matrix(ncol=6, nrow=length(wv)))
          colnames(wv_tmp_df) <- c("time", "value", "subject", "trial.type", "response", "date")
          wv_tmp_df$time <- seq(-0.85, length.out = length(wv), by = sample_rate)
          wv_tmp_df$value <- wv
          wv_tmp_df$subject <- subject
          wv_tmp_df$trial.type <- "CS2"
          wv_tmp_df$response <- blinks[1]
          wv_tmp_df$date <- date
          
          trial_wave <- rbind(trial_wave, wv_tmp_df)
        }
      }
    }
  }
  setwd(home_dir)
}
# finally save as CSVs for further analysis
write.csv(response_df, "response dataframe.csv", row.names = F)
write.csv(trial_wave, "trial waveforms.csv", row.names = F)


