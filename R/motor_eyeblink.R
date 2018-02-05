# human motor movement eyeblink script
# analyse each button press for the incidence of CR
# onset relative to button press
# peak latency - likely to be airpuff

# motorBlink function
source("functions/motorBlink.R")

# setwd to folder containing second order files
setwd(choose.dir())
home_dir <- getwd()

# function to apply in sapply for times of stim
stimTime <- function(x){
  as.numeric(strsplit(x, "\t")[[1]][1])
}

# list all the files to analyse
files <- list.files()

# empty data frame for response classification
response_df <- data.frame()

# subject dataframe
sub_df <- as.data.frame(matrix(ncol=15, nrow=0))
colnames(sub_df) <- c("ID", "group", "recording.type","n.airpuff", "n.button", "n.A", "n.B", "n.C", "n.D",
                      "n.BTest", "n.CTest", "avg.RT", "min.RT", "max.RT", "avg.ITI")

# trial waveform dataframe
trial_wv_df <- as.data.frame(matrix(nrow=0, ncol=7))
colnames(trial_wv_df) <- c("time", "value", "response", "id", "group", "homepad", "trial.n")

# iterate through the files
for(f in files){
  tmp_sub_df <- as.data.frame(matrix(ncol=15, nrow=1))
  colnames(tmp_sub_df) <- c("ID", "group", "recording.type","n.airpuff", "n.button", "n.A", "n.B", "n.C", "n.D",
                            "n.BTest", "n.CTest", "avg.RT", "min.RT", "max.RT", "avg.ITI")
  
  # cat which files the script is on 
  cat(paste0("working on ", which(f==files), " out of ", length(files), "\n"))
  
  # extract information from the filename
  file_info <- strsplit(f, split = ".", fixed = T)[[1]][1]
  file_info <- unlist(strsplit(file_info, split = "_"))
  session_type <- file_info[1]
  date <- file_info[2]
  id_tag <- paste0(file_info[2], file_info[3])
  recording_type <- gsub(file_info[4], "MT", "Mech Tra")
  
  # fill tmp_sub_df with id and group(experimental/control)
  tmp_sub_df$ID <- id_tag
  tmp_sub_df$group <- session_type
  tmp_sub_df$recording.type <- recording_type
  
  # read in waveform
  start_point <- which(grepl(recording_type, readLines(f)))[2] + 4
  # positions of "CHANNEL", waveform goes to next one - 2
  chan_pos <- which(grepl("CHANNEL", readLines(f)))
  # which one is the right index to stop at
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  # read in waveform
  waveform <- as.numeric(readLines(f)[start_point : next_chan_pos])
  
  # create a waveform dataframe
  waveform_data <- data.frame(time = seq(0, length.out = length(waveform), by = 0.001), waveform = waveform)
  
  # homepad leaving time vector
  # home pad is mentioned twice, second one is 2 lines before start of time vector 
  start_point <- which(grepl("Home Pad", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  homepad_time <- as.numeric(readLines(f)[start_point : next_chan_pos])
  
  # read in the keyboard comment vector
  start_point <- which(grepl("Keyboard", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  keyboard_time <- character(0)
  for(i in start_point : next_chan_pos){
    time <- strsplit(readLines(f)[i], "\t")[[1]][1]
    letter <- strsplit(strsplit(readLines(f)[i], "\t")[[1]][2], "")[[1]][2]
    if(letter %in% LETTERS){
      keyboard_time <- append(keyboard_time, list(c(time, letter)))
    }
  }
  
  # read in vector of US to get an idea of the size UR and number of airpuffs overall
  start_point <- which(grepl("DigMark", readLines(f)))[2] + 2
  if(is.na(start_point)){
    start_point <- which(grepl("US", readLines(f)))[2] + 2
  }
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  us_time <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  
  # a button time vector(
  start_point <- which(grepl("A button", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  a_time <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  
  # b button time vector
  start_point <- which(grepl("B button", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  b_time <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  
  # c button
  start_point <- which(grepl("C button", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  c_time <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  
  # d button
  start_point <- which(grepl("D button", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  d_time <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  
  # b test
  start_point <- which(grepl("B test", readLines(f)))[2] + 2
  next_chan_pos <- min(chan_pos[chan_pos > start_point]) - 2
  b_test <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  b_test <- na.omit(b_test)
  
  # c test
  start_point <- which(grepl("C test", readLines(f)))[2] + 2
  next_chan_pos <- length(readLines(f))
  c_test <- unlist(lapply(readLines(f)[start_point : next_chan_pos], stimTime))
  c_test <- na.omit(c_test)
  
  # button dataframe
  button_df <- data.frame(button = c(rep("A", length(a_time)), rep("B", length(b_time)), rep("C", length(c_time)),
                                          rep("D", length(d_time)), rep("BTest", length(b_test)), rep("CTest", length(c_test))),
                           time = c(a_time, b_time, c_time, d_time, b_test, c_test))
  button_df <- button_df[order(button_df$time), ]
  button_df$press.N <- seq(1, nrow(button_df), by = 1)
  button_df$RT <-  NA
  for(i in 1:nrow(button_df)){
    button_df$RT[i] <- button_df$time[i] - max(homepad_time[homepad_time < button_df$time[i]])
    if(button_df$RT[i] > 3){
      button_df$RT[i] <- Inf
    }
  }
  button_df$id <- id_tag
  button_df$group <- session_type
  
  # number of airpuff and sub data
  tmp_sub_df$n.airpuff <- length(us_time)
  tmp_sub_df$n.button <- length(c(a_time, b_time, c_time, d_time, b_test, c_test))
  tmp_sub_df$n.A <- length(a_time)
  tmp_sub_df$n.B <- length(b_time)
  tmp_sub_df$n.C <- length(c_time)
  tmp_sub_df$n.D <- length(d_time)
  tmp_sub_df$n.BTest <- length(b_test)
  tmp_sub_df$n.CTest <- length(c_test)
  tmp_sub_df$avg.ITI <- abs(mean(button_df$time[1 : nrow(button_df)-1] - button_df$time[2 : nrow(button_df)]))
  tmp_sub_df$avg.RT <- mean(button_df$RT[!is.infinite(button_df$RT)])
  tmp_sub_df$min.RT <- min(button_df$RT[!is.infinite(button_df$RT)])
  tmp_sub_df$max.RT <- max(button_df$RT[!is.infinite(button_df$RT)])
  
  # pass dataframes to function to determine CR incidence
  dataframes <- motorBlink(waveform = waveform_data, button = button_df, us = us_time)
  
  # bind dataframes
  response_df <- rbind(response_df, dataframes[[1]])
  
  sub_df <- rbind(sub_df, tmp_sub_df)
  
  trial_wv_df <- rbind(trial_wv_df, dataframes[[2]])
}

response_df$reinforced.cut <- cut(response_df$press.N, breaks = c(1,20,120,140,200,220,260), include.lowest = T)
response_df$blocks <- cut(response_df$press.N, breaks = c(1, seq(20, 260, 20)), include.lowest = T)

trial_wv_df$reinforced.cut <- cut(trial_wv_df$trial.n, breaks = c(1,20,120,140,200,220,260), include.lowest = T)
trial_wv_df$blocks <- cut(trial_wv_df$trial.n, breaks = c(1, seq(20, 260, 20)), include.lowest = T)

setwd("..")

write.csv(response_df, file = "motor eyeblink conditioning responses.csv", row.names = F)
write.csv(trial_wv_df, file = "trial waveform.csv", row.names = F)

