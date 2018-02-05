# assess whether cs2 blinks occur
# change in CR onset in CS2-CS1 trials
# glm of CS2 prop
# averages of responses on each trial type
# dataframe of trial responses, latency, type etc
# dataframe of waveform aligned to CS2 onset: time, value, trial type, date, trial n, response category etc
source("functions/blinkResponse.R")

# setwd to folder containing second order files
setwd(choose.dir())
home_dir <- getwd()

# list of sub directories, contains all the text files
folders <- list.dirs(recursive = F)

for (f in folders){
  # set to the first folder of text files, each folder should paired conditioning protocol, can add unpaired later
  setwd(f)
  cat(paste0("working on ", f, "\n"))
  
  #list files to choose NMR and markers
  files <- list.files()
  
  #need to determine how many files for this session
  number_of_files <- NULL
  for (i in 1:length(files)){
    number_of_files <- append(number_of_files, strsplit(files[i], "_")[[1]][5])
  }
  
  #determines the number of unique file00N.txt
  unique_files <- length(unique(number_of_files))
  
  for (file_n in 1:unique_files){
    #find correct files for that session and get filenames
    correct_files <- files[grepl(paste0("file00", file_n, ".txt"), files)]
    
    NMRfile <- correct_files[grepl("NMR", correct_files)]
    
    CS1file <- correct_files[grepl("CS1", correct_files)]
    
    CS1testfile <- correct_files[grepl("primarytest", correct_files)]
    
    CS2alonefile <- correct_files[grepl("CS2alone", correct_files)]
    
    CS2pairedfile <- correct_files[grepl("CS2paired", correct_files)]
    
    #set date and subject
    date <- strsplit(correct_files[1], "_")[[1]][1]
    subject <- strsplit(correct_files[1], "_")[[1]][4]
    
    #point in NMRfile that has "START" and also temporal sampling rate
    start_point <- which(grepl("START", readLines(NMRfile)))
    sample_rate <- as.numeric(strsplit(readLines(NMRfile)[start_point], split = "\t" )[[1]][3])
    
    #generate NMR vector to deal with gaps and remember start times
    NM <- readLines(NMRfile)[(start_point + 1):length(readLines(NMRfile))]
    NM <- NM[NM!=""]
    
    gaps <- NM[grepl("GAP", NM)]
    
    #if no gaps then length of grepl should be zero
    if(length(gaps) != 0){
      ref_table <- as.data.frame(matrix(ncol = 2,nrow = length(gaps)))
      colnames(ref_table) <- c("Index","Time")
      
      
      for(i in 1:length(gaps)){
        ref_table$Index[i] <- which(NM==gaps[i])
        ref_table$Time[i] <- as.numeric(strsplit(gaps[i], split="\t")[[1]][2])
        NM<-NM[NM!=gaps[i]]
      }
      
      NM <- as.numeric(NM)
      
      waveform_data <- as.data.frame(matrix(nrow=length(NM),ncol=2))
      colnames(waveform_data) <- c("NMR.position", "Time")
      
      waveform_data$NMR.position <- NM
      
      ref_table <- rbind(c(1, 0), ref_table)
      ref_table <- rbind(ref_table, c(nrow(waveform_data), NA))
      
      gap_start <- numeric(0)
      gap_end <- numeric(0)
      for(i in 1:(nrow(ref_table) - 1)){
        if (i == (nrow(ref_table) - 1)){
          waveform_data$Time[ref_table$Index[i] : ref_table$Index[i + 1]] <- seq(ref_table$Time[i], 
                                                                               length.out = nrow(waveform_data[ref_table$Index[i] : ref_table$Index[i+1],]),
                                                                               by = sample_rate)
        } else {
          waveform_data$Time[ref_table$Index[i]:(ref_table$Index[i+1]-1)] <- seq(ref_table$Time[i], 
                                                                                                 length.out = ref_table$Index[i+1]-ref_table$Index[i],
                                                                                                 by=sample_rate)
          gap_start <- append(gap_start, waveform_data$Time[ref_table$Index[i+1]-1])
          gap_end <- append(gap_end, ref_table$Time[i+1])
        }
      }
    } else {
      NM <- as.numeric(NM)
      
      waveform_data <- as.data.frame(matrix(nrow=length(NM),ncol=2))
      colnames(waveform_data) <- c("NMR.position","Time")
      
      waveform_data$NMR.position <- NM
      waveform_data$Time <- seq(0, sample_rate*(nrow(waveform_data) - 1), sample_rate)		
    }
    
	# some stimulus markers are in gaps of sampling --> remove
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
      
    #read in marker data, first the CS1
    CS1US.marker <- NULL
    for (i in 17:(length(readLines(CS1file)))){
      CS1US.marker <- append(CS1US.marker, as.numeric(strsplit(readLines(CS1file)[i], "\t")[[1]][1]))
    }
    
    #remove any NAs
    CS1US.marker <- CS1US.marker[!is.na(CS1US.marker)]
    # remove any outside waveform time
    if(length(gaps) != 0){
      CS1US.marker <- inGaps(markers = CS1US.marker, start = gap_start, end = gap_end)
    }
    
    #now read in CS2 marker data
    CS2CS1.marker <- NULL
    for (i in 17:(length(readLines(CS2pairedfile)))){
      CS2CS1.marker <- append(CS2CS1.marker, as.numeric(strsplit(readLines(CS2pairedfile)[i], "\t")[[1]][1]))
    }
    
    CS2CS1.marker <- CS2CS1.marker[!is.na(CS2CS1.marker)]
    if(length(gaps) != 0){
      CS2CS1.marker <- inGaps(markers = CS2CS1.marker, start = gap_start, end = gap_end)
    }
    
    #now read in CS2 alone marker data
    CS2alone.marker <- NULL
    for(i in 17:(length(readLines(CS2alonefile)))){
      CS2alone.marker <- append(CS2alone.marker, as.numeric(strsplit(readLines(CS2alonefile)[i], "\t")[[1]][1]))
    }
    CS2alone.marker <- CS2alone.marker[!is.na(CS2alone.marker)]
    if(length(gaps) != 0){
      CS2alone.marker <- inGaps(markers = CS2alone.marker, start = gap_start, end = gap_end)
    }
    
    #now read in CS1test marker, if there is one, isn't for first three
    if(length(CS1testfile) != 0){
      CS1test.marker <- NULL
      for(i in 17:(length(readLines(CS1testfile)))){
        CS1test.marker <- append(CS1test.marker, as.numeric(strsplit(readLines(CS1testfile)[i],"\t")[[1]][1]))
      }
      CS1test.marker <- CS1test.marker[!is.na(CS1test.marker)]
      if(length(gaps) != 0){
        CS1test.marker <- inGaps(markers = CS1test.marker, start = gap_start, end = gap_end)
      }
    }
    
    # make dataframes for each marker
    # CS1US
    dfs <- blinkDataframes(markers = CS1US.marker, waveform = waveform_data, trial.type = "CS1US",
                           subject = subject, date = date, session.n = file_n, sample.rate = sample_rate)
    
    if(exists("trial_response_df")){
      trial_response_df <- rbind(trial_response_df, dfs[[1]])
      
      trial_waveform_df <- rbind(trial_waveform_df, dfs[[2]])
      
    } else{
      trial_response_df <- dfs[[1]]
      
      trial_waveform_df <- dfs[[2]]  
    }
    
    # CS2CS1
    dfs <- blinkDataframes(markers = CS2CS1.marker, waveform = waveform_data, trial.type = "CS2CS1",
                          subject = subject, date = date, session.n = file_n, sample.rate = sample_rate)
    
    trial_response_df <- rbind(trial_response_df, dfs[[1]])
    
    trial_waveform_df <- rbind(trial_waveform_df, dfs[[2]])
    
    # CS2 test trials
    dfs <- blinkDataframes(markers = CS2alone.marker, waveform = waveform_data, trial.type= "CS2",
                           subject = subject, date = date, session.n = file_n, sample.rate = sample_rate)
    trial_response_df <- rbind(trial_response_df, dfs[[1]])
    
    trial_waveform_df <- rbind(trial_waveform_df, dfs[[2]])
    
    
    # CS1 tests if exist
    if(length(CS1testfile) != 0){
      dfs <- blinkDataframes(markers = CS1test.marker, waveform = waveform_data, trial.type= "CS1",
                             subject = subject, date = date, session.n = file_n, sample.rate = sample_rate)
      trial_response_df <- rbind(trial_response_df, dfs[[1]])
      
      trial_waveform_df <- rbind(trial_waveform_df, dfs[[2]])
    }
  }
  setwd(home_dir)
}
# save as csv
write.csv(trial_response_df, row.names = F, file = paste0(subject, " - all trials response.csv"))
write.csv(trial_waveform_df, row.names = F, file = paste0(subject, " - all trials waveform.csv"))

rm(list=ls())
