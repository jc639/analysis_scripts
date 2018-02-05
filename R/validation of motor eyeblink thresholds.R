# comparison of different thresholds and standard deviation
# error % and average onset difference

# 3d scatterplot library
library(scatterplot3d)

# set to the directory with the validation
setwd(choose.dir())

# list the files
files <- list.files()

# human check file
validation <- read.csv(grep("check", files, value=T))
validation$response <- gsub("yes", "CR", validation$response)
validation$response <- gsub("no", "No_CR", validation$response)

# remove the validation
files <- files[-(grep("check", files))]

# create df
error_df <- as.data.frame(matrix(nrow=length(files), ncol=4))
colnames(error_df) <- c("threshold", "SD", "correct", "onset_diff")

for(i in 1:length(files)){
  file_split <- strsplit(strsplit(files[i], ".csv")[[1]], " ")[[1]]
  error_df$threshold[i] <- file_split[3]
  error_df$SD[i] <- file_split[7]
  
  test <- read.csv(files[i])
  
  error <- numeric(0)
  onset_diff <- numeric(0)
  
  for(j in 1:nrow(validation)){
    CR <- test$CR[test$id==validation$id[j] & (test$time > validation$time[j] - 0.2) & (test$time < validation$time[j] + 0.2)]
    onset <- test$onset[test$id==validation$id[j] & (test$time > validation$time[j] - 0.2) & (test$time < validation$time[j] + 0.2)]
    if(!length(CR)){
      next
    }
    
    if(CR == "not_flat"){
      next
    }
    error <- c(error, CR==validation$response[j])
    
    if(CR==validation$response[j]){
      onset_diff <- c(onset_diff, onset - validation$onset[j])
    }
  }
  error_df$correct[i] <- mean(error) * 100
  error_df$onset_diff[i] <- mean(onset_diff, na.rm = T)
}

error_df$threshold <- as.numeric(error_df$threshold)
error_df$SD <- as.numeric(error_df$SD)


pdf(file = "correct percentage plot.pdf", width=9, height = 6)
correct_plot <- scatterplot3d(x=error_df$threshold, y=error_df$SD, z=error_df$correct,
                              highlight.3d = T, type = 'h', pch=16,
                              angle=120, xlab="blink threshold (average UR size * threshold)", ylab = "Number of standard deviations from baseline",
                              zlab = "% matching human judgement")
dev.off()

pdf(file = "correct percentage plot2.pdf", width=9, height = 6)
correct_plot2 <- scatterplot3d(x=error_df$threshold, y=error_df$SD, z=error_df$correct,
                               highlight.3d = T, type = 'h', pch=16,
                               angle=40, xlab="blink threshold (average UR size * threshold)", ylab = "Number of standard deviations from baseline",
                               zlab = "% matching human judgement")

dev.off()

pdf(file = "onset plot.pdf", width=9, height = 6)
onset_plot <- scatterplot3d(x=error_df$threshold, y=error_df$SD, z=abs(error_df$onset_diff)*1000,
                            highlight.3d = T, type = 'h', pch=16,
                            angle=120, xlab="blink threshold (average UR size * threshold)", ylab = "Number of standard deviations from baseline",
                            zlab = "absolute average difference in onset (ms)")
dev.off()

pdf(file = "onset plot2.pdf", width=9, height = 6)
onset_plot2 <- scatterplot3d(x=error_df$threshold, y=error_df$SD, z=abs(error_df$onset_diff)*1000,
                           highlight.3d = T, type = 'h', pch=16,
                           angle=40, xlab="blink threshold (average UR size * threshold)", ylab = "Number of standard deviations from baseline",
                           zlab = "absolute average difference in onset (ms)")


dev.off()




