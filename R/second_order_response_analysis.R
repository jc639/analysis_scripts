# bring all coherent data together
# CSVs for all trial responses in one directory
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(lme4)
library(nlme)
# standard error function
se <- function(x){
  x <- na.omit(x)
  sqrt(var(x)/length(x))
} 

# choose that directory
setwd(choose.dir())

# list files to iterate over and add day identifier
csv <- list.files()

trial_response_df <- as.data.frame(matrix(nrow=0, ncol=11))
colnames(trial_response_df) <- c("response", "trial.type", "subject", "trial.n", "session.n", "date", "onset",
                                 "peak", "amplitude", "day", "session")
# iterate through CSVs
for(f in csv){
  temp_df <- read.csv(f)
  unique_days <- unique(temp_df$date)
  temp_df$day <- sapply(temp_df$date, function(x) which(unique_days==x))
  temp_df$session <- ifelse(temp_df$subject=="Ben", temp_df$day, temp_df$session)
  trial_response_df <- rbind(trial_response_df, temp_df)
}

# unique subjects
uni_sub <- unique(trial_response_df$subject)

# Ben also has CS1
if("Ben"%in%uni_sub){
  trial_response_df <- trial_response_df[trial_response_df$trial.type!="CS1",]
}

# simple onset diagram of latencies of stimuli, just a geom rect with x axis above
onset_diagram <- ggplot()+geom_rect(aes(xmin = -0.35, xmax=0.06, ymin=-1.2, ymax=-0.2), fill="#619CFF", alpha=0.8)+
  geom_rect(aes(xmin=0, xmax=0.42, ymax=-1.2, ymin=-2.2), fill= "#F8766D", alpha=0.8)+
  geom_rect(aes(xmin=0.35, xmax=0.42, ymin=-3.3, ymax=-2.2), fill="grey", alpha=0.8)+
  theme_bw() + theme(plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border= element_blank(),
                    axis.line.x = element_line(color="black"),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())+
  geom_text(aes(x=-0.14, y=-0.7,  label="CS2"), size=6) +
  geom_text(aes(x=0.21, y=-1.7, label="CS1"), size=6) +
  geom_text(aes(x=0.385, y=-2.7, label="US"), size=6) +
  xlab("onset relative to CS1 (ms)") +
  scale_x_continuous(limits = c(-0.4, 0.42), breaks = c(-0.35, 0, 0.35), labels = function(x)x*1000, position = "top",
                     expand=c(0,0))+scale_y_continuous(expand=c(0,0))

# count of responses
response_counts <- as.data.frame(matrix(nrow=0, ncol=4))
colnames(response_counts) <- c("trial.type", "response", "count", "subject")
total_trial_n <- as.data.frame(matrix(nrow=0, ncol=3))
colnames(total_trial_n) <- c("trial.type", "total", "subject")

# counts of early responses, chi-sq test against same period before (350ms-alpha)CS1-US vs CS2CS1
# make histograms of onsets of all trial types overlaid
# filter by subject
for(s in uni_sub){
  df <- trial_response_df %>%
  filter(subject == s)
  
  # group by by trial type and repsonse
  tmp_resp_count <- df %>%
    group_by(trial.type, response) %>%
    summarise(count=n())
  tmp_resp_count$subject <- s
  
  # indicator of blink onset before CS1 starts
  df$cs2on.before <- ifelse(df$trial.type=="CS2CS1" & df$onset < 0, 1, 0)
  # remove spurious onsets of those classed as not flat for hist
  df$onset <- ifelse(df$response=="not flat baseline" & is.infinite(df$onset), NA, df$onset)
  # indicator of not flat before CS1 in period same as CS2-CS1
  df$cs1on.before <- ifelse(df$trial.type=="CS1US" &
                              df$response=="not flat baseline" & df$onset > 0.185, 1, 0)
  # bind counts into data frame
  tmp_resp_count$response <- as.character(tmp_resp_count$response)
  tmp_resp_count <- rbind(as.data.frame(tmp_resp_count), c("CS2CS1", "early onset", sum(df$cs2on.before, na.rm = T), s))
  tmp_resp_count <- rbind(tmp_resp_count, c("CS1US", "early onset", sum(df$cs1on.before, na.rm = T), s))
  tmp_resp_count$count <- as.numeric(tmp_resp_count$count)
  
  # bind temporary into response counts data frame
  response_counts <- rbind(response_counts, tmp_resp_count)
  
  # count df
  tmp_tot_count <- df %>%
    group_by(trial.type) %>%
    summarise(total=n())
  tmp_tot_count$subject <- s
  total_trial_n <- rbind(total_trial_n, as.data.frame(tmp_tot_count))
  
  # count of early for CS1, CS2
  cs1early <- sum(df$cs1on.before, na.rm = T)
  cs2early <- sum(df$cs2on.before, na.rm = T)
  cs1total <- tmp_tot_count$total[tmp_tot_count$trial.type=="CS1US"] - cs1early
  cs2total <- tmp_tot_count$total[tmp_tot_count$trial.type=="CS2CS1"] - cs2early
  
  # Chi matrix
  chi_m <- as.table(rbind(c(cs1total, cs2total), c(cs1early, cs2early)))
  dimnames(chi_m) <- list(timing=c("Not Early", "Early"),
                          trial=c("CS1-US", "CS2-CS1"))
  # chisq test, no expected <5, so not Fisher
  Xsq <- chisq.test(chi_m)
  chi_table <- tableGrob(chi_m)
  ggsave(paste0(s, " - counts of early onsets.pdf"), chi_table)
  
  # results in text file
  sink(file = paste0(s, " - results of chi sq test of counts of early onset.txt"))
  print(Xsq)
  sink()
  
  # histogram of onset latency and peak latency for each subject
  resp_types <- c("onset", "peak")
  for(r in resp_types){
    
    # filter to cs1 df and only CR, no CR
    cs1_df <- filter(df, trial.type=="CS1US", response%in%c("CR", "NO CR"))
    cs1_hist <- ggplot(cs1_df, aes(x=cs1_df[r]))+geom_histogram(binwidth = 0.02, fill="#F8766D", colour="#F8766D", alpha=0.8)+
      theme_classic() + scale_y_continuous(expand=c(0,0)) +  scale_x_continuous(limits = c(-0.4, 0.42), breaks = c(-0.35, 0, 0.35),
                                                                                labels = function(x)x*1000, expand = c(0,0))+
      xlab("") + ggtitle("CS1-US trials")
    
    # cs2cs1 filter
    cs2cs1_df <- filter(df, trial.type=="CS2CS1", response%in%c("CR", "NO CR"))
    cs2cs1_hist <- ggplot(cs2cs1_df, aes(x=cs2cs1_df[r]))+geom_histogram(binwidth = 0.02, fill="#619CFF", colour="#619CFF", alpha=0.8)+
      theme_classic() + scale_y_continuous(expand=c(0,0)) +  scale_x_continuous(limits = c(-0.4, 0.42), breaks = c(-0.35, 0, 0.35),
                                                                                labels = function(x)x*1000, expand = c(0,0))+
      xlab("") + ggtitle("CS2-CS1 trials")
    
    # cs2 filter
    cs2_df <- filter(df, trial.type=="CS2", response%in%c("CR", "NO CR"))
    cs2_hist <- ggplot(cs2_df, aes(x=cs2_df[r]))+geom_histogram(binwidth = 0.02, fill="#00BA38", colour="#00BA38", alpha=0.8)+
      theme_classic() + scale_y_continuous(expand=c(0,0)) +  scale_x_continuous(limits = c(-0.4, 0.42), breaks = c(-0.35, 0, 0.35),
                                                                                labels = function(x)x*1000, expand = c(0,0))+
      xlab(paste0("blink ", r, " relative to CS1 (ms)")) + ggtitle("CS2 test trials")
    
    # arrange into a grid
    grid_plot <- arrangeGrob(onset_diagram, cs1_hist, cs2cs1_hist, cs2_hist, 
                             ncol=1, heights=unit(c(1.5,3,3,3), rep("in", 4)), widths=unit(6, "in"))
    ggsave(paste0(r," latency histograms for ", s, ".pdf"), grid_plot, width=8, height = 12, units = "in")
  }
}

# onset_df for latency dot plots
onset_df <- trial_response_df %>%
  filter(trial.type%in%c("CS1US", "CS2CS1"), response=="CR")

# filter by subject
for(s in uni_sub){
  df <- onset_df %>%
    filter(subject==s)
  
  # others have consistently 25 trials in a session, Ben has 30
  if (s != "Ben"){
    cs1_df <- filter(df, trial.type=="CS1US", trial.n <= 25)
    cs2_df <- filter(df, trial.type=="CS2CS1", trial.n <= 25)
  } else {
    cs1_df <- filter(df, trial.type=="CS1US", trial.n <= 30)
    
    cs2_df <- filter(df, trial.type=="CS2CS1", trial.n <= 30)
  }
  
  # cs1 dotplot
  cs1_dotplot <- ggplot(cs1_df, aes(x=trial.n, y=onset)) + geom_point(colour="#F8766D", alpha=0.5) + theme_classic() +
    scale_y_continuous(limits = c(-0.36,0.36), labels = function(y)y*1000, breaks = c(-0.35, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS1-US")
  
  # cs2cs1
  cs2_dotplot <- ggplot(cs2_df, aes(x=trial.n, y=onset)) + geom_point(colour="#619CFF", alpha=0.5) + theme_classic() +
    scale_y_continuous(limits = c(-0.36,0.36), labels = function(y)y*1000, breaks = c(-0.35, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS2-CS1")

  # probabilty of early response in CS2CS1
  prob_df <- cs2_df %>%
    group_by(trial.n) %>%
    summarise(total=n(),
              early=length(onset[onset < 0]),
              prob = early/total)
  
  # plot of the probabilty of early onset for trial numbers within session  
  prob_plot <- ggplot(prob_df, aes(x=trial.n, y=prob)) + geom_line(colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("trial number within a session")
  
  # add an indicator column to cs2_df, where onset < 0 = 1, and onset > 0 = 0
  cs2_df$early <- ifelse(cs2_df$onset < 0, 1, 0)
  
  # create single variable logistic regression and predictions
  single_var_glm <- glm(early~trial.n, cs2_df, family = binomial(link = "logit"))
  prob_df$pred <- predict(single_var_glm, newdata = prob_df[, c("trial.n")], type = "response")
  
  # create prediction plot based on logistic regression line and observed probabilities
  pred_plot <- ggplot(prob_df) + geom_point(aes(x=trial.n, y=prob), colour="#619CFF", size=1.5) + 
    geom_line(aes(x=trial.n, y=pred), colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10)) +
    scale_y_continuous(limits=c(0,1)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("trial number within a session")

  # arrange the plots into a grid, 1 column and 3 rows - this one has prob_plot
  grid_plot <- arrangeGrob(cs1_dotplot, cs2_dotplot, prob_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  
  # save the grid plot
  ggsave(filename = paste0(s, " - onset latency for trial numbers within sessions1.pdf"), grid_plot,
         height = 12, width = 8, units = "in")
  
  # same as before but this one contains pred_plot
  grid_plot <- arrangeGrob(cs1_dotplot, cs2_dotplot, pred_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  ggsave(filename = paste0(s, " - onset latency for trial numbers within sessions2.pdf"), grid_plot,
         height = 12, width = 8, units = "in")
  
  # put the trial number glm stats into a text file
  sink(paste0(s, " - trial number single variable glm on early onset.txt"))
  print(summary(single_var_glm))
  sink()
  
  # same but for session
  cs1_dotplot <- ggplot(cs1_df, aes(x=session, y=onset)) + geom_point(colour="#F8766D", alpha=0.5) + theme_classic() +
    scale_y_continuous(limits = c(-0.36,0.36), labels = function(y)y*1000, breaks = c(-0.35, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS1-US")
  
  # cs2 plot
  cs2_dotplot <- ggplot(cs2_df, aes(x=session, y=onset)) + geom_point(colour="#619CFF", alpha=0.5) + theme_classic() +
    scale_y_continuous(limits = c(-0.36,0.36), labels = function(y)y*1000, breaks = c(-0.35, 0, 0.35)) +
    scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5)) +
    xlab("") + ylab("blink onset relative to CS1 (ms)") +
    geom_hline(aes(yintercept=0),linetype=2) + ggtitle("CS2-CS1")
  
  # probability of early onset by session
  prob_df <- cs2_df %>%
    group_by(session) %>%
    summarise(total=n(),
              early=length(onset[onset < 0]),
              prob = early/total)
  
  # probability line plot for session onset for cs2cs1
  prob_plot <- ggplot(prob_df, aes(x=session, y=prob)) + geom_line(colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("session number")
  
  # single variable - session - glm logistic regression 
  single_var_glm <- glm(early~session, data=cs2_df, family = binomial(link = "logit"))
  prob_df$pred <- predict(single_var_glm, newdata = prob_df[, c("session")], type = "response")
  
  # prediction plot
  pred_plot <- ggplot(prob_df) + geom_point(aes(x=session, y=prob), colour="#619CFF", size=1.5) + 
    geom_line(aes(x=session, y=pred), colour="#619CFF", size=1) +
    theme_classic() + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5)) +
    ylab("Probability of early onset (< 0 ms)") + xlab("session number") + scale_y_continuous(limits = c(0,1))
  
  # grid arrange with prob_plot first
  grid_plot <- arrangeGrob(cs1_dotplot, cs2_dotplot, prob_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  
  ggsave(filename = paste0(s, " - onset latency for sessions1.pdf"), grid_plot,
         height = 12, width = 8, units = "in")
  
  # same but with prediction 
  grid_plot <- arrangeGrob(cs1_dotplot, cs2_dotplot, pred_plot, 
                           ncol=1, heights=unit(c(3,3,3), rep("in", 4)), widths=unit(6, "in"))
  
  ggsave(filename = paste0(s, " - onset latency for sessions2.pdf"), grid_plot,
         height = 12, width = 8, units = "in")
  
  # save the session glm in a text
  sink(paste0(s, " - session number single variable glm on early onset.txt"))
  print(summary(single_var_glm))
  sink()

}

# cs2 onset df for GLM, convert < 0 onset -> 1(=success)
cs2_onset <- onset_df %>%
  filter(trial.type=="CS2CS1")
cs2_onset$early <- ifelse(cs2_onset$onset < 0, 1, 0)

# taken the probability as whole, don't think this is the correct approach
if(!"Ben"%in%uni_sub){
  simple_early_glm <- glm(early~trial.n+session+subject, data = cs2_onset, family = binomial(link = "logit"))
  
  int_early_glm <- glm(early~trial.n*session*subject, data = cs2_onset, family = binomial(link = "logit"))
  
  random_eff_glm <- glmer(early~trial.n+session+(1|subject), data=cs2_onset, family = binomial(link = "logit"))
  
  sink("GLM results on probability of early onset.txt")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  cat("\tSIMPLE ADDITIVE BINOMIAL GLM OF EARLY ONSET IN CSCS1 TRIALS\n")
  cat(paste0(rep("*", 80), collapse = ""),"\n")
  print(summary(simple_early_glm))
  cat("\n\n\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  cat("\tBINOMIAL GLM OF EARLY ONSET IN CSCS1 TRIALS WITH INTERACTIONS\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  print(summary(int_early_glm))
  cat("\n\n\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  cat("\tMIXED MODEL BINOMIAL GLM OF EARLY ONSET IN CSCS1 TRIALS\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  print(summary(random_eff_glm))
  sink()
  
  # for each subject iterate through and produce plots that take into account interaction between the session and trial n
  for (s in uni_sub){
    simple_sub_glm <- glm(early~trial.n+session, data = cs2_onset[cs2_onset$subject==s,], family = binomial(link = "logit"))
    int_sub_glm <- glm(early~trial.n*session, data = cs2_onset[cs2_onset$subject==s,], family = binomial(link = "logit"))
    
    df <- onset_df %>%
      filter(subject==s)
    
    # others have consistently 25 trials in a session, Ben has 30
    if (s != "Ben"){
      cs1_df <- filter(df, trial.type=="CS1US", trial.n <= 25)
      cs2_df <- filter(df, trial.type=="CS2CS1", trial.n <= 25)
    } else {
      cs1_df <- filter(df, trial.type=="CS1US", trial.n <= 30)
      
      cs2_df <- filter(df, trial.type=="CS2CS1", trial.n <= 30)
    }
    
    # save the results of the glm
    sink(paste0(s, " - GLM of early response probabilities.txt"))
    print(summary(simple_sub_glm))
    cat(paste(rep("*", 80), collapse = ""))
    cat("\nInteraction terms added\n")
    print(summary(int_sub_glm))
    sink()
    
    # dataframe to pass to predict 
    trial_new_data <- data.frame(trial.n = rep(1:25,2), session = c(rep(1,25), rep(25,25)))
    session_new_data <-data.frame(trial.n = c(rep(1,25), rep(25, 25)), session = rep(1:25, 2))
    
    # add predictions
    trial_new_data$simple <- predict(simple_sub_glm, trial_new_data, type = "response")
    trial_new_data$int <- predict(int_sub_glm, trial_new_data[,c(1,2)], type="response")
    trial_new_data$session.n <- c(rep("1st session", 25), rep("25th session", 25))
    
    # add predictions
    session_new_data$simple <- predict(simple_sub_glm, session_new_data, type="response")
    session_new_data$int <- predict(int_sub_glm, session_new_data[,c(1,2)], type="response")
    session_new_data$trial.number <- c(rep("1st trial", 25), rep("25th trial", 25)) 
    
    # generate plots, simple plots
    pred_simple_session <- ggplot(session_new_data, aes(x=session, y=simple, group=trial.number)) + 
      geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=trial.number), size=3, colour="#619CFF") +
      xlab("session number") + ylab("predicted probability of early response") + guides(shape=guide_legend("trial number")) +
      theme_classic() + ggtitle("simple glm") + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5))
    pred_simple_trial <- ggplot(trial_new_data, aes(x=trial.n, y=simple, group=session.n)) + 
      geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=session.n), size=3, colour="#619CFF") +
      xlab("trial number") + ylab("predicted probability of early response") + guides(shape=guide_legend("session number")) +
      theme_classic() + ggtitle("simple glm") + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10))
    
    # plots for interaction glm
    pred_int_session <- ggplot(session_new_data, aes(x=session, y=int, group=trial.number)) + 
      geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=trial.number), size=3, colour="#619CFF") +
      xlab("session number") + ylab("predicted probability of early response") + guides(shape=guide_legend("trial number")) +
      theme_classic() + ggtitle("interaction glm") + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5))
    pred_int_trial <- ggplot(trial_new_data, aes(x=trial.n, y=int, group=session.n)) + 
      geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=session.n), size=3, colour="#619CFF") +
      xlab("trial number") + ylab("predicted probability of early response") + guides(shape=guide_legend("session number")) +
      theme_classic() + ggtitle("interaction glm") + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10))
    
    # save plots as pdf
    ggsave(pred_simple_session, filename = paste0(s, " - simple predicted probabilities for early onset in session.pdf"), height = 4, width = 6)
    ggsave(pred_int_session, filename = paste0(s, " - interaction predicted probabilities for early onset in session.pdf"), height = 4, width = 6)
    ggsave(pred_simple_trial, filename = paste0(s, " - simple predicted probabilities for early onset in trial.pdf"), height = 4, width = 6)
    ggsave(pred_int_trial, filename = paste0(s, " - interaction predicted probabilities for early onset in trial.pdf"), height = 4, width = 6)
  }
} else {
  # same but for Ben
  simple_early_glm <- glm(early~trial.n+session, data = cs2_onset, family = binomial(link = "logit"))
  
  int_early_glm <- glm(early~trial.n*session, data = cs2_onset, family = binomial(link = "logit"))
  
  sink("GLM results on probability of early onset.txt")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  cat("\tSIMPLE ADDITIVE BINOMIAL GLM OF EARLY ONSET IN CSCS1 TRIALS\n")
  cat(paste0(rep("*", 80), collapse = ""),"\n")
  print(summary(simple_early_glm))
  cat("\n\n\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  cat("\tBINOMIAL GLM OF EARLY ONSET IN CSCS1 TRIALS WITH INTERACTIONS\n")
  cat(paste0(rep("*", 80), collapse = ""), "\n")
  print(summary(int_early_glm))
  sink()
  
  trial_new_data <- data.frame(trial.n = rep(1:30,2), session = c(rep(1,30), rep(25,30)))
  session_new_data <-data.frame(trial.n = c(rep(1,30), rep(30, 30)), session = rep(1:30, 2))
  
  trial_new_data$simple <- predict(simple_early_glm, trial_new_data, type = "response")
  trial_new_data$int <- predict(int_early_glm, trial_new_data[,c(1,2)], type="response")
  trial_new_data$session.n <- c(rep("1st session", 30), rep("25th session", 30))
  
  session_new_data$simple <- predict(simple_early_glm, session_new_data, type="response")
  session_new_data$int <- predict(int_early_glm, session_new_data[,c(1,2)], type="response")
  session_new_data$trial.number <- c(rep("1st trial", 30), rep("25th trial", 30)) 
  
  pred_simple_session <- ggplot(session_new_data, aes(x=session, y=simple, group=trial.number)) + 
    geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=trial.number), size=3, colour="#619CFF") +
    xlab("session number") + ylab("predicted probability of early response") + guides(shape=guide_legend("trial number")) +
    theme_classic() + ggtitle("simple glm") + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5))
  pred_simple_trial <- ggplot(trial_new_data, aes(x=trial.n, y=simple, group=session.n)) + 
    geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=session.n), size=3, colour="#619CFF") +
    xlab("trial number") + ylab("predicted probability of early response") + guides(shape=guide_legend("session number")) +
    theme_classic() + ggtitle("simple glm") + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10))
  
  pred_int_session <- ggplot(session_new_data, aes(x=session, y=int, group=trial.number)) + 
    geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=trial.number), size=3, colour="#619CFF") +
    xlab("session number") + ylab("predicted probability of early response") + guides(shape=guide_legend("trial number")) +
    theme_classic() + ggtitle("interaction glm") + scale_x_continuous(limits = c(1, max(cs1_df$session)), breaks=seq(5, max(cs1_df$session), 5))
  pred_int_trial <- ggplot(trial_new_data, aes(x=trial.n, y=int, group=session.n)) + 
    geom_line(colour="#619CFF", size=1) + geom_point(aes(shape=session.n), size=3, colour="#619CFF") +
    xlab("trial number") + ylab("predicted probability of early response") + guides(shape=guide_legend("session number")) +
    theme_classic() + ggtitle("interaction glm") + scale_x_continuous(limits = c(1, max(cs1_df$trial.n)), breaks=seq(10, max(cs1_df$trial.n), 10))
  
  ggsave(pred_simple_session, filename = paste0(uni_sub, " - simple predicted probabilities for early onset in session.pdf"), height = 4, width = 6)
  ggsave(pred_int_session, filename = paste0(uni_sub, " - interaction predicted probabilities for early onset in session.pdf"), height = 4, width = 6)
  ggsave(pred_simple_trial, filename = paste0(uni_sub, " - simple predicted probabilities for early onset in trial.pdf"), height = 4, width = 6)
  ggsave(pred_int_trial, filename = paste0(uni_sub, " - interaction predicted probabilities for early onset in trial.pdf"), height = 4, width = 6)
}

# show as grouped data
group_onset_sess <- onset_df %>%
  group_by(session, trial.type) %>%
  summarise(avg=mean(onset),
            se=se(onset))

if(!"Ben"%in%uni_sub){
  group_onset_trial <- onset_df %>%
    filter(trial.n <= 25) %>%
    group_by(trial.n, trial.type) %>%
    summarise(avg=mean(onset),
              se=se(onset))
} else {
  group_onset_trial <- onset_df %>%
    group_by(trial.n, trial.type) %>%
    summarise(avg=mean(onset),
              se=se(onset))
}
 # grouped data
sess_plot <- ggplot(group_onset_sess, aes(x=session, y=avg, colour=trial.type)) + geom_point() +
  geom_errorbar(aes(x=session, ymax=avg+se, ymin=avg-se, colour=trial.type), width=0.25) +
  geom_line() + theme_classic() + scale_x_continuous(limits=c(0.5,max(group_onset_sess$session)+0.5),
                                                     breaks=seq(10, max(group_onset_sess$session),10)) +
  scale_y_continuous(limits=c(-0.36, 0.36), breaks=c(-0.35, 0, 0.35), labels=function(x)x*1000) +
  geom_hline(aes(yintercept=0), linetype=2) + ylab("blink onset relative to CS1 (ms)") +
  xlab("session number") + scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=F) +
  ggtitle("All subjects")



trial_plot <- ggplot(group_onset_trial, aes(x=trial.n, y=avg, colour=trial.type)) + geom_point() +
  geom_errorbar(aes(x=trial.n, ymax=avg+se, ymin=avg-se, colour=trial.type), width=0.25) +
  geom_line() + theme_classic() + scale_x_continuous(limits=c(0.5,max(group_onset_trial$trial.n)+0.5),
                                                     breaks=seq(5, max(group_onset_trial$trial.n),5)) +
  scale_y_continuous(limits=c(-0.36, 0.36), breaks=c(-0.35, 0, 0.35), labels=function(x)x*1000) +
  geom_hline(aes(yintercept=0), linetype=2) + ylab("blink onset relative to CS1 (ms)") +
  xlab("trial number within sessions") + scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) +
  ggtitle("")

# save
grid_plot <- arrangeGrob(sess_plot, trial_plot,  
                         ncol=2, heights=unit(3,"in"), widths=unit(c(6,6), c("in", "in")))
ggsave(filename="grouped onset latency.pdf", grid_plot, height=5, width=15)

# same but for each session
if(length(uni_sub)!=1){
  for(s in uni_sub){
    onset_sess <- onset_df %>%
      filter(subject==s) %>%
      group_by(session, trial.type) %>%
      summarise(avg=mean(onset),
                se=se(onset))
    
    onset_trial <- onset_df %>%
      filter(subject==s, trial.n <= 25) %>%
      group_by(trial.n, trial.type) %>%
      summarise(avg=mean(onset),
                se=se(onset))
    
    sub_sess_plot <- ggplot(onset_sess, aes(x=session, y=avg, colour=trial.type)) + geom_point() +
      geom_errorbar(aes(x=session, ymax=avg+se, ymin=avg-se, colour=trial.type), width=0.25) +
      geom_line() + theme_classic() + scale_x_continuous(limits=c(0.5,max(onset_sess$session)+0.5),
                                                         breaks=seq(10, max(onset_sess$session),10)) +
      scale_y_continuous(limits=c(-0.36, 0.36), breaks=c(-0.35, 0, 0.35), labels=function(x)x*1000) +
      geom_hline(aes(yintercept=0), linetype=2) + ylab("blink onset relative to CS1 (ms)") +
      xlab("session number") + scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=F) +
      ggtitle(paste0("Subject ", which(uni_sub==s)))
    
    
    sub_trial_plot <- ggplot(onset_trial, aes(x=trial.n, y=avg, colour=trial.type)) + geom_point() +
      geom_errorbar(aes(x=trial.n, ymax=avg+se, ymin=avg-se, colour=trial.type), width=0.25) +
      geom_line() + theme_classic() + scale_x_continuous(limits=c(0.5,max(onset_trial$trial.n)+0.5),
                                                         breaks=seq(5, max(onset_trial$trial.n),5)) +
      scale_y_continuous(limits=c(-0.36, 0.36), breaks=c(-0.35, 0, 0.35), labels=function(x)x*1000) +
      geom_hline(aes(yintercept=0), linetype=2) + ylab("blink onset relative to CS1 (ms)") +
      xlab("trial number within sessions") + scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) +
      ggtitle("")
    
    grid_plot <- arrangeGrob(sub_sess_plot, sub_trial_plot,  
                             ncol=2, heights=unit(3,"in"), widths=unit(c(6,6), c("in", "in")))
    
    ggsave(filename = paste0(s, " - average onset latency line graph.pdf"), grid_plot, height=5, width=15)
  }
}

# alternative to glm, anova and linear model
onset_anova <- aov(onset~trial.type+as.factor(trial.n)+as.factor(session), data=onset_df)

random_anova <- anova(lme(onset~trial.type+as.factor(trial.n), random = ~1|subject, data=onset_df))

onset_lm <- lm(onset~trial.type+trial.n+session, data=onset_df)

random_lm <- lme(onset~trial.type + trial.n + session, random=~1|subject, data = onset_df)

# anova text file
sink(file = "ANOVA AND LM RESULTS ON ONSET.txt")
cat(paste0(rep("*", 80), collapse = ""), "\n")
cat("\tANOVA OF DIFFERENCE IN MEANS FOR VARIOUS FACTORS\n")
cat(paste0(rep("*", 80), collapse = ""),"\n")
print(summary(onset_anova))
cat("\n\n\n")
cat(paste0(rep("*", 80), collapse = ""), "\n")
cat("\tRANDOM EFFECTS ANOVA TO ACCOUNT FOR SUBJECT\n")
cat(paste0(rep("*", 80), collapse = ""),"\n")
print(random_anova)
cat("\n\n\n")
cat(paste0(rep("*", 80), collapse = ""), "\n")
cat("\tLINEAR REGRESSION FOR ONSET\n")
cat(paste0(rep("*", 80), collapse = ""),"\n")
print(summary(onset_lm))
cat("\n\n\n")
cat(paste0(rep("*", 80), collapse = ""), "\n")
cat("\tRANDOM EFFECTS LINEAR REGRESSION TO ACCOUNT FOR SUBJECT\n")
cat(paste0(rep("*", 80), collapse = ""),"\n")
print(summary(random_lm))
cat("\n\n\n")
sink()

# add 1/0 to trial response df
trial_response_df$success <- ifelse(trial_response_df$response=="CR", 1, 0)
# now plots for trial success
CRincidence_df <- trial_response_df %>%
  filter(response%in%c("CR", "NO CR")) %>%
  group_by(subject, session, trial.type) %>%
    summarise(success = length(response[response=="CR"]),
              total_trials = n(),
              percentage = (success/total_trials) * 100)

# Tarquin has one session that has only one good trial =100%, more likely that he didnt so adjusted to be lower
CRincidence_df$percentage <- ifelse(CRincidence_df$total_trials == 1, (CRincidence_df$success/5)*100, CRincidence_df$percentage)
for(s in uni_sub){
  sub_cr_sess <- filter(CRincidence_df, subject==s, trial.type=="CS2") %>%
    select(subject, session, percentage)
  sub_cr_trial <- trial_response_df %>%
    filter(response%in%c("CR", "NO CR"), subject==s, trial.type=="CS2", trial.n <= ifelse(s=="Ben", 6, 5)) %>%
    group_by(trial.n) %>%
    summarise(success = length(response[response=="CR"]),
              total_trials = n(),
              percentage = (success/total_trials) * 100)
  
  # CS2 response probability by trial n
  single_trial_glm <- glm(success~trial.n, 
                          data = filter(trial_response_df, subject==s, response%in%c("CR", "NO CR"), 
                                        trial.type=="CS2"), family = binomial(link = "logit"))
  # CS2 response probability by session
  single_sess_glm <- glm(success~session, 
                         data = filter(trial_response_df, subject==s, response%in%c("CR", "NO CR"), 
                                       trial.type=="CS2"), family = binomial(link = "logit"))
  # predictions
  sub_cr_sess$pred <- predict(single_sess_glm, newdata = sub_cr_sess[, c("session")], type = "response")*100
  sub_cr_trial$pred <- predict(single_trial_glm, newdata = sub_cr_trial[, c("trial.n")], type = "response")*100
  
  # plots with observed and predicted
  pred_cr_sessplot <- ggplot(sub_cr_sess) + geom_point(aes(x=session, y=percentage), colour="#00BA38", size=1.5) + 
    geom_line(aes(x=session, y=pred), colour="#00BA38") + scale_y_continuous(limits=c(0, 100)) +
    theme_classic() + ylab("percentage of CRs") + 
    scale_x_continuous(limits = c(1, max(sub_cr_sess$session)), breaks=seq(5, max(sub_cr_sess$session), 5))
  
  # plots for observed and predicted for trial n
  pred_cr_trialplot <- ggplot(sub_cr_trial) + geom_point(aes(x=trial.n, y=percentage), colour="#00BA38", size=1.5) + 
    geom_line(aes(x=trial.n, y=pred), colour="#00BA38", size=1) + scale_y_continuous(limits=c(0, 100)) +
    theme_classic() + ylab("percentage of CRs") + xlab("trial number") +
    scale_x_continuous(limits = c(1, max(sub_cr_trial$trial.n)), breaks=seq(1, max(sub_cr_trial$trial.n), 1))
  
  # save
  ggsave(pred_cr_sessplot, filename = paste0(s, " - CRs on CS2 test by session and predicted.pdf"),
         width = 6, height = 4)
  ggsave(pred_cr_trialplot, filename = paste0(s, " - CRs on CS2 test by trial and predicted.pdf"),
         width = 6, height = 4)
  sink(paste0(s, " - single variable glm for CS2 test CRs.txt"))
  print(summary(single_sess_glm))
  cat("\n\n\n\n\n")
  print(summary(single_trial_glm))
  sink()
}

# block the percentage of responses
block_start <- seq(1, max(CRincidence_df$session), 3)
block_end <- seq(3, max(CRincidence_df$session), 3)
block_start <- block_start[1:min(c(length(block_start), length(block_end)))]
block_end <- block_end[1:min(c(length(block_start), length(block_end)))]

blocked_df <- as.data.frame(matrix(ncol=6, nrow=length(block_start)*length(unique(CRincidence_df$trial.type))*length(uni_sub)))
colnames(blocked_df) <- c("subject", "block", "trial.type", "subtrial", "avgP", "se")

j <- 1
for(s in uni_sub){
  for(i in 1:length(block_start)){
    df <- CRincidence_df %>%
      filter(subject==s, session%in%block_start[i]:block_end[i])
    for(tt in unique(df$trial.type)){
      blocked_df$subject[j] <- s
      blocked_df$trial.type[j] <- tt
      blocked_df$block[j] <- i
      blocked_df$subtrial[j] <- paste0(s, tt)
      blocked_df$avgP[j] <- mean(df$percentage[df$trial.type==tt], na.rm=T)
      blocked_df$se[j] <- se(df$percentage[df$trial.type==tt])
      j <- j + 1
    }
  }
}

blocked_df <- blocked_df[!is.na(blocked_df$subject),]
for(i in 1:nrow(blocked_df)){
  blocked_df$subNUM[i] <- paste0("Subject ", which(blocked_df$subject[i]==uni_sub))
}
# plots of CR percentage by block
perc_plot <- ggplot(blocked_df, aes(x=block, y=avgP)) +
  geom_line(aes(group=subtrial, colour=trial.type), alpha=0.6, size=0.8) +
  geom_errorbar(aes(x=block, ymax=avgP+se, ymin=avgP-se, colour=trial.type), width=0.25, alpha=0.6, size=0.8) +
  theme_classic() + facet_grid(.~subNUM) + ylab("percentage of CRs") +
  xlab("block number") + scale_colour_manual(values = c("#F8766D", "#00BA38","#619CFF"), guide=guide_legend(title = "trial type"))

for(subs in  unique(blocked_df$subNUM)){
  sub_perc_plot <- ggplot(filter(blocked_df, subNUM==subs), aes(x=block, y=avgP)) +
    geom_line(aes(group=subtrial, colour=trial.type), alpha=0.6, size=0.8) +
    geom_errorbar(aes(x=block, ymax=avgP+se, ymin=avgP-se, colour=trial.type), width=0.25, alpha=0.6, size=0.8) +
    theme_classic() + ylab("percentage of CRs") +
    xlab("block number") + scale_colour_manual(values = c("#F8766D", "#00BA38","#619CFF"), guide=guide_legend(title = "trial type"))
  ggsave(filename = paste0(subs, " - percentage of CRs in blocks by trial type.pdf"),  width=6, height=4)
}

if(!"Ben"%in%uni_sub){
  ggsave(filename = "percentage of CRs in each trial type.pdf", perc_plot, width=4*length(uni_sub), height=1*length(uni_sub))
}




# glm of CS2 test
for(s in uni_sub){
  simple_CS2_glm <- glm(success~trial.n+session, data=filter(trial_response_df, response%in%c("CR", "NO CR"),
                                                      trial.type=="CS2", subject==s),
                 family = binomial(link = "logit"))
  int_CS2_glm <- glm(success~trial.n*session, data=filter(trial_response_df, response%in%c("CR", "NO CR"),
                                                             trial.type=="CS2", subject==s),
                        family = binomial(link = "logit"))
  sink(paste0(s, " - glm on CS2 responses over session and trial number.txt"))
  print(summary(simple_CS2_glm))
  cat("\n")
  cat(paste(rep("*",80), collapse = ""))
  cat("\nInteraction terms model\n")
  print(summary(int_CS2_glm))
  sink()
}

# add cuts to trial response
trial_response_df$cut <- cut(trial_response_df$session, breaks=seq(0, max(trial_response_df$session), 3))
amplitude_df <- trial_response_df %>%
  filter(trial.type%in%c("CS2CS1", "CS1US"), response=="CR", trial.n <= ifelse("Ben"%in%uni_sub, 30, 25)) %>%
  group_by(subject, cut, trial.type) %>%
  summarise(avg_amp = mean(amplitude, na.rm=T),
            se_amp = se(amplitude))

# remove na so sessions outside cut, so can iterate through block
amplitude_df <- na.omit(amplitude_df)
amplitude_df$block <- NA
for(i in 1:nrow(amplitude_df)){
  amplitude_df$block[i] <- which(unique(amplitude_df$cut) == amplitude_df$cut[i])
}

# response amplitude
for(s in uni_sub){
  sub_amp <- filter(amplitude_df, subject==s)
  amplitude_plot <- ggplot(sub_amp, aes(x=block, y=avg_amp, colour=trial.type)) + geom_point() +
    geom_errorbar(aes(x=block, ymax=avg_amp+se_amp, ymin=avg_amp-se_amp, colour=trial.type), width=0.25) +
    geom_line(size=1) + theme_classic() +
    scale_y_continuous(limits=c(min(sub_amp$avg_amp)-max(sub_amp$se_amp)+.1, max(sub_amp$avg_amp)+max(sub_amp$se_amp)+.1), breaks=seq(round(min(sub_amp$avg_amp)), round(max(sub_amp$avg_amp)))) + 
    ylab("blink amplitude (mm)") + xlab("block number") + scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type"))
  ggsave(amplitude_plot, filename = paste0(s, " - amplitude of response by session.pdf"), width=6, height = 4)
}


#generate glm for CS1US and CS2CS1 responses and also linear models for ampltiude
for(s in uni_sub){
  sub_df <- trial_response_df %>%
    filter(subject == s, trial.type%in%c("CS1US", "CS2CS1"), response%in%c("CR", "NO CR"))
 
  perc_df <- sub_df %>%
    group_by(session, trial.type) %>%
    summarise(success = length(response[response=="CR"]),
              total_trials = n(),
              percentage = (success/total_trials) * 100)
  
  perc_glm_simple <- glm(success~trial.type+session, data = sub_df, family = binomial(link = "logit"))
  perc_glm_int <- glm(success~trial.type*session, data = sub_df, family = binomial(link = "logit"))
  amplitude_lm <- lm(amplitude~trial.type+session, data = sub_df)
  int_amplitude_lm <- lm(amplitude~trial.type*session, data=sub_df)
  
  perc_df$simple <- predict(perc_glm_simple, newdata = perc_df[, c("session", "trial.type")], type = "response") * 100
  perc_df$int <- predict(perc_glm_int, newdata = perc_df[, c("session", "trial.type")], type = "response") * 100
  
  amp_df <- data.frame(session = perc_df[, c("session")], trial.type = perc_df[, c("trial.type")],
                       simple = predict(amplitude_lm, newdata = perc_df[, c("session", "trial.type")], type = "response"),
                       int = predict(int_amplitude_lm, newdata = perc_df[, c("session", "trial.type")], type = "response"))
  
  sink(paste0(s , ' - glm test of primary CRs and linear model of amplitude.txt'))
  cat("\n")
  print(summary(perc_glm_simple))
  cat("\n\n\n\n\n\n")
  print(summary(perc_glm_int))
  cat(paste(rep("*", 80), collapse = ""))
  cat("\nLM RESULTS\n")
  print(summary(amplitude_lm))
  cat("\n\n\n\n\n")
  print(summary(int_amplitude_lm))
  sink()
  
  simple_glm_plot <- ggplot(perc_df) + geom_point(aes(x = session, y = percentage, colour = trial.type), alpha=0.6) +
    geom_line(aes(x = session, y = simple, colour = trial.type), size=1) + theme_classic() +
    scale_y_continuous(limits=c(0, 100)) + ylab("percentage of CRs") +
    scale_x_continuous(limits=c(0.5,max(perc_df$session)+0.5), breaks=seq(10, max(perc_df$session),10)) +
    scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type"))
  
  int_glm_plot <- ggplot(perc_df) + geom_point(aes(x = session, y = percentage, colour = trial.type), alpha=0.6) +
    geom_line(aes(x = session, y = int, colour = trial.type), size=1) + theme_classic() +
    scale_y_continuous(limits=c(0, 100)) + ylab("percentage of CRs") +
    scale_x_continuous(limits=c(0.5,max(perc_df$session)+0.5), breaks=seq(10, max(perc_df$session),10)) +
    scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type"))
  
  lm_plot <- ggplot() + geom_point(data = sub_df, aes(x = session, y = amplitude, colour=trial.type), alpha=0.5) +
    geom_line(data = amp_df, aes(x=session, y=simple, colour=trial.type), size=1) +
    scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) + ylab("blink amplitude (mm)") +
    scale_x_continuous(limits=c(0.5,max(perc_df$session)+0.5), breaks=seq(10, max(perc_df$session),10)) + theme_classic() +
    facet_grid(.~trial.type)
  
  int_lm_plot <- ggplot() + geom_point(data = sub_df, aes(x = session, y = amplitude, colour=trial.type), alpha=0.5) +
    geom_line(data = amp_df, aes(x=session, y=int, colour=trial.type), size=1) +
    scale_colour_manual(values = c("#F8766D", "#619CFF"), guide=guide_legend(title = "trial type")) + ylab("blink amplitude (mm)") +
    scale_x_continuous(limits=c(0.5,max(perc_df$session)+0.5), breaks=seq(10, max(perc_df$session),10)) + theme_classic() +
    facet_grid(.~trial.type)
  
  ggsave(simple_glm_plot, filename = paste0(s, " - CS1US CS2CS1 simple glm observed and predict.pdf"),
         width=6, height = 4)
  ggsave(int_glm_plot, filename = paste0(s, " - CS1US CS2CS1 interaction glm observed and predict.pdf"),
         width=6, height = 4)
  ggsave(lm_plot, filename = paste0(s, " - CS1US CS2CS1 simple amplitude lm observed and predict.pdf"),
         width=6, height = 4)
  ggsave(int_lm_plot, filename = paste0(s, " - CS1US CS2CS1 interaction amplitude lm observed and predict.pdf"),
         width=6, height = 4)
}

