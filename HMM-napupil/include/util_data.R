# causal filter


# 
filter_causal <- function(y, win_s, padded = FALSE) {
  
  N <- length(y)
  
  if (padded) {
    y_smooth <- numeric(N)
    # y_pad <- c(rep(0, win_s), y)
    y_pad <- c(rep(mean(y), win_s), y)
    
    for (t in 1:N) {
      y_smooth[t] <- sum(y_pad[t + (1:win_s)])/win_s; 
    }
  } else {
    N_smooth <- N - win_s + 1
    y_smooth <- c()
    for (t in 1:N_smooth) {
      y_smooth <- c(y_smooth, sum(y[t:(t+win_s-1)])/win_s)
    }
    y_smooth <- c(rep(NA, win_s  - 1), y_smooth)
  }
  
  return(y_smooth)
}

# - 
detrend_linear <- function(x) {
  N = length(x)
  X = matrix(1, nrow = N, ncol = 2)
  X[, 2] =  1:N
  
  b = solve(t(X) %*% X) %*% t(X) %*% x
  trend = X %*% b
  y = x - trend
  
  return(y)
}


scale_by_trials <- function(data) 
{
  
  N_trials <- nrow(data)
  N_samples <- ncol(data)
  
  mu_trials <- apply(data, 1, mean, na.rm = T)
  sd_trials <- apply(data, 1, sd, na.rm = T)
  
  out <- matrix(NA, nrow = N_trials, ncol = N_samples)
  for (n in 1:N_trials) {
    out[n, ] <- (data[n, ] - mu_trials[n])/sd_trials[n]
  }
  return(out)
}



# - auxiliary for plotting 
get_last_NA <- function(win_s, filter_shifted = TRUE, padded = FALSE) {
  if(filter_shifted) {
    if (padded) {
      last_NA <- 0
    } else {
      last_NA <- win_s - 1
      # double check if it's win_s or win_s - 1, non andiamo a fede
    }
  } else {
    last_NA <-  floor(win_s/2)
  }
  return(last_NA)
}



# - ?
smooth_trials_mod <- function(dataset, win_s, 
                              filter_shifted = TRUE, padded = FALSE) 
{
  
  N_trials <- nrow(dataset)
  N_samples <- ncol(dataset)
  out <- matrix(NA, nrow(dataset), ncol(dataset))
  for(i in 1:N_trials) {
    if (any(is.na(dataset[i, ]))) {
      out[i, ] <-rep(NA, N_samples)
    } else {
      if (filter_shifted) {
        out[i, ] <- filter_causal(dataset[i, ], win_s, padded)
      } else {
        out[i, ] <- filter(dataset[i, ], rep(1/win_s, win_s), method = "convolution")
      }
    }
  }
  return(out)
}



get_neuromodulators <- function(ID, info_vars, snippet = "stimulus", 
                                substrat = "allTrials", 
                                groups = NULL, win_s = 5, filter_shifted = TRUE, 
                                padded = FALSE, detrend = FALSE, plt = T) 
{
  info_id <- info_vars$id
  info_id_aux <- info_vars$id_aux
  disease <- info_vars$disease
  disease_cat <- disease[info_id == ID]
  
  path_nm_temp <- paste(getwd(), "data", snippet, "nm/", sep="/")
  path_behav <- paste(getwd(), "data", snippet, "behavior/", sep="/")
  
  # start & end 
  N_all = 61
  
  if (filter_shifted) {
    start = ifelse(padded, 1, win_s)
    end = N_all
  } else {
    if (win_s>0) {
      start = 1 + floor(win_s/2)
      end = N_all - floor(win_s/2)
    } else {
      start = 1; end = N_all;
    }
  }
  
  N <- length(start:end)
  
  if (ID == 123) {
    IDs_tmp = c(1, 2, 3)
    out <- get_data_allSubjects(IDs_tmp, snippet = snippet)
    includeTrials <- out$includeTrials
    idx_includeTrials <- which(includeTrials == 1)
    NE_all <- out$NE[idx_includeTrials, ]
    Pupil_all <- out$Pupil[idx_includeTrials, ]
    Oddball <- out$Oddball[idx_includeTrials, ]
    Valence <- out$Valence[idx_includeTrials, ]
    Arousal <- out$Arousal[idx_includeTrials, ]

    # do something (get_data_allSubjects)
  } else {
    # --- load data
    includeTrials <- as.vector(readMat(paste(path_nm_temp, "ID", 
                                             ID, "_includeTrials.mat", sep =""))$includeTrials)
    idx_includeTrials <- which(includeTrials == 1)
    
    # nm
    NE_all <- readMat(paste(path_nm_temp, "ID", 
                            ID, "_NE.mat", sep =""))$NE[idx_includeTrials, ]
    Pupil_all <- readMat(paste(path_nm_temp, "ID", 
                               ID, "_pupil.mat", sep =""))$pupil[idx_includeTrials, ]
    
    #  behevioral
    Oddball <- readMat(paste(path_behav, "ID", 
                             ID, "_Oddball.mat", sep =""))$Oddball[idx_includeTrials, ]
    Valence <- readMat(paste(path_behav, "ID", 
                             ID, "_Valence.mat", sep =""))$Valence[idx_includeTrials, ]
    Arousal <- readMat(paste(path_behav, "ID", 
                             ID, "_Arousal.mat", sep =""))$Arousal[idx_includeTrials, ]
    
  }
  

  # - behavioral - relabeling
  Oddball[which(Oddball == 0)] = "Checkerboard"
  Oddball[which(Oddball == 1)] = "Oddball"
  
  Valence[which(Valence == -1)] = "Negative"
  Valence[which(Valence == 0)] = "Neutral"
  Valence[which(Valence == 1)] = "Positive"
  
  Arousal[which(Arousal == -1)] = "Low"
  Arousal[which(Arousal == 0)] = "Neutral"
  Arousal[which(Arousal == 1)] = "High"
  
  # smoothing each trials NE
  if (win_s > 0) {
    NE_all <- smooth_trials_mod(NE_all, win_s, filter_shifted, padded)
    Pupil_all <- smooth_trials_mod(Pupil_all, win_s, filter_shifted, padded)
  }
  
  # scaling by trials (necessary when comparing across experiments)
  NE_all <- scale_by_trials(NE_all)
  Pupil_all <- scale_by_trials(Pupil_all)
  
  # excluding NAs (trials)
  include <- rep(FALSE, dim(NE_all)[1])
  idx_include <- unique(which(!is.na(as.data.frame(NE_all)),
                              arr.ind = T)[, 1])
  include[idx_include] <- TRUE
  
  
  # number of trials
  ntrials_cat <- length(which(include == TRUE))
  
  
  if (substrat != "allTrials") {
    
    include_sub <- list()
    n_groups <- length(groups)
    
    if (substrat == "Oddball") {
      temp <- Oddball
      for (jj in 1:n_groups) {
        include_sub[[jj]] <- (temp == groups[jj] & include)
      }
    } else if (substrat == "Valence") {
      temp <- Valence
      for (jj in 1:n_groups) {
        include_sub[[jj]] <- (temp == groups[jj] & include)
      }
    } else if (substrat == "Arousal") {
      temp <- Arousal
      for (jj in 1:n_groups) {
        include_sub[[jj]] <- (temp == groups[jj] & include)
      }
    }
    
    
    # cross designed experiment
    if (substrat == "Oddball_x_Arousal") {
      include_sub[[1]] <- (Oddball == "Oddball") & (Arousal == "Low") & include
      #include_sub[[2]] <- (Oddball == "Oddball") & (Arousal == "Neutral") & include
      include_sub[[2]] <- (Oddball == "Oddball") & (Arousal == "High") & include
    } 
    
    # - substratifying observations
    obs_groups <- vector(mode="list")
    n_trials_groups <- numeric(n_groups)
    
    # - means groups
    for (jj in 1:n_groups) {
      n_trials_groups[jj] <- length(which(include_sub[[jj]] == TRUE))
      obs_groups[[jj]] <- matrix(NA, nrow = N, ncol = 2)
      obs_groups[[jj]][, 1] <- colMeans(NE_all[include_sub[[jj]], start:end])
      obs_groups[[jj]][, 2] <- colMeans(Pupil_all[include_sub[[jj]], start:end], na.rm = T)
    }
    
    
    # detrend 
    if (detrend) {
      for (jj in 1:n_groups) {
        for (ii in 1:2) {
          obs_groups[[jj]][, ii] = detrend_linear(obs_groups[[jj]][, ii])
        }
      }
    }
    
    # - obtaining common scale
    sd_scale <- rep(0, 2) # one for each neurotransimmeter (or pupil)? 
    for(jj in 1:n_groups) {
      sd_scale <- sd_scale + apply(obs_groups[[jj]], 2, sd, na.rm = T)
    }
    sd_scale <- sd_scale/n_groups
    for (jj in 1:n_groups) {
      for (ii in 1:2) {
        obs_groups[[jj]][, ii] = obs_groups[[jj]][, ii]/sd_scale[ii]
      }
    }
    
    
  } else {
    
    obs <- matrix(NA, nrow = N, ncol = 2)
    obs[, 1] <- colMeans(NE_all[include, start:end])
    obs[, 2] <- colMeans(Pupil_all[include, start:end], na.rm = T)
    
    # detrend 
    if (detrend) {
      for (ii in 1:2) {
        obs[, ii] = detrend_linear(obs[, ii])
      }
    }
    
    # scale it 
    obs[, 1] <- scale(obs[, 1])[,]
    obs[, 2] <- scale(obs[, 2])[,]
    
    
  }
  
  if (plt) {
    
    # - figure parameters 
    pars_fig <- list()
    pars_fig$cex = 1.5
    pars_fig$pch = 20
    ticks <- seq(from = 1, to = N_all, by = 10)
    x <- round(tail(ticks, 1)/10)/2
    labels <- as.character(seq(from=-x, to = x, 
                               len = length(ticks)))
    line_stimulus <- floor(N_all/2) + 1
    # last_NA <- ifelse(causal, 0, floor(win_s/2))
    last_NA <- get_last_NA(win_s, filter_shifted, padded)
    
    if (substrat != "allTrials") {
      
      # - figure parameters 
      pars_fig <- list()
      pars_fig$cex = 1.5
      pars_fig$pch = 20
      ticks <- seq(from = 1, to = N_all, by = 10)
      x <- round(tail(ticks, 1)/10)/2
      labels <- as.character(seq(from=-x, to = x, 
                                 len = length(ticks)))
      line_stimulus <- floor(N_all/2) + 1
      
      
      # if substrats != "allTrials"
      
      # - 
      par(mfrow = c(1, n_groups), mar = c(5,7,4,2) + 0.1, 
          mai = c(0.6, 0.35, 0.5, 0.25))
      
      # - plot groups
      for(jj in 1:n_groups) {
        
        plot(1:N + last_NA, obs_groups[[jj]][, 1], cex = pars_fig$cex, 
             pch = pars_fig$pch, type = "o", 
             xaxt = "n",
             xlim = c(1, N_all),
             ylab = "Neurotransmitters", xlab = "Time (seconds)",
             col = "blue", ylim = c(-3.5, 3.5))
        lines(1:N + last_NA, obs_groups[[jj]][, 2],cex = pars_fig$cex, 
              pch = pars_fig$pch,col = "black", lwd = 2)
        abline(v = line_stimulus, lwd = 1, lty = "dotted")
        abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
        axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
        legend("bottomleft", pch = c(rep(pars_fig$pch, 1
                                         ), NA),
               legend = c("NE", "Pupil"), 
               col = c("blue", "black"), cex = 0.6, 
               lty = c(rep(NA, 1), 1))
        title(paste(groups[jj],  " (n: ", n_trials_groups[jj], ")", 
                    sep = ""), line = 0.5)
      }
      
      mtext(paste(ifelse(snippet == "stimulus", "Stimulus, ", "Stimulus, "),
                  "ID: ", info_id_aux[which(info_id == ID)],  
                  "  (by ", substrat, ")               ", 
                  #get_label_smoothing(win_s, filter_shifted, padded),
                  sep = ""), side = 3, line = -1.2, outer = T, 
            cex = 0.7)
      
    } else {
      # plot all trials 
      par(mfrow = c(1, 1))
      plot(1:N + last_NA, obs[, 1], cex = pars_fig$cex, 
           pch = pars_fig$pch, type = "o", 
           xaxt = "n",
           xlim = c(1, N_all),
           ylab = "Neurotransmitters", xlab = "Time (seconds)",
           col = "blue", ylim = c(-3.5, 3.5))
      lines(1:N + last_NA, obs[, 2],cex = pars_fig$cex, 
            pch = pars_fig$pch,col = "black", lwd = 2)
      abline(v = line_stimulus, lwd = 1, lty = "dotted")
      abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
      axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
      legend("bottomleft", pch = c(rep(pars_fig$pch, 1), NA),
             legend = c("NE", "Pupil"), 
             col = c("blue", "black"), cex = 0.6, 
             lty = c(rep(NA, 1), 1))
      
      title(paste(ifelse(snippet == "stimulus", "Stimulus, ", "Stimulus, "),
                  "ID: ", info_id_aux[which(info_id == ID)],  
                  "  (All Trials)               ",
                  #get_label_smoothing(win_s, filter_shifted, padded),
                  sep = ""), line = -1.2, outer = T, 
            cex = 0.7)
      
    }
  }
  # return observations
  if (substrat != "allTrials") {
    return(obs_groups)
  } else {
    return(obs)
  }
}


get_data_allSubjects <- function(IDs = c(1, 2, 3), snippet = "stimulus")
{
  path_nm_temp <- paste(getwd(), "data", snippet, "nm/", sep="/")
  path_behav <- paste(getwd(), "data", snippet, "behavior/", sep="/")
  
  # - neuromodulators
  
  includeTrials_subjects <- as.vector(readMat(paste(path_nm_temp, "ID", 
                                                    1, "_includeTrials.mat", sep =""))$includeTrials)
  NE_all_subjects <- readMat(paste(path_nm_temp, "ID", 
                                   1, "_NE.mat", sep =""))$NE
  Pupil_all_subjects <- readMat(paste(path_nm_temp, "ID", 
                                      1, "_pupil.mat", sep =""))$pupil
  
  
  for (ii in 2:length(IDs)) {
    includeTrials_temp <- as.vector(readMat(paste(path_nm_temp, "ID", 
                                                  IDs[ii],
                                                  "_includeTrials.mat", sep =""))$includeTrials)
    NE_temp <- readMat(paste(path_nm_temp, "ID", 
                             IDs[ii], "_NE.mat", sep =""))$NE
    Pupil_temp <- readMat(paste(path_nm_temp, "ID", 
                                IDs[ii], "_pupil.mat", sep =""))$pupil
    
    includeTrials_subjects <- c(includeTrials_subjects, includeTrials_temp)
    NE_all_subjects <- rbind(NE_all_subjects, NE_temp)
    Pupil_all_subjects <- rbind(Pupil_all_subjects, Pupil_temp)
  }
  
  # - behavioral
  Oddball_subjects <- readMat(paste(path_behav, "ID", 
                                    1, "_Oddball.mat", sep =""))$Oddball
  Valence_subjects <- readMat(paste(path_behav, "ID", 
                                    1, "_Valence.mat", sep =""))$Valence
  Arousal_subjects <- readMat(paste(path_behav, "ID", 
                                    1, "_Arousal.mat", sep =""))$Arousal 
  
  for (ii in 2:length(IDs)) {
    Oddball_subjects <- rbind(Oddball_subjects, readMat(paste(path_behav, "ID", 
                                                              IDs[ii], "_Oddball.mat", sep =""))$Oddball)
    Valence_subjects <- rbind(Valence_subjects, readMat(paste(path_behav, "ID", 
                                                              IDs[ii], "_Valence.mat", sep =""))$Valence)
    Arousal_subjects <- rbind(Arousal_subjects, readMat(paste(path_behav, "ID", 
                                                              IDs[ii], "_Arousal.mat", sep =""))$Arousal)
  }
  
  out <- list()
  
  out$includeTrials <- includeTrials_subjects
  out$NE <- NE_all_subjects
  out$Pupil <- Pupil_all_subjects
  out$Oddball <- Oddball_subjects
  out$Valence <- Valence_subjects
  out$Arousal <- Arousal_subjects
  
  return(out)
}

