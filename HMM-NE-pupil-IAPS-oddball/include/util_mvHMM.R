N_all = 61

c_light = "gray95"
c_light_highlight = "gray90"
c_mid = "gray85"
c_mid_highlight = "gray80"
c_dark = "gray75"

cols_grey = c("gray95", "gray90", "gray85", "gray80", "gray75")
cols_green = c("palegreen", "palegreen1", "palegreen2", "palegreen3", 
               "springgreen4")
cols_salmon = c("lightsalmon", "lightsalmon1", "lightsalmon2", "lightsalmon3", 
                "lightsalmon4")

cols_grey = c("gray95", "gray90", "gray85", "gray80", "gray75")
cols_green = c("palegreen", "palegreen1", "palegreen2", "palegreen3", 
               "springgreen4")
cols_salmon = c("lightsalmon", "lightsalmon1", "lightsalmon2", "lightsalmon3", 
                "lightsalmon4")
cols_gold = c("gold", "gold1", "gold2", "gold3", "gold4")

cols_all = list()
cols_all[[1]] <- cols_green 
cols_all[[2]] <- cols_salmon
cols_all[[3]] <- cols_grey
cols_all[[4]] <- cols_gold


# - 
mvHMM_performance <- function(fit, obs, K_max = 4) 
{
  sims <- rstan::extract(fit)
  N <- nrow(obs)
  
  # collect samples
  mu_sims <- sims$mu
  tau_sims <- sims$tau
  L_Omega_sims <- sims$L_Omega
  Sigma_sims <- get_Sigma_sims(sims)
  gamma_sims <- sims$gamma
  corr_sims <- get_corr_sims(sims)
  
  # bayes estimates
  n_MCMC <- dim(mu_sims)[1]
  K <- dim(mu_sims)[2]
  D <- ncol(obs)
  mu_hat <- lapply(1:K, function(kk) {
    apply(mu_sims[, kk, ],  2, mean)})
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  gamma_hat <- apply(gamma_sims,  c(2, 3), mean)
  corr_hat <- lapply(1:K, function(kk) {mean(corr_sims[[kk]])})
  
  
  # quantiles
  mu_info <- lapply(1:K, function(kk) {
    round(apply(mu_sims[, kk, ],  2, quantile, probs = c(.05, 0.5, .95)), 2)})
  corr_info <- lapply(1:K, function(kk) {round(quantile(corr_sims[[kk]], 
                                                        probs = c(.05, 0.5, .95)), 2)})
  # (fill empty; auxilary due to plot)
  if (K != K_max) {
    for (jj in (K + 1):K_max) {
      mu_info[[jj]] <- matrix("", 3, 2)
      corr_info[[jj]] <- rep("", 3)
    }
  }
  
  
  allprobs <- mvGauss_emission(obs, mu_hat, Sigma_hat) 
  
  llk <- mvHMM_forward_backwards(obs, allprobs, gamma_hat)$llk
  np <- (K-1)*K + K*D + K*((D*(D+1))/2)
  
  AIC <- round(-2*(llk - np), 2)
  BIC <- round(-2*llk+np*log(N), 2)
  logml <- round(bridge_sampler(fit, silent = T)$logml, 2)
  
  
  bayes_est <- list(mu = mu_hat, Sigma = Sigma_hat, 
                    gamma = gamma_hat, 
                    corr = corr_hat)
  posterior_quantiles <- list(mu = mu_info, 
                              corr = corr_info)
  
  return(list(AIC = AIC, BIC = BIC, logml = -logml,
              bayes_est = bayes_est,
              posterior_quantiles = posterior_quantiles )) 
}




# ? 
mvHMM_forward_backwards <- function(obs, allprobs, gamma)
{
  K <- nrow(gamma)
  N <- nrow(obs)
  
  lalpha <- lbeta <- matrix(NA, K, N)
  gamma_init <- solve(t(diag(K)-gamma+1),rep(1,K))
  foo <- gamma_init*allprobs[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha[, 1] <- log(foo) + lscale
  
  for (i in 2:N) {
    foo <- foo%*%gamma*allprobs[i,]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[, i] <- log(foo) +lscale
  }
  llk <- lscale
  
  lbeta[, N] <- rep(0, K)
  foo <- rep(1/K, K)
  lscale <- log(K)
  
  for (i in (N-1):1) {
    foo <- gamma%*%(allprobs[i+1,]*foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  list(la = lalpha, lb = lbeta, llk = llk)
}



# (stan) - retrieve Sigma sims from L_Omega and tau 
# (see stan users guide Section 1.13 Mvt Priors..)
get_Sigma_sims <- function(sims)
{
  L_Omega.sims <- sims$L_Omega  
  tau.sims <- sims$tau
  n_MCMC <- dim(L_Omega.sims)[1]
  K <-  dim(L_Omega.sims)[2]
  D <- dim(L_Omega.sims)[3]
  
  Sigma.sims <- array(NA, dim = c(n_MCMC, K, D, D))
  
  for (tt in 1:n_MCMC) {
    L_O.all <- L_Omega.sims[tt, , ,]
    tau.all <- tau.sims[tt, , ]
    for (kk in 1:K) {
      Omega <- L_O.all[kk, , ] %*% t(L_O.all[kk, , ])
      tau <- tau.all[kk, ]
      Sigma.sims[tt, kk, ,] <- quad_form_diag(Omega, tau)
    }
  }
  return(Sigma.sims)
}

# ?
get_corr_sims <- function(sims)
{
  
  Sigma_sims <- get_Sigma_sims(sims)
  n_MCMC <- dim(Sigma_sims)[1]
  K <- dim(Sigma_sims)[2]
  corr_sims <- list()
  for (kk in 1:K) {
    temp <- numeric(n_MCMC)
    for (tt in 1:n_MCMC) {
      S <- Sigma_sims[tt, kk, , ]
      temp[tt] <- S[1, 2]/prod(sqrt(diag(S)))
    }
    corr_sims[[kk]] <- temp
  }
  return(corr_sims)
}

# ?
get_tau_Omega <- function(Sigma) 
{
  D <- nrow(Sigma)
  tau <- sqrt(diag(Sigma))
  Omega <- matrix(NA, D, D)
  
  for (i in 1:D) {
    for (j in 1:D) {
      Omega[i, j] <- Sigma[i, j]/(tau[i]*tau[j])
    }
  }
  return(list(tau=tau, 
              Omega=Omega))
}

# ?
quad_form_diag <- function(Sigma, tau) {
  diag(tau) %*% Sigma %*% diag(tau)
}


# ? 
mvGauss_emission <- function(obs, mu, Sigma) 
{
  N <- nrow(obs)
  K <- length(mu)
  allprobs <- matrix(NA, N, K)
  for (j in 1:K) {
    allprobs[, j] = dmvnorm(obs, mu[[j]], Sigma[[j]])
  }
  allprobs <- ifelse(!is.na(allprobs), allprobs, 1)
  allprobs
}

# - ?
mvHMM_get_predictive <- function(fit, obs, ndraw = 100)
{
  sims <- rstan::extract(fit)
  niter <- nrow(sims$mu)
  if(ndraw > niter) {stop("ndraw is larger than the number of  mcmc iterations:
                           choose a value between 1 and niter")}
  D <- ncol(obs)
  N <- nrow(obs)
  y_hat <- array(NA, dim = c(ndraw, N, D))
  
  mu_sims <- sims$mu
  gamma_sims <- sims$gamma
  # gamma_init_sims <- sims$gamma_init
  Sigma_sims <- get_Sigma_sims(sims)
  
  idxs <- sample(1:niter, ndraw, replace = FALSE)
  i <- 1
  cat(" ... sampling from posterior predictive... \n")
  for (t in idxs) {
    if((i %% (ndraw/10)) == 0) {
      cat(" ...", as.integer((i/ndraw)*100), "% \n")
    }
    mu_temp <- lapply(1:K, function(k, t) {mu_sims[t, k, ]}, t)
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    gamma_temp <- gamma_sims[t, , ]
    # gamma_init_temp <- gamma_init_sims[t, ]
    z_star <- mvHMM_viterbi(obs, mu_temp, Sigma_temp, 
                            gamma_temp,  draw = TRUE) 
    y_hat[i, , ] <- t(sapply(1:N, function(n) {rmvnorm(1, mu_temp[[z_star[n]]], 
                                                       Sigma_temp[[z_star[n]]])}))
    i <- i + 1
  }
  
  mu_hat <- lapply(1:K, function(kk) {
    apply(mu_sims[, kk, ], 2, mean)})
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  gamma_hat <- apply(gamma_sims, c(2, 3), mean)
  # gamma_init_hat <- apply(gamma_init_sims, 2, mean)
  
  z_hat <- mvHMM_viterbi(obs, mu_hat, Sigma_hat, 
                         gamma_hat, draw = FALSE) 
  
  return(list(y_hat = y_hat, 
              z_hat = z_hat))
}


# generate most likely state sequence using Viterbi algorithm 
# (see Zucchini et al., p.82)
mvHMM_viterbi <- function(obs, mu, Sigma, gamma, draw = FALSE) 
{
  N <- nrow(obs)
  K <- length(mu)
  
  allprobs <- sapply(1:K, function(k, mu, Sigma) {
    dmvnorm(obs, mu[[k]], Sigma[[k]])}, mu = mu, Sigma = Sigma)
  allprobs <- ifelse(!is.na(allprobs), allprobs, 1)
  
  gamma_init <- solve(t(diag(K)-gamma+1),rep(1,K))
  
  xi <- matrix(0, N, K)
  foo <- gamma_init * allprobs[1, ]
  xi[1, ] <- foo/sum(foo)
  
  for (t in 2:N) {
    foo <- apply(xi[t-1, ] * gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo/sum(foo)
  }
  z_star <- numeric(N)
  
  # - Resampling
  if (draw == TRUE) {
    z_star[N] <- sample(1:K, size = 1, prob = xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- sample(1:K, size = 1, prob = gamma[, z_star[t+1]] * xi[t, ])
    }
  } else { # Maximizing
    z_star[N] <- which.max(xi[N, ])
    for (t in (N-1):1) {
      z_star[t] <- which.max(gamma[, z_star[t+1]] * xi[t, ])
    }
  }
  
  z_star
}

# - ? 
mvHMM_plot_predictive_univ <- function(obs, y_hat, z_col, plt_pars, snippet, 
                                       win_s, filter_shifted = TRUE, padded = FALSE) 
{
  plt_pch <- plt_pars$pch
  xlabel <- plt_pars$xlab
  ylabel  <- plt_pars$ylab
  
  N <- length(obs)
  ticks <- seq(from = 1, to = N_all, by = 10)
  x <- round(tail(ticks, 1)/10)/2
  labels <- as.character(seq(from=-x, to = x, 
                             len = length(ticks)))
  line_stimulus <- floor(N_all/2) + 1
  last_NA <- get_last_NA(win_s, filter_shifted, padded)
  
  plot(1:N + last_NA, obs,  xaxt = "n", ylab = ylabel, xlab = xlabel, 
       xlim = c(1, N_all),
       ylim = c(min(obs) - 0.3, max(obs) + 0.3), cex.lab = 1.1)
  
  # credible intervals posterior predictive
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(y_hat[, t], probs = probs))
  polygon(c(1:N, rev(1:N)) + last_NA, c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(1:N, rev(1:N)) + last_NA, c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(1:N, rev(1:N)) + last_NA, c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(1:N, rev(1:N)) + last_NA, c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(1:N + last_NA, cred[5,], col=c_dark, lwd=2)
  points(1:N + last_NA, obs, pch = plt_pch, cex = 1.7, col = z_col)
  
  abline(v = line_stimulus, lwd = 1, lty = "dotted")
  abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
  axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
}

# - ? 
mvHMM_plot_predictive_joint <- function(obs, y_hat, z_hat, 
                                        plt_pars_joint, snippet, win_s, 
                                        filter_shifted = TRUE, padded = FALSE) 
{
  z_col <- plt_pars_joint$zcol 
  N <- nrow(obs)
  D <- ncol(obs)
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  
  
  layout(matrix(c(1, 1, 1, 1, 
                  1, 1, 1, 1, 
                  2, 2, 2, 2, 
                  2, 2, 2, 2), byrow = T, nrow=4, ncol=4))
  par(mai = c(0.6, 0.7, 0.1, 0.3))
  
  
  for (d in 1:D) {
    plt_pars <- list()
    plt_pars$pch <- plt_pars_joint$pch[d]
    plt_pars$xlab <- plt_pars_joint$xlab[d]
    plt_pars$ylab <-  plt_pars_joint$neuromods[d]
    mvHMM_plot_predictive_univ(obs[, d], y_hat[, , d],
                               z_aux, plt_pars, snippet, win_s, 
                               filter_shifted, padded)
  }
}

# - (only bivariate here for now)

mvHMM_get_correlation <- function(sims, obs)
{
  
  mu_sims <- sims$mu
  Sigma_sims <- get_Sigma_sims(sims)
  gamma_sims <- sims$gamma
  
  N <- nrow(obs)
  n_MCMC <- dim(gamma_sims)[1]
  K <- dim(gamma_sims)[2]
  
  corr_seq_sims <- matrix(NA, nrow=n_MCMC, ncol=N)
  
  # maybe make inner function for this routine
  corr_sims <- list()
  for (kk in 1:K) {
    temp <- numeric(n_MCMC)
    for (tt in 1:n_MCMC) {
      S <- Sigma_sims[tt, kk, , ]
      temp[tt] <- S[1, 2]/prod(sqrt(diag(S)))
    }
    corr_sims[[kk]] <- temp
  }
  
  cat(" ... sampling from posterior predictive... \n")
  for (t in 1:n_MCMC) {
    if((t %% (n_MCMC/10)) == 0) {
      cat(" ...", as.integer((t/n_MCMC)*100), "% \n")
    }
    mu_temp <- lapply(1:K, function(k, t) {mu_sims[t, k, ]}, t)
    Sigma_temp <- lapply(1:K, function(k, t) {Sigma_sims[t, k, , ]}, t)
    gamma_temp <- gamma_sims[t, , ]
    z_hat <- mvHMM_viterbi(obs, mu_temp, Sigma_temp,
                           gamma_temp, draw = TRUE)
    for (n in 1:N) {
      k <- z_hat[n]
      corr_seq_sims[t, n] <- corr_sims[[k]][t]
    }
  }
  
  return(corr_seq_sims)
}


# - 
mvHMM_plot_correlation <- function(corr_sims, z_hat, z_col, ylabel, cols, add = F, 
                                   only_mean = F, snippet, substrat, plt_state_probs = T, 
                                   plt_labels = T, win_s, filter_shifted = TRUE, padded = FALSE)
{
  N <- dim(corr_sims)[2]
  z_col <- plt_pars_joint$zcol
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  
  ticks <- seq(from = 1, to = N_all, by = 10)
  x <- round(tail(ticks, 1)/10)/2
  labels <- as.character(seq(from=-x, to = x, 
                             len = length(ticks)))
  line_stimulus <- floor(N_all/2) + 1
  last_NA <- get_last_NA(win_s, filter_shifted, padded)
  
  if (!add) {
    plot(1:N + last_NA, type = "n", xaxt = "n", ylab = ylabel, 
         xlim = c(1, N_all),
         ylim = c(-1, 1), cex.lab = 1.1, xlab = "")
    if (plt_state_probs) {
      for (j in 1:(length(cp_loc) - 1)) {
        rect(cp_loc[j] + last_NA, -1, cp_loc[j+1] +last_NA, 1, 
             col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
      }
    }
  }
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(corr_sims[, t], probs = probs))
  
  if (!only_mean) {
    polygon(c(1:N, rev(1:N)) +last_NA, c(cred[1,], rev(cred[9,])),
            col = scales::alpha(cols[1], 0.2), border = NA)
    polygon(c(1:N, rev(1:N)) +last_NA, c(cred[2,], rev(cred[8,])),
            col = scales::alpha(cols[2], 0.2), border = NA)
    polygon(c(1:N, rev(1:N)) +last_NA, c(cred[3,], rev(cred[7,])),
            col = scales::alpha(cols[3], 0.2), border = NA)
    polygon(c(1:N, rev(1:N)) +last_NA, c(cred[4,], rev(cred[6,])),
            col = scales::alpha(cols[4], 0.2), border = NA)
  }
  lines(1:N +last_NA, cred[5,], col= scales::alpha(cols[5], 0.9), lwd=3)
  
  if (plt_labels == T) {
    abline(v = line_stimulus, lwd = 1, lty = "dotted")
    axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
  }
  abline(v = line_stimulus, lwd = 1, lty = "dotted")
  abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
  
  
  if (add) {
    legend("topleft", legend = groups, 
           lwd = rep(3, length(groups)), 
           col = c(cols_all[[1]][5],cols_all[[2]][5], cols_all[[3]][5]),
           cex = rep(0.9, length(groups)),
           bg = "transparent")
  }
}


# ?
mvHMM_stateprobs <- function(sims, obs) 
{
  K <- dim(sims$mu)[2]
  L_Omega_sims <- sims$L_Omega
  Sigma_sims <- get_Sigma_sims(sims)
  mu_hat <- lapply(1:K, function(kk) {
    apply(sims$mu[, kk, ], 2, mean)})
  Sigma_hat <- lapply(1:K, function(kk) {
    apply(Sigma_sims[, kk, , ], c(2, 3), mean)})
  gamma_hat <- apply(sims$gamma, c(2, 3), mean)
  # gamma_init_hat <- apply(sims$gamma_init, 2, mean)
  
  K <- nrow(gamma_hat)
  N <- nrow(obs)
  
  allprobs <- mvGauss_emission(obs, mu_hat, Sigma_hat)
  fb <- mvHMM_forward_backwards(obs, allprobs, gamma_hat) 
  la <- fb$la
  lb <- fb$lb
  llk <- fb$llk
  
  stateprobs <- matrix(NA, ncol = N, nrow = K)
  for (n in 1:N) stateprobs[, n] <- exp(la[, n] + lb[, n] - llk)
  return(stateprobs)
}


# - ? 
mvHMM_plot_stateprobs <- function(sims, z_col, obs,
                                  plt_labels = T,
                                  snippet, 
                                  most_likely_only = F, z_hat = NULL, win_s, 
                                  filter_shifted = TRUE, padded = FALSE) 
{
  K <- ncol(sims$mu)
  N <- nrow(obs)
  time <- 1:N
  
  ticks <- seq(from = 1, to = N_all, by = 10)
  x <- round(tail(ticks, 1)/10)/2
  labels <- as.character(seq(from=-x, to = x, 
                             len = length(ticks)))
  line_stimulus <- floor(N_all/2) + 1
  last_NA <- get_last_NA(win_s, filter_shifted, padded)
  
  
  state_probs <- mvHMM_stateprobs(sims, obs)
  plot(1:N, type = "n", xaxt = "n", ylab = "", yaxt = 'n',
       cex = 0.5, las = 1,
       xlim = c(1, N_all),
       ylim = c(0,1), cex.lab = 1.1, xlab = "")
  
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  cp_loc <- c(1, which(diff(as.numeric(factor(z_aux))) != 0), N)
  
  if (most_likely_only) { 
    for (j in 1:(length(cp_loc) - 1)) { 
      rect(cp_loc[j] + last_NA, -3, cp_loc[j+1] + last_NA, 3, 
           col =scales::alpha(z_aux[cp_loc[j+1]], 0.8), border = F)
    }
  } else {
    
    plot_p <- matrix(NA,nrow=N-1,ncol=2*K)
    a <- 0
    for (n in 2:N) {
      
      for (j in 1:K) {
        plot_p[n-1, (j*2-1)] <- state_probs[j, n-1]
        plot_p[n-1, j*2] <- state_probs[j, n]
        if	(j==1){col_states<- z_col[1]}					
        if	(j==2){col_states<- z_col[2]}					
        if	(j==3){col_states<- z_col[3]}	
        if	(j==4){col_states<- z_col[4]}	
        
        # (at some point need to make function to 
        # generalize for all K, like this is redundant)
        if	(j==1){				
          point_1<-a
          point_2<-point_1+plot_p[n-1,(j*2-1)]
          point_4<-a
          point_3<-point_4+plot_p[n-1,(j*2)]	}
        
        if	(j==2){				
          point_1<-a+plot_p[n-1,(j-1)*2-1]
          point_2<-point_1+plot_p[n-1,(j*2-1)]
          point_4<-a+plot_p[n-1,(j-1)*2]
          point_3<-point_4+plot_p[n-1,(j*2)]	}
        
        if	(j==3){				
          point_1<-a+plot_p[n-1,(j-2)*2-1]+plot_p[n-1,(j-1)*2-1]
          point_2<-point_1+plot_p[n-1,(j*2-1)]
          point_4<-a+plot_p[n-1,(j-2)*2]+plot_p[n-1,(j-1)*2]
          point_3<-point_4+plot_p[n-1,(j*2)]}
        if (j==4) {
          point_1 <- a+ plot_p[n-1,(j-3)*2-1] + plot_p[n-1,(j-2)*2-1]+plot_p[n-1,(j-1)*2-1]
          point_2 <- point_1+plot_p[n-1,(j*2-1)]
          point_4 <- a + plot_p[n-1,(j-3)*2] + plot_p[n-1,(j-2)*2]+plot_p[n-1,(j-1)*2]
          point_3 <- point_4+plot_p[n-1,(j*2)]}
        
        polygon(c(time[n-1],time[n-1],time[n],time[n]) + last_NA,
                c(point_1,point_2,point_3,point_4),col=scales::alpha(col_states, 0.8), border=NA)
        lines(c(time[n-1],time[n]) + last_NA,c(point_2,point_3), col=scales::alpha(col_states, 0.5))
      }
    }
    
    abline(v = line_stimulus, lwd = 1, lty = "dotted")
    axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
    
    if (plt_labels) {
      title(ylab='Prob State',)
      abline(v = line_stimulus, lwd = 1, lty = "dotted")
      axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
    } else {
      title(ylab='')
    }
  }
}



# mvHMM_plot_mean 
mvHMM_plot_mean <- function(obs, y_hat, z_col, plt_pars, add = F, plt_states = T, 
                            col_lines = "black", plt_obs = F, snippet, plt_means = T, 
                            win_s, filter_shifted = TRUE, padded = FALSE) 
{
  plt_pch <- plt_pars$pch
  plt_lty <- plt_pars$lty
  xlabel <- plt_pars$xlab
  ylabel  <- plt_pars$ylab
  ylimit <- plt_pars$ylim
  cex <- plt_pars$cex
  
  # - 
  N <- length(obs)
  last_NA <- get_last_NA(win_s, filter_shifted, padded)
  
  
  # ticks <- seq(from = 1, to = N, by = 10)
  # x <- round(tail(ticks, 1)/10)/2
  # labels <- as.character(seq(from=-x, to = x, 
  #                            len = length(ticks)))
  # line_stimulus <- round(N/2) + 1
  
  # credible intervals posterior predictive
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(y_hat[, t], probs = probs))
  cp_loc <- c(1, which(diff(as.numeric(factor(z_col))) != 0), N)
  
  if (!add) {
    if (plt_means) {
      plot(1:N, type = "n", xaxt = "n", ylab = ylabel, xlab = xlabel, 
           cex = 1.4, las = 1,
           xlim = c(1, N_all),
           ylim = ylimit, cex.lab = 1.5)
      lines(1:N + last_NA, cred[5, ], lty = plt_lty, lwd = 3)
      # plot(1:N + last_NA, cred[5, ], type = "l", xaxt = "n", 
      #      ylab = ylabel, xlab = xlabel, col = col_lines, lwd = 2,
      #      cex.lab = 1.5, cex = 1.4, ylim = ylimit, lty = plt_lty)
    } else {
      plot(1:N, type = "n", xaxt = "n", ylab = ylabel, 
           xlim = c(1, N_all),
           ylim = ylimit, cex.lab = 1.1, xlab = "")
      # plot.new(ylim = ylimit, cex.lab = 1.5, cex = 1.4, 
      #          ylab = ylabel, xlab = xlabel)
    }
    if (plt_states) {
      for (j in 1:(length(cp_loc) - 1)) {
        rect(cp_loc[j] + last_NA, -3, cp_loc[j+1] + last_NA, 3, 
             col =scales::alpha(z_col[cp_loc[j+1]], 0.8), border = F)
      }
    }
  }
  
  if(plt_means) {
    lines(1:N + last_NA, cred[5,], col=col_lines, lwd=3, lty = plt_lty)
  }
  
  if (plt_obs) {
    points(1:N  + last_NA, obs, col=scales::alpha("white", .9),
           pch=plt_pch, cex= (cex + 0.2))
    points(1:N + last_NA, obs, pch = plt_pch,
           cex = cex, col =  scales::alpha("black", .7))
  }
  
  # abline(v = line_stimulus, lwd = 1, lty = "dotted")
  # abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
  # axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
  # 
}



# - ? 
bivariateHMM_plot_mean <- function(obs, y_hat, z_hat, z_col, neuromods, 
                                   snippet,
                                   col_lines = "black", 
                                   plt_obs = F, plt_states = T, 
                                   add = F, plt_means = T, win_s, 
                                   filter_shifted = TRUE, padded = FALSE) 
{
  N <- nrow(obs)
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  plt_pars <- list()
  plt_pars$pch <- 16
  plt_pars$lty <- "solid"
  plt_pars$cex <- 1.8
  plt_pars$xlab <- ""
  plt_pars$ylab <-  ""
  plt_pars$ylim <- c(-3, 3)
  
  # mvHMM_plot_mean_univ(obs[, 1], y_hat[, , 1], z_aux, plt_pars, 
  #                      NA_smoothing)
  mvHMM_plot_mean(obs[, 1], y_hat[, , 1], z_aux, plt_pars, 
                  plt_obs = plt_obs, 
                  plt_states = plt_states,
                  col_lines = col_lines,
                  plt_means = plt_means, 
                  snippet = snippet, 
                  add = add, win_s, 
                  filter_shifted = filter_shifted,
                  padded = padded)
  # here probably need to change add = T, and the colors of the line, 
  
  plt_pars$pch <- 18
  plt_pars$cex <- 2.0
  plt_pars$lty <- "dotdash"
  
  
  # mvHMM_plot_mean_univ(obs[, 2], y_hat[, , 2], z_aux, plt_pars, NA_smoothing, 
  #                      add = T)
  mvHMM_plot_mean(obs[, 2], y_hat[, , 2], z_aux, plt_pars, 
                  plt_obs = plt_obs, 
                  plt_states = plt_states,
                  col_lines = col_lines,
                  plt_means = plt_means,
                  snippet = snippet, 
                  add = T, win_s, 
                  filter_shifted = filter_shifted,
                  padded = padded)
  
  if (!add) {
    
    ticks <- seq(from = 1, to = N_all, by = 10)
    x <- round(tail(ticks, 1)/10)/2
    labels <- as.character(seq(from=-x, to = x,
                               len = length(ticks)))
    line_stimulus <- floor(N_all/2) + 1
    last_NA <- get_last_NA(win_s, filter_shifted, padded)
    
    
    title(ylab=str_c(neuromods, 
                     collapse = ", "), cex.lab = 1.2, 
          xlim = c(1, N_all))
    abline(v = line_stimulus, lwd = 1, lty = "dotted")
    abline(v = line_stimulus + 10, lwd = 1, lty = "dotted")
    axis(side = 1, ticks, labels, cex.axis = 1.2, tck=-0.04)
    
    if (plt_means) {
      legend("topright", 
             legend = c(paste("Mean(", neuromods[1], ")", sep =""), 
                        paste("Mean(", neuromods[2], ")", sep ="")),
             lty = c("solid", "dotdash"),
             lwd = c(1.5, 1.5),
             cex = c(0.7, 0.7),
             bg = "transparent",
             col = c("black", "black"))
    } else {
      legend("topright", 
             legend = c(neuromods[1], neuromods[2]),
             pch = c(16, 18),
             cex = c(0.9, 0.9),
             bg = "transparent",
             col = c(scales::alpha("black", 0.5),
                     scales::alpha("black", 0.5))
      )
    }
  }
}

# - ? 
mvHMM_plot_joint <- function(obs, z_hat, plt_pars) 
{
  N <- nrow(obs)
  labels <- plt_pars$neuromods
  z_col <- plt_pars$zcol
  z_aux <- sapply(1:N, function(n) z_col[z_hat[n]])
  plot(obs[, 1], obs[, 2], col = z_aux, type = "p", pch = 19, 
       cex = 1.9, xlab = labels[1], ylab=labels[2], cex.lab = 1.2)
  for (t in 1:(N-1)) {
    arrows(obs[t, 1], obs[t, 2], obs[t+1, 1], obs[t+1, 2], 
           length = 0.05, col = "black") 
  }
}

# - 
mvHMM_select_best_model <- function(ID, snippets, substrats, 
                                    groups, folder_analysis, allTrials = F, 
                                    Ks = c(3, 4))
{
  
  
  if (!allTrials) {
    n_models <- length(unlist(substrats_groups)) * length(snippets)
    df_summary <- tibble(Snippet=character(), 
                         Substrat=character(),
                         Group=character(), 
                         K_best=integer(), 
                         Warnings = logical())
    bayes_summary <- list(n_models) # there four different time series eventually
    jj = 1 # auxiliary index
    
    for (snippet in snippets) {
      for (substrat in substrats) {
        
        cat("SUBSTRAT: ", substrat, " ---------------------------------------", "\n")
        
        
        groups = substrats_groups[[substrat]]
        
        for (group in groups) {
          #cat("GROUP: ", group, "%%%%%%%%%%%%%%%%%%%%%-", "\n")
          
          AICs <- rep(Inf, length(Ks) + 1)
          BICs <- rep(Inf, length(Ks) + 1)
          LogMLs <- rep(Inf, length(Ks) + 1)
          warnings <- rep(NA, length(Ks) +1)
          
          for (K in Ks) {
            path_main <- paste("ID", ID, "/", snippet, "_by",  
                               substrat, "_", group,  "_K", K, 
                               "_performance.RData", sep ="")
            path <- paste(getwd(), folder_analysis, path_main, sep="/")
            load(path)
            
            estimates <- perf$bayes_est
            AICs[K] <- perf$AIC
            BICs[K] <- perf$BIC
            LogMLs[K] <- perf$logml
            
            warnings[K] <- perf$issues
          }
          
          # ---- selecting best model (using AIC) --
          #K_best <- which(AICs == min(AICs))
          #K_best <- which(BICs == min(BICs))
          K_best <- which(LogMLs == min(LogMLs))
          #K_best <- 4
          df_summary <- df_summary %>% add_row(tibble_row(Snippet = snippet,
                                                          Substrat = substrat,
                                                          Group = group, 
                                                          K_best = K_best, 
                                                          Warnings = warnings[K_best]))
          
          path_main <- paste("ID", ID, "/", snippet, "_by",  
                             substrat, "_", group,  "_K", K_best, 
                             "_performance.RData", sep ="")
          path <- paste(getwd(), folder_analysis, path_main, sep="/")
          load(path)
          estimates <- perf$bayes_est
          bayes_summary[[jj]] <- list(mu = estimates$mu, 
                                      corr = estimates$corr)
          jj = jj + 1
        }
      }
    }
    return(list(df = df_summary, est = bayes_summary))
    
  } else {
    n_models <- length(snippets)
    
    df_summary <- tibble(Snippet=character(), 
                         K_best=integer(), 
                         Warnings = logical())
    bayes_summary <- list(n_models) # there four different time series eventually
    jj = 1
    
    for(snippet in snippets) {
      
      AICs <- rep(Inf, length(Ks) + 1) 
      BICs <- rep(Inf, length(Ks) + 1)
      LogMLs <- rep(Inf, length(Ks) + 1)
      warnings<- rep(NA, length(Ks) + 1)
      
      for (K in Ks) {
        path_main <- paste("ID", ID, "/", snippet, "_allTrials_K", K, 
                           "_performance.RData", sep ="")
        path <- paste(getwd(), folder_analysis, path_main, sep="/")
        load(path)
        AICs[K] <- perf$AIC
        BICs[K] <- perf$BIC
        LogMLs[K] <- perf$logml
        warnings[K] <- perf$issues
      }
      
      # ---- selecting best model (using AIC) --
      #K_best <- which(AICs == min(AICs))
      #K_best <- which(BICs == min(BICs))
      K_best <- which(LogMLs == min(LogMLs))
      # K_best <- 4
      
      df_summary <- df_summary %>% add_row(tibble_row(Snippet = snippet,
                                                      K_best = K_best, 
                                                      Warnings = warnings[K_best]))
      path_main <- paste("ID", ID, "/", snippet, "_allTrials_K", K_best,
                         "_performance.RData", sep ="")
      path <- paste(getwd(), folder_analysis,path_main, sep="/")
      load(path)
      estimates <- perf$bayes_est
      bayes_summary[[jj]] <- list(mu = estimates$mu, 
                                  corr = estimates$corr)
      jj = jj + 1
    }
    return(list(df = df_summary, est = bayes_summary))
  }
}


