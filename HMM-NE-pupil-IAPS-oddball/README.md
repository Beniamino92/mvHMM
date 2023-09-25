# HMM-NE-pupil-IAPS-oddball

`/HMM-NE-pupil-IAPS-oddball/` contains the stan software (as well as R utilities) to model pupil and norepinephrine (NE) measurements using a Bayesian hidden Markov model (HMM) as detailed in Bang et. al (2023), "Noradrenaline tracks emotional modulation of attention in human amygdala", published in Current Biology. 

## Overview

We investigate the relationship between the pupil and NE measurements. The HMM is modeling the joint pupil and NE measurements as a sequence of visits to a finite number
of hidden states, with each state represented as a bivariate Gaussian distribution with state-specific means and a state-specific variance-covariance matrix for the pupil and NA measurements. 
The number of latent states, K, for an HMM was selected by inspecting the posterior predictive fits and, in a more principled way, by calculating the ratio of marginal likelihoods from different models using the R package `bridgesampling`, whose compatibility with Stan makes it straightforward to estimate the marginal likelihood directly from a Stan output. 


This software is illustrated in `tutorial_HMM-NE-pupil.Rmd`, where we show results from three different scenarios: 

1) Trials with oddball and high arousal
2) Trials with oddball and low arousal.
3) All trials irrespective of stimulus type. 


## Example - mvHMM for analysis of pupil and norepinephrine

We provide a snapshot of  `tutorial_HMM-NE-pupil.Rmd` for using our (stan) software in R 


* Get and plot time series of neuromodulator and pupil, for selected substratification 
```r
get_neuromodulators(ID = 123, info_vars, snippet = "stimulus",
                    substrat =  "Oddball_x_Arousal",
                    groups = c("Oddball_ArousLow","Oddball_ArousHigh"), 
                    win_s = 5,
                    detrend = T, plt = T)

```
<p align="center">
<img src="https://github.com/Beniamino92/mvHMM/blob/main/figures/data_substrat.png" width="400" heigth="170"/> 
</p>

* Run MCMC sampler
```r
K = 3 # number of states (selected among competing models using bridgesampling)
fit <- stan(file = stan_path,
            data = list(N = nrow(obs_group), D = 2, 
                                  K = K, y = obs_group, 
                                  mu_loc = 0,
                                  mu_scale = mu_scale, 
                                  alpha_0 = matrix(c(rep(1, K*K)), 
                                                   nrow = K, ncol = K, byrow = TRUE),
                                  tau_loc = 1, tau_scale = tau_scale,
                                  Omega_shape = 1),  
            seed = 123, 
            chains = 1, iter = n_MCMC, cores = 1)
```

* Get posterior predictive + plot
```r
mvHMM_predictive <- mvHMM_get_predictive(fit, obs_group, ndraw = 200)
z_hat <- mvHMM_predictive$z_hat
y_hat <- mvHMM_predictive$y_hat
mvHMM_plot_predictive_joint(obs_group, y_hat,
                            z_hat, plt_pars_joint,
                            snippet = snippet, win_s = win_s)
```
<p align="center">
<img src="https://github.com/Beniamino92/mvHMM/blob/main/figures/OddballLowArousal_postpred-1.png" width="400" heigth="170"/> 
</p>

* Get time-varying posterior correlation + plot
```r
corr_sims <- mvHMM_get_correlation(sims, obs_group)
mvHMM_plot_correlation(corr_sims, z_hat, plt_pars_joint$zcol,
                       corr_label, cols_grey, snippet = "stimulus",
                       substrat = "allTrials",
                       plt_state_probs = T, win_s = win_s)
```
<p align="center">
<img src="https://github.com/Beniamino92/mvHMM/blob/main/figures/OddballLowArousal_correlation-1.png" width="400" heigth="170"/> 
</p>

* Get time-varying state probabilities + plot
```r
mvHMM_plot_stateprobs(sims, plt_pars_joint$zcol,
                      obs_group, snippet = "stimulus", win_s = 5)
```
<p align="center">
<img src="https://github.com/Beniamino92/mvHMM/blob/main/figures/OddballLowArousal_stateprob-1.png" width="400" heigth="170"/> 
</p>

