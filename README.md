# mvHMM

mvHMM contains the stan software (as well as R utilities) to model pupil and norepinephrine (NA) measurements using a Bayesian hidden Markov model (HMM) as detailed in Bang et. al (2023), "Noradrenaline tracks emotional modulation of attention in human amygdala", published in Current Biology. 

## Overview

We investigate the relationship between the pupil and NA measurements. The HMM is modeling the joint pupil and NA measurements as a sequence of visits to a finite number
of hidden states, with each state represented as a bivariate Gaussian distribution with state-specific means and a state-specific variance-covariance matrix for the pupil and NA measurements. 
The number of latent states, K, for an HMM was selected by inspecting the posterior predictive fits and, in a more principled way, by calculating the ratio of marginal likelihoods from different models using the R package `bridgesampling`, whose compatibility with Stan makes it straightforward to estimate the marginal likelihood directly from a Stan output. 


This software is illustrated in `tutorial.Rmd`, where we show results from three different scenarios: 

1) Trials with oddball and high arousal
2) Trials with oddball and low arousal.
3) All trials irrespective of stimulus type. 


## Example - mvHMM for analysis of pupil and norepinephrine

Here is an example of using our (stan) software in R from `tutorial.Rmd`

```r
get_neuromodulators(ID = 123, info_vars, snippet = "stimulus",
                    substrat =  "Oddball_x_Arousal",
                    groups = c("Oddball_ArousLow","Oddball_ArousHigh"), 
                    win_s = 5,
                    detrend = T, plt = T)

```

<p align="center">
<img src="https://github.com/Beniamino92/mvHMM/blob/master/figures/data_substrat.png" width="500" heigth="250"/> 
</p>
