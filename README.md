# mvHMM

mvHMM contains the stan software (as well as R utilities) to model pupil and norepinephrine (NA) estimates using a Bayesian hidden Markov model (HMM) as detailed in Bang et. al (2023), "Noradrenaline tracks emotional modulation of attention in human amygdala", published in Current Biology. 

## Overview

We investigate the relationship between the pupil and NA estimates, where we fitted an HMM. The HMM modelled the joint pupil and NA estimates as a sequence of visits to a finite number
of hidden states, with each state represented as a bivariate Gaussian distribution with state-specific means and a state-specific variance-covariance matrix for the pupil and NA estimates.

This software is illustrated in `tutorial.Rmd`, where we show results from three scenarios: (i) all trials irrespective of stimulus type; (ii) trials with oddball and high arousal; (iii) trials with oddball and low arousal. The number of latent states, K, for an HMM was selected by inspecting the posterior predictive fits and, in a more principled way, by calculating the ratio of marginal likelihoods from different models using the R package `bridgesampling`, whose compatibility with Stan makes it straightforward to estimate the marginal likelihood directly from a Stan output. 

* Main Inference Script:
```
tutorial.Rmd
```

## Example - mvHMM

Here is an example of using our (stan) software in R from `tutorial.Rmd`.
