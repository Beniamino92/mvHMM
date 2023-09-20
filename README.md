# mvHMM
Bayesian Multivariate HMM for modeling neuromodulators

mvHMM contains the stan software (as well as R utilities) to model time series and sequential data using a Bayesian HMM that is a reformulation of any given HSMM as detailed in "Bayesian Approximations to Hidden Semi-Markov Models for Telemetric Monitoring of Physical Activity" (2022) by B.Hadj-Amar, J.Jewson, M.Fiecas
We implemented a bivariate HMM based on a discrete latent sequence that partitions the pupil and norepineprhine estimates into different, potentially recurring regimes.
