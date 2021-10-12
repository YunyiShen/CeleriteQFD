## Simultaneous Detrending and Flare Detection using Gaussian Process with Hidden Markov Model

A working repo using Gaussian Process as trend model with a Hidden Markov flare detection.

### Models implemented
In current stage, several methods for simultaneous detrending and flare detection with photometric data are implemented in [Stan](https://mc-stan.org/), a state-of-the-art platform for statistical modeling. 

- [Vanilla HMM](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/Vanilla): this is the usual HMM model with two normal states, which is not very effective
- [Gaussian Process](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/GP): this is the usual GP model with 0 mean, can be used to detrend the data
- [Celerite](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/Celerite/celerite.stan): adopted from Python package [exoplanet-celerite2](https://github.com/exoplanet-dev/celerite2), one could use this to detrend the data and use 3-sigma rule. Check [the interface](https://github.com/YunyiShen/AstroHMMs/tree/master/celerite2/celerite2.hpp) out if you are interested in using celerite2 in Stan environment and thus need an interface between C++ and stan. Cholesky decomposition of several recommended kernels were implemented that ready to use. 



- [CeleriteHMM](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/CeleriteHMM): this type of models combined Celerite2 and HMMs. We used the Celerite as trend model. For flare, we have two states namely a quiet state and a flare state while quiet state the data follow normal with mean of the Celerite trend and flare state follow the trend added by some random variable from different distribution, e.g. Laplace and power-law.

- [QFDexN](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Morphology/QFD/QFDexN.stan), this model relies more on the morphology, similar to CeleriteHMM but we have three states, namely **Q**uiet-**F**iring-**D**ecay, and for **F** and **D** state, we model the data generating process as a AR(1) model with an exponential modified normal (exN) random walk for **F**iring and an exponential decay for **D**ecay, this works well as pure firing model.
- [CeleriteQFDexN](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Morphology/QFD/), the celerite combined version of QFDexN, With different type of celerite kernels implemented, currently not very stable for time series longer than 2,000 probably due to stan. 

### Injection recovery

Some amount of injection-recovery simulation was done with TIC-131799991's time series with injected (rather small) Kepler flares. 

One example run, compared with Celerite and 1-3-sigma rule

![](https://github.com/YunyiShen/AstroHMMs/raw/master/Res/Injection_recover/small_flares/QFD_example.jpg)

Result for 100 injection recovery simulations, evaluated by sensitivity and specificity:

![](https://github.com/YunyiShen/AstroHMMs/raw/master/Res/Injection_recover/small_flares/injection-recoveryz-smallflare.jpg)


### Real data

An example run from TIC 131799991

![](https://github.com/YunyiShen/AstroHMMs/raw/master/Res/CeleriteQFD/131799991_16400-17400/det.png)

The code to repeat this result is in `./tests/Short_time_series_tests/celeriteQFDexN.R`, it requires `rstan` (currently does not work for `cmdstan` due to the order it include C++ header files).


