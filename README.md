## Simultaneous Detrending and Flare Detection using Gaussian Process with Hidden Markov Model

A working repo using Gaussian Process as trend model with a Hidden Markov flare detection.


In current stage, several methods for simultaneous detrending and flare detection with photometric data are implemented in [Stan](https://mc-stan.org/), a state-of-the-art platform for statistical modeling. 

- [Vanilla HMM](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/Vanilla): this is the usual HMM model with two normal states, which is not very effective
- [Gaussian Process](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/GP): this is the usual GP model with 0 mean, can be used to detrend the data
- [Celerite](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/Celerite): adopted from Python package [exoplanet-celerite2](https://github.com/exoplanet-dev/celerite2), one could use this to detrend the data and use 3-sigma rule
- [CeleriteHMM](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Prototypes/CeleriteHMM): this type of models combined Celerite2 and HMMs. We used the Celerite as trend model. For flare, we have two states namely a quiet state and a flare state while quiet state the data follow normal with mean of the Celerite trend and flare state follow the trend added by some random variable from different distribution, e.g. Laplace and power-law.

And several of them are experimental 

- [CeleriteQFD](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Morphology/QFD), this model relies more on the morphology, similar to CeleriteHMM but we have three states, anmely **Q**uiet-**F**iring-**D**ecay, and for **F** and **D** state, we model the data generating process as a AR(1) model without slope (for **F**) or intercept (for **D**), this impose a linear increase and an exponential derease of the brightness. 
- [QFDexN](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Morphology/QFD), similar to QFD, we have Firing process follow a exponentially modified normal random walk of the previous step
- [CeleriteQFDexN](https://github.com/YunyiShen/AstroHMMs/tree/master/Stan/Morphology/QFD), the celerite combined version of QFDexN, currently not working very well as of June 2021.