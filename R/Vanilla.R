library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

rawdata <- read.csv(list.files("./Data", full.names = TRUE)[1])

timeseries <- rawdata$PDCSAP_FLUX[-1] # first entry is NA...
observed <- 1 * (!is.na(timeseries))

timeseries[is.na(timeseries)] <- timeseries[1]

T <- length(timeseries)

star_stan <- list(T=T,
                    K=2, 
                    y = timeseries,
                    observed = observed, 
                    alpha = c(1/2,1/2),
                    mu0=0,
                    lambda = 1e-10,
                    shape = 0.001,
                    rate = 0.001)

fit <- stan(file = './Stan/Vanilla.stan', data = star_stan, 
                iter = 3000, control = list(adapt_delta = .8))
plot(fit)