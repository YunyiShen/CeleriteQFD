library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

set.seed(42)
t <- seq(0,50,.1)
f <- sin(5*t/(2*pi))
#plot(t,f)
y <- 5*f + rnorm(length(t),0,1)
#y <- rnorm(length(t),0,1)
N <- length(t)
n_fire <- 50
firing <- sample( N, n_fire)
y[firing] <- f[firing] + rexp(n_fire, 0.05)
plot(t,y)

star_data <- list(N=N, t = t, y = y,
                     sigma_prior = c(-2,2),
                     period_prior = c(-1,3),
                     Q0_prior = c(-2,2),
                     dQ_prior = c(-2,2),
                     f_prior = c(1e-6,1-1e-6),
                     alpha = c(1,1), # transition
                     mu0 = 0, # usual mean
                     lambda = 0.001, # prior for usual
                     noise_prior = c(0.01,0.01),
                     rate = 0.0001, 
                     diag = 0*t,
                     eps_pare = 1)

modellaplace <- stan_model(file = './Stan/celerite_HMM_laplace.stan', 
            model_name = "celeritHMM", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitlaplace <- sampling(modellaplace, data = star_data)
fit2laplace <- optimizing(modellaplace, data = star_data, verbose = T, iter = 10000)
summ_fitlaplace <- summary(fitlaplace)
plot(summ_fit[[1]][1:501,1])
plot(fit2laplace$par[1:501])


modelpowerlaw <- stan_model(file = './Stan/celerite_HMM_powerlaw.stan', 
            model_name = "celeritHMM", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fit2powerlaw <- optimizing(modelpowerlaw, data = star_data, verbose = T, iter = 10000)
plot(y-fit2powerlaw$par[1:501]) 
plot(fit2powerlaw$par[1023:1519])
