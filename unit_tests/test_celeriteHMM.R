library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = T)

set.seed(42)
t <- seq(0,50,.5)
f <- sin(4*t/(2*pi))
#plot(t,f)
y <- 5*f + rnorm(length(t),0,1)
#y <- rnorm(length(t),0,1)
N <- length(t)
n_fire <- 10
firing <- sample( N, n_fire)
y[firing] <- f[firing] + rexp(n_fire, 0.05)
plot(t,y)
points(t[firing], y[firing], col = "red")

firing_or_not <- y/y
firing_or_not[firing] <- 2

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
                     eps_neg = 1e-2)

modellaplace <- stan_model(file = './Stan/Prototypes/CeleriteHMM/celerite_HMM_laplace.stan', 
            model_name = "celeritHMMlaplace", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitlaplace <- sampling(modellaplace, data = star_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)
#fit2laplace <- optimizing(modellaplace, data = star_data, verbose = T, iter = 10000)
summ_fitlaplace <- summary(fitlaplace)
plot(summ_fitlaplace[[1]][1:101+217,1])
plot(summ_fitlaplace[[1]][1:101+117,1], ylim = c(1,2))
state_estlaplace <- summ_fitlaplace[[1]][1:101+117,1]
lines(firing_or_not)

# plot
par(mfrow = c(1,2))
plot(t,y, main = "ground truth")
points(t[firing], y[firing], col = "red")
plot(t,y, main = "celeriteHMM-Laplace firing")
points(t[state_estlaplace>=1.5], y[state_estlaplace>=1.5], col = "red")


modelpowerlaw <- stan_model(file = './Stan/Prototypes/CeleriteHMM/celerite_HMM_powerlaw.stan', 
            model_name = "celeritHMMpowerlaw", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))
fitpowerlaw <- sampling(modelpowerlaw, data = star_data, iter = 2000,control = list(max_treedepth=15))



modelexp <- stan_model(file = './Stan/Prototypes/CeleriteHMM/celerite_HMM_exp.stan', 
            model_name = "celeritHMMexp", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitexp <- sampling(modelexp, data = star_data, iter = 2000,control = list(max_treedepth=15))
summ_fitexp <- summary(fitexp)
plot(summ_fitexp[[1]][1:101 + 319,1])
plot(summ_fitexp[[1]][1:101 + 117,1], ylim = c(1,2))
lines(firing_or_not)


modelgumble <- stan_model(file = './Stan/Prototypes/CeleriteHMM/celerite_HMM_gumble.stan', 
            model_name = "celeritHMMgumble", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitgumble<- sampling(modelgumble, data = star_data,control = list(adapt_delta = 0.99))
summ_fitgumble <- summary(fitgumble)
plot(summ_fitgumble[[1]][1:101,1])
plot(summ_fitgumble[[1]][1:101+117,1], ylim = c(1,2))
lines(firing_or_not)
