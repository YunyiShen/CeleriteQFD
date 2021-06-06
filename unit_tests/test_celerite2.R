library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = F)

t <- seq(0,50,.1)
f <- sin(5*t/(2*pi))
#plot(t,f)
y <- 10*f + rnorm(length(t),0,5)
plot(t,y)
N <- length(t)

celerite_data <- list(N=N, t = t, y = y,
                     sigma_prior = c(-2,2),
                     period_prior = c(-1,3),
                     Q0_prior = c(-2,2),
                     dQ_prior = c(-2,2),
                     f_prior = c(1e-6,1-1e-6),
                     err_prior = c(1,0.01),
                     diag = 0*t)

modelcelerite <- stan_model(file = './Stan/celerite.stan', 
            model_name = "celerit", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fit <- sampling(modelcelerite, data = celerite_data)
fit2 <- optimizing(modelcelerite, data = celerite_data, verbose = T, iter = 10000)
par_opt <- fit2$par
summ_fit <- summary(fit)
plot(summ_fit[[1]][1:501,1])
plot(y-summ_fit[[1]][1:501,1])
plot(fit2$par[1:501])
plot(y-fit2$par[1:501])
plot(fit2$par[1023:1519])
