library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 3000
set.seed(42)
res <- simuQFD(N = N, incrm_firing = 8, 
                theta_quiet = c(0.99,0.01),
                theta_firing = c(0.3, 0.7),
                theta_decay = c(0.1,0,0.9),rate_decay = 0.8, 
                sigma_decay = 0.5, 
                sigma_firing = 2, sigma_quiet = 2)

tt <- seq(0,50,length.out = N)
f <- 5 * sin(6*tt/(2*pi))
plot(res$timeseries)
simu_signal <- f+res$timeseries
plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")

QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                sigma_prior = c(-2,2),
                period_prior = c(-1,3),
                Q0_prior = c(-2,2),
                dQ_prior = c(-2,2),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,0.1,1),
                mu0_quiet = 0,
                lambda_quiet = 1e-4,
                gamma_quiet = c(0.01,0.01),
                #alpha_incre_firing = 0.15,
                mu0_incre_firing = 0,
                sigma_incre_firing = 1e3,
                gamma_firing = c(0.01,0.01),
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                gamma_decay = c(0.01,0.01),
                diag = rep(0,N)
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFD.stan', 
            model_name = "celeritQFD", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))
vbQFD <- vb(modelQFD, data = QFD_data,tol_rel_obj = 0.0005)

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)
