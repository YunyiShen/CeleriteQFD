library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 300
set.seed(42)
res <- simuQFDexN(N = N, rate_firing = .2, 
                theta_quiet = c(0.95,0.05),
                theta_firing = c(0.4, 0.6),
                theta_decay = c(0.15,0,0.85),rate_decay = 0.7, 
                sigma_noise = 1)

tt <- seq(0,(50/300)*N,length.out = N)
f <- 5 * sin(3*tt/(2*pi))
plot(res$timeseries)
simu_signal <- f+res$timeseries
plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")



QFDsimple_data <- list(N=N,
                y = res$timeseries,
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = 1e-4,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3
                )


modelQFDsimple <- stan_model(file = './Stan/Morphology/QFD/QFDexN.stan', 
            model_name = "QFTexN")

optQFDsimple <- optimizing(modelQFDsimple, QFDsimple_data)
vbQFDsimple <- vb(modelQFDsimple, data = QFDsimple_data,tol_rel_obj = 0.0005)
fitQFDsimple <- sampling(modelQFDsimple, data = QFDsimple_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)



QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                sigma_prior = c(-3,3),
                period_prior = c(-2,7),
                Q0_prior = c(-3,3),
                dQ_prior = c(-3,3),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = 1e-4,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(0,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite.stan', 
            model_name = "celerit", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

vbcelerite <- vb(modelcelerite, data = QFD_data, tol_rel_obj = 0.001)
fitcelerite <- sampling(modelcelerite, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)
summ_celerite <- summary(fitcelerite)
plot(summ_celerite[[1]][1:300 + 310,1])

vbQFD <- vb(modelQFD, data = QFD_data,tol_rel_obj = 0.001)
summ_vbQFD <- summary(vbQFD)
plot(summ_vbQFD[[1]][1:300 + 322,1])
fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)

fitQFD_array <- as.array(fitQFD)
hist(fitQFD_array[,1,313])



