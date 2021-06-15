library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 300
set.seed(12345)
res <- simuQFDexN(N = N, rate_firing = .1, 
                theta_quiet = c(0.95,0.05),
                theta_firing = c(0.2, 0.8),
                theta_decay = c(0.2,0,0.8),rate_decay = 0.6, 
                sigma_noise = 1)

tt <- seq(0,(50/300)*N,length.out = N)
f <- 5 * sin(2*tt/(2*pi))
plot(tt,f)
plot(res$timeseries)
simu_signal <- f+res$timeseries
plot(simu_signal)
points(which(res$state==2),simu_signal[res$state==2], col = "red")
points(which(res$state==3),simu_signal[res$state==3], col = "blue")

QFD_data_simple <- list(N=N, t = tt,
                y = res$timeseries,
                sigma_prior = c(-5,5),
                period_prior = c(-2,2),
                Q0_prior = c(-5,5),
                dQ_prior = c(-5,5),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = .01,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = -1,
                sigma_rate_decay = 1e3,
                diag = rep(0,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                sigma_prior = c(-5,5),
                period_prior = c(-2,2),
                Q0_prior = c(-5,5),
                dQ_prior = c(-5,5),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,1),
                alpha_firing = c(1,1),
                alpha_decay = c(1,1,1),
                mu0_quiet = 0,
                lambda_quiet = .01,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(0,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

modelQFDsimple <- stan_model(file = './Stan/Morphology/QFD/QFDexN2.stan', 
            model_name = "QFTexN2")

optQFDsimple <- optimizing(modelQFDsimple, QFD_data_simple)
temp <- res$timeseries[-1]
viterbi <- optQFDsimple$par[1:299 + 13]
plot(temp)
points(which(viterbi==2),temp[viterbi==2], col = "red")
points(which(viterbi==3),temp[viterbi==3], col = "blue")



modelQFDexN2 <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN2.stan', 
            model_name = "celeritQFTexN2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))



vbQFD <- vb(modelQFDexN2, data = QFD_data,tol_rel_obj = 0.001)
summ_vbQFD <- summary(vbQFD)
plot(summ_vbQFD[[1]][1:299 + 322,1])
(summ_vbQFD[[1]][300 : 322,1])
fitQFD <- sampling(modelQFDexN2, data = QFD_data,control = list(adapt_delta = 0.999, max_treedepth=15), iter = 2000)
summQFD <- summary(fitQFD)
plot(summQFD[[1]][1:299 + 322,1]) 

par(mfrow = c(3,1))
plot(summ_vbQFD[[1]][1:299 + 322,1]) 
plot(simu_signal[-1])
plot(res$state)


fitQFD_array <- as.array(fitQFD)
hist(fitQFD_array[,1,1006])



