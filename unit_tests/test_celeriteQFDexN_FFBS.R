library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

source("./R/sampling_model.R")
N <- 300
set.seed(12345)
res <- simuQFDexN(N = N, rate_firing = .1, 
                theta_quiet = c(0.95,0.05),
                theta_firing = c(0.15, 0.85),
                theta_decay = c(0.2,0,0.8),rate_decay = 0.6, 
                sigma_noise = 1)

tt <- seq(0,10,length.out = N)
f <- 5*sin(3 * tt)
simu_signal <- f+res$timeseries

QFD_data <- list(N=N, t = tt,
                y = simu_signal,
                sigma_prior = c(-2,5),
                #Q0_prior = c(2,4),
                Q0_prior = c(-2,4),# this is key, we need to set quality to be not too small
                dQ_prior = c(-2,4),
                period_prior = c(-3,3),
                #period_prior = c(0,3),
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
                diag = rep(1e-6,N),
                err_prior = c(0.01, 0.01) # just for celerite, same as gamma_noise
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN_FFBS.stan', 
            model_name = "celeritQFTexNFFBS", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))


fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000)
vbQFD <- vb(modelQFD, data = QFD_data)

QFD_vb <- as.data.frame(vbQFD)
trend_QFD_vb <- summary(vbQFD)[[1]][1:N + (N + 22),1]
FFBS_raw_vb <- QFD_vb[,1:(N-1) + (3*N + 3 * (N-1) + 22)]
FFBS_max_vb <- apply(FFBS_raw_vb,2,majority)


plot(simu_signal, main = "QFD")
lines(trend_QFD_vb, col = "#d400ff",lwd=3.0)
points(which(FFBS_max_vb==2)+1,simu_signal[which(FFBS_max_vb==2)+1], col = "red",lwd=3.0)
points(which(FFBS_max_vb==3)+1,simu_signal[which(FFBS_max_vb==3)+1], col = "blue",lwd=3.0)



QFD_samples <- as.data.frame(fitQFD)
trend_QFD <- summary(fitQFD)[[1]][1:N + (N + 22),1]
FFBS_raw <- QFD_samples[,1:(N-1) + (3*N + 3 * (N-1) + 22)]
FFBS_max <- apply(FFBS_raw,2,majority)

plot(simu_signal, main = "QFD")
lines(trend_QFD, col = "#d400ff",lwd=3.0)
points(which(FFBS_max==2)+1,simu_signal[which(FFBS_max==2)+1], col = "red",lwd=3.0)
points(which(FFBS_max==3)+1,simu_signal[which(FFBS_max==3)+1], col = "blue",lwd=3.0)

