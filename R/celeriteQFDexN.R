library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper

rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16400:17400,c("TIME","PDCSAP_FLUX")]
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])
N <- nrow(rawdata)
plot(rawdata)
tt <- 50 * (rawdata[,1] - min(rawdata[,1]))/(range(rawdata[,1])[2]-range(rawdata[,1])[1]) # sort of normalize the time to avoid super short period

QFD_data <- list(N=N, t = tt,
                y = rawdata[,2],
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
                diag = rep(1e-6,N)
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 3000,init_r = 5)
summQFD <- summary(fitQFD)


plot(tt, rawdata[,2], col = "blue")
lines(tt, summQFD[[1]][1:N+2*N+21, 1], type = "l")

par(mfrow = c(2,1))
plot(tt[-1], summQFD[[1]][1:(N-1) + (N + 22),1])
plot(tt[-1], rawdata[-1,2])
lines(tt, summQFD[[1]][1:N+2*N+21, 1], type = "l", col = "red")

plot(tt, rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1])

residual <- rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1]
plot(tt, residual)
flare3s <- which(residual >= (mean(residual) + 3 * sd(residual)))
points(tt[flare3s], residual[flare3s], col = "red")

save.image("../res164-174-cQFDexN.RData")
load("../res164-174-cQFDexN.RData")

QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]

Viterbi_max <- apply(Viterbi_raw,2,majority)

pdf("./Res/CeleriteQFD/131799991_16400-17400/det.pdf", width = 10, height = 6)
plot(rawdata)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = 1.5)
dev.off()
