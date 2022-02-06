library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper

# run QFD
rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[16000:17400,c("TIME","PDCSAP_FLUX")]
rawdata <- read.csv("./Data/tess2018206045859-s0001-0000000031381302-0120-s_lc.csv")[6000:8000,c("TIME","PDCSAP_FLUX")]

# handling missing data is currently not implemented, but only little are missing so I will just omit it for now

rawdata[,2] <- rawdata[,2] - mean(rawdata[,2], na.rm = TRUE)
observed <- (!is.na(rawdata[,2])) * 1
rawdata[is.na(rawdata[,2]),2] <- 0
N <- nrow(rawdata)
plot(rawdata)

QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                observed = observed,
                sigma_prior = c(-8,8),
                #Q0_prior = c(2,4),
                Q0_prior = c(-8,8),# this is key, we need to set quality to be not too small
                dQ_prior = c(-8,8),
                #period_prior = c(-3,3),
                period_prior = c(-8,8),
                f_prior = c(1e-6,1-1e-6),
                alpha_quiet = c(1,.1), 
                alpha_firing = c(1,1),
                alpha_decay = c(1,.1,1),
                mu0_quiet = 0,
                lambda_quiet = .01,
                gamma_noise = c(0.01,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(1e-6,N)
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN-missing-handling.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000,init_r = 15, chains = 2)
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

#save.image("../res164-174-cQFDexN.RData")
load("../res164-174-cQFDexN.RData")

## visualize
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

# celerite along

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite-missing-handling.stan', 
            model_name = "celerit2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

celeritedata <- QFD_data
celeritedata$err_prior <- c(0.01,0.01)

fitcelerite <- sampling(modelcelerite, data = celeritedata,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 2000,init_r = 2, chains = 2)
summcelerite <- summary(fitcelerite)
celerite_trend <- summcelerite[[1]][1:N + (N+23),1]
residual <- rawdata[,2] - celerite_trend

flares3sigma <- residual >= (mean(residual) + 3 * sd(residual))

save.image("../res031381302_6000-8000-cQFDexN.RData")
load("../res031381302_6000-8000-cQFDexN.RData")

pdf("./Res/CeleriteQFD/031381302_6000-8000/det_compare2.pdf", width = 7, height = 7.5)
par(mfrow = c(2,1))
plot(rawdata, main = "Proposed HMM", ylab = "Centered Flux", xlab = "Time")

points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+21, 1], col = "#d400ff",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = 1.5)

plot(rawdata, main = "sigma-clipping", ylab = "Centered Flux", xlab = "Time")

points(rawdata[flares3sigma,], col = "red",lwd=3.0)
lines(rawdata[,1], summcelerite[[1]][1:N + (N+23), 1], col = "#d400ff",lwd=3.0)
legend("topleft", legend = c("Potential Flares","Trend"), 
                lty = c(NA,1), pch = c(1,NA), col = c("red","#d400ff"),
                cex = 1.5)
dev.off()

detrended_data <- data.frame(TIME = rawdata$TIME, detrended = rawdata[,2]-summQFD[[1]][1:N+2*N+21, 1])

pdf("./Res/CeleriteQFD/131799991_16400-17400/det_flares.pdf", width = 10, height = 10)
par(mfrow = c(2,1))
plot(detrended_data[90:140,],type = "l")
points(detrended_data[90:140,])
points(detrended_data[110,],col = "red")
points(detrended_data[111:120,], col = "blue")
plot(detrended_data[500:550,],type = "l") 
points(detrended_data[500:550,])
points(detrended_data[c(520,527),],col = "red")
points(detrended_data[c(521:526,528:533),], col = "blue")
dev.off()
