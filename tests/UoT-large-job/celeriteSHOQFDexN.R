library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper

decoding <- "viterbi" # or FFBS

# run QFD
rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[10000:17400,c("TIME","PDCSAP_FLUX")]
# handling missing data is currently not implemented, but only little are missing so I will just omit it for now
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])
N <- nrow(rawdata)
#plot(rawdata)

QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                S0_prior = c(-10,10),
                w0_prior = c(-10,10),
                Q_prior = c(-10,10),
                alpha_quiet = c(1,.01), 
                alpha_firing = c(1,1),
                alpha_decay = c(1,.01,.1),
                mu0_quiet = 0,
                lambda_quiet = 10,
                gamma_noise = c(1,0.01),
                mu0_rate_firing = 0,
                sigma_rate_firing = 1e3,
                mu0_rate_decay = 0,
                sigma_rate_decay = 1e3,
                diag = rep(1e-3,N),
                err_prior = c(.01,.01)
                )

modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteSHOQFDexN.stan', 
            model_name = "celeritSHOQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))                       

set.seed(42)
fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=10), chain = 2,iter = 2000, thin = 1,init_r = 15)

summQFD <- summary(fitQFD)

trend_QFD <- summQFD[[1]][1:N+2*N+8, 1]

plot(rawdata)
lines(rawdata[,1], trend_QFD, type = "l", col = "blue")

par(mfrow = c(2,1))
plot(tt[-1], summQFD[[1]][1:(N-1) + (N + 9),1])
plot(tt[-1], rawdata[-1,2])
lines(tt, trend_QFD, type = "l", col = "red")

plot(tt, rawdata[,2]-trend_QFD)

residual <- rawdata[,2]-trend_QFD
plot(tt, residual)
flare3s <- which(residual >= (mean(residual) + 3 * sd(residual)))
points(tt[flare3s], residual[flare3s], col = "red")


QFD_samples <- as.data.frame(fitQFD)
Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 19)]
Viterbi_max <- apply(Viterbi_raw,2,majority)


residual <- rawdata[,2]-trend_QFD

plot(rawdata)
lines(rawdata[,1], summQFD[[1]][1:N+2*N+8, 1], col = "#d400ff",lwd=3.0)
points(rawdata[which(Viterbi_max==2)+1,], col = "red",lwd=3.0)
points(rawdata[which(Viterbi_max==3)+1,], col = "blue",lwd=3.0)
legend("topleft", legend = c("Firing","Decay","Trend"), 
                lty = c(NA,NA,1), pch = c(1,1,NA), col = c("red","blue","#d400ff"),
                cex = 1.5)

save.image("./Res/UoT-large-job/res100-174-cSHOQFDexN.RData")

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celeriteSHO.stan', 
            model_name = "celeritSHO", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))
fitcelerite <- sampling(modelcelerite, data = QFD_data, chains = 2)
