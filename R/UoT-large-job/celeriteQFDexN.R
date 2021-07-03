library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper

decoding <- "viterbi" # or FFBS

# run QFD
rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[10000:17400,c("TIME","PDCSAP_FLUX")]
# handling missing data is currently not implemented, but only little are missing so I will just omit it for now
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])
N <- nrow(rawdata)
plot(rawdata)
tt <- 350 * (rawdata[,1] - min(rawdata[,1]))/(range(rawdata[,1])[2]-range(rawdata[,1])[1]) # sort of normalize the time to avoid super short period

QFD_data <- list(N=N, t = tt,
                y = rawdata[,2],
                sigma_prior = c(-2,5),
                Q0_prior = c(-2,4),
                dQ_prior = c(-2,4),
                period_prior = c(-3,3),
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

if(decoding=="FFBS"){
    modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN_FFBS.stan', 
                model_name = "celeritQFTexNFFBS", 
                allow_undefined = TRUE,
                includes = paste0('\n#include "', 
                                 file.path(getwd(), 
                                 'celerite2/celerite2.hpp'), '"\n'))
}

if(decoding=="viterbi"){
    modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN.stan', 
                model_name = "celeritQFTexN", 
                allow_undefined = TRUE,
                includes = paste0('\n#include "', 
                                 file.path(getwd(), 
                                 'celerite2/celerite2.hpp'), '"\n'))

}
                             

set.seed(42)
fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 8000, thin = 1,init_r = 15)

summQFD <- summary(fitQFD)



QFD_samples <- as.data.frame(fitQFD)

trend_QFD <- summQFD[[1]][1:N + (N + 22),1]
FFBS_raw <- QFD_samples[,1:(N-1) + (3*N + 3 * (N-1) + 22)]
FFBS_max <- apply(FFBS_raw,2,majority)


residual <- rawdata[,2]-trend_QFD

save.image("./Res/UoT-large-job/res100-174-cQFDexN.RData")
