library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper

decoding <- "viterbi" # or FFBS

# run QFD
rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[15400:17400,c("TIME","PDCSAP_FLUX")]
# handling missing data is currently not implemented, but only little are missing so I will just omit it for now
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2] - mean(rawdata[,2])
N <- nrow(rawdata)
#plot(rawdata)
QFD_data <- list(N=N, t = rawdata[,1],
                y = rawdata[,2],
                sigma_prior = c(-8,8),
                Q0_prior = c(-8,8),
                dQ_prior = c(-8,8),
                period_prior = c(-8,8),
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
                err_prior = c(0.01,0.01)
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
fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.99, max_treedepth=15), iter = 4000, thin = 1,init_r = 15, chains = 2)

summQFD <- summary(fitQFD)



QFD_samples <- as.data.frame(fitQFD)

trend_QFD <- summQFD[[1]][1:N + (2*N + 21),1]
FFBS_raw <- QFD_samples[,1:(N-1) + (N + 22)]
FFBS_max <- apply(FFBS_raw,2,majority)

plot(rawdata)
plot(rawdata[,1], trend-QFD)


residual <- rawdata[,2]-trend_QFD

save.image("./Res/UoT-large-job/res100-174-cQFDexN.RData")


modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite.stan', 
                       model_name = "celerit2", 
                       allow_undefined = TRUE,
                       includes = paste0('\n#include "', 
                                         file.path(getwd(), 
                                                   'celerite2/celerite2.hpp'), '"\n'))


fitcelerite <- sampling(modelcelerite, data = QFD_data, chains = 2, init_r = 5)
summ_celerite <- summary(fitcelerite)
plot(summ_celerite[[1]][1:N + (N+23),1])
