library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./R/misc.R") # some helper
source("./R/simuFlares.R")

rawdata <- read.csv("./Data/tess2019006130736-s0007-0000000131799991-0131-s_lc.csv")[10000:11500,c("TIME","PDCSAP_FLUX")]
rawdata <- na.omit(rawdata)
rawdata[,2] <- rawdata[,2]-mean(rawdata[,2])

res_path <- "./Res/Injection_recover/"
file_num <- 1
xm <- 50
alpha <- 1
n_inj <- 5
t_half <- .00005
res_file <- paste0(res_path,"inj_rec_",file_num,".csv")

n_rep <- 5
i_res <- 1

n_res <- n_rep * 2


modelQFD <- stan_model(file = './Stan/Morphology/QFD/CeleriteQFDexN.stan', 
            model_name = "celeritQFTexN", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))

modelcelerite <- stan_model(file = './Stan/Prototypes/Celerite/celerite.stan', 
            model_name = "celerit2", 
            allow_undefined = TRUE,
            includes = paste0('\n#include "', 
                             file.path(getwd(), 
                             'celerite2/celerite2.hpp'), '"\n'))


res_df <- data.frame(matrix(NA, n_res, 13))
colnames(res_df) <- c("n_inj","t_half","xm","alpha","i_rep","method","n_detected","n_injected","TP","FP","FN","SEN","SPC")

set.seed(12345)
for(i in 1:n_rep){
    keplerflare_sim <- kepler_flare(rawdata[,1], t_half, n_inj,rPareto, xm, alpha)
    injected <- rawdata
    injected[,2] <- injected[,2] + keplerflare_sim$flare

    N <- nrow(injected)

    QFD_data <- list(N=N, t = injected[,1],
                y = injected[,2],
                sigma_prior = c(-8,8),
                Q0_prior = c(-8,8),
                dQ_prior = c(-8,8),
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
    fitQFD <- sampling(modelQFD, data = QFD_data,control = list(adapt_delta = 0.9, max_treedepth=15), iter = 3000,init_r = 15, chains = 2)
    QFD_samples <- as.data.frame(fitQFD)
    Viterbi_raw <- QFD_samples[,1:(N-1) + (N + 22)]
    Viterbi_max <- apply(Viterbi_raw,2,majority)
    res_df[i_res,1:5] <- c(n_inj,t_half,xm,alpha,i)
    res_df$method[i_res] <- "QFD"
    res_df[i_res,7:13] <- error_decoding(Viterbi_max, keplerflare_sim$states)
    i_res <- i_res + 1
    write.csv(res_df, res_file)

    celeritedata <- QFD_data
    celeritedata$err_prior <- c(0.01,0.01)

    fitcelerite <- sampling(modelcelerite, data = celeritedata,control = list(adapt_delta = 0.9, max_treedepth=15), iter = 3000,init_r = 2, chains = 2)
    summcelerite <- summary(fitcelerite)
    celerite_trend <- summcelerite[[1]][1:N + (N+23),1]
    residual <- rawdata[,2] - celerite_trend

    flares3sigma <- (residual >= (mean(residual) + 3 * sd(residual))) + 1
    res_df[i_res,1:5] <- c(n_inj,t_half,xm,alpha,i)
    res_df$method[i_res] <- "1-3-sigma"
    res_df[i_res,7:13] <- error_decoding(flares3sigma, keplerflare_sim$states)
    i_res <- i_res + 1
    write.csv(res_df, res_file)
}
