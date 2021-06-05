library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("./R/misc.R")
T <- 500
mu <- c(-3,3)
sigma <- c(1,1)
theta <- matrix(c(0.9,.1,.99,.01),2,2,byrow = T)
p <- c(0.5,0.5)
K <- 2

simu_data <- simuHMM(T, mu, sigma, theta, p)
plot(simu_data)
simu_data[5:7] <- NA
observed <- 1 * (!is.na(simu_data))
simu_data[is.na(simu_data)] <- simu_data[1]

simu_stan <- list(T=T,
                    K=2, 
                    y = simu_data,
                    observed = observed, 
                    alpha = c(1,1),
                    mu0=0,
                    lambda = 1e-6,
                    shape = 0.001,
                    rate = 0.001)

fit <- stan(file = './Stan/Vanilla.stan', data = simu_stan, 
                iter = 3000, control = list(adapt_delta = .999))
plot(fit)
